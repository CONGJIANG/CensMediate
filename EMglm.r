library(data.table) # Assuming 'data' is a data.table
n_obs <- 500
(data <- Study_dgp_Mcensor(n_obs = n_obs, lod = 1.0))
sum(data$study_dat$censored)/n_obs *100

head(data$study_dat)
head(data$study_dat_full)

dataset <- data$study_dat
dataset



sample_truncated_proposal <- function(n, data, LOD, beta_m) {
  meanL <- as.matrix(data[c("L1", "L2", "L3")]) %*% beta_m[c("L1", "L2", "L3")]
  rtruncnorm(n, a = 0, b = LOD, mean = beta_m["A"] * data$A + meanL, sd = beta_m["sd"])
}

# Function to compute f(M)
compute_f_M <- function(data, LOD, beta_m) {
  # Calculate the means for each data point
  meanL <- as.matrix(data[, .(L1, L2, L3)]) %*% beta_m[c("L1", "L2", "L3")]
  means <- beta_m["A"] * data$A + meanL
  # Compute the probability density for each M
  f_M <- dtruncnorm(data$M, a = 0, b = LOD, mean = means, sd = beta_m["sd"])
  return(f_M)
}





########################################
beta_m <- c(A = 0.5, L1 = 0.3, L2 = 0.2, L3 = 0.1, sd = 1)
LOD_value <- 1.0  # Set the LOD value
library(truncnorm)

# Define the E-step
Impute_step <- function(data, beta_m, S, LOD = LOD_value) {
  # Duplicate and impute LOD rows
  imputed_rows_list <- lapply(1:nrow(data[data$M == "LOD", ]), function(i) {
    replicated_row <- data[data$M == "LOD", ][rep(i, S), , drop = FALSE]
    replicated_row$M <- sample_truncated_proposal(
      n = S,
      data = as.data.frame(replicated_row),
      LOD = LOD,
      beta_m = beta_m
    )
    return(replicated_row)
  })
  
  # Combine non-LOD data with imputed LOD rows
  data <- rbind(data[data$M != "LOD", ], do.call(rbind, imputed_rows_list))
  # Convert all character columns to numeric or factor as needed
  data <- data %>%
    mutate(across(where(is.character), as.numeric))
  # Ensure the outcome and weights are numeric
  data$censored <- as.numeric(data$censored)
  return(data)
}

(data_imputed <- Impute_step(data = dataset, beta_m, S = 10))


mod_update_glm <- function(data) {
  # Check and assign w_norm if NULL
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  # Censoring Modeling (non parametric assumption, include Y? )
  cen_mod <- glm(censored ~ L1 + L2 + L3 + A + M, data = data, family = binomial, weights = w_norm)
  cen_prob <- predict(cen_mod, type = "response")
  # Outcome Modeling
  out_mod <- glm(Y ~ L1 + L2 + L3 + A + M, data = data, family = binomial, weights = w_norm)
  out_prob <- predict(out_mod, type = "response")
  # Return predictions and coefficients
  return(list(
    cen_prob = cen_prob, 
    out_prob = out_prob,
    cen_coef = coef(cen_mod),
    out_coef = coef(out_mod)
  ))
}

mod_pred <- mod_update_glm(data_imputed)
mod_pred$cen_coef

# mod_pred are from mod_update_glm funciton, which inlucdes both out_prob and cen_prob
Wet_step <- function(data_imputed, mod_pred, beta_m_new, beta_m_0, LOD) {
  # Update probabilities
  data_imputed$cen_prob <- mod_pred$cen_prob
  data_imputed$out_prob <- mod_pred$out_prob
  
  # Filter censored data
  censored_data <- data_imputed[data_imputed$censored == TRUE, ]
  
  # Compute new and initial weights for censored data
  cenwet_new <- compute_f_M(censored_data, LOD = LOD, beta_m_new)
  cenwet_int <- compute_f_M(censored_data, LOD = LOD, beta_m_0)
  
  # Handle non-censored data weights
  noncen <- rep(1, nrow(data_imputed) - nrow(censored_data))
  
  # Update weights (ensure correct reference)
  data_imputed$w <- data_imputed$cen_prob * data_imputed$out_prob * c(noncen, cenwet_new / cenwet_int)
  
  # Normalize weights by ID
  data_imputed$w_norm <- ave(data_imputed$w, data_imputed$id, FUN = function(x) x / sum(x))
  
  # Remove temporary columns
  data_weted <- data_imputed[, !names(data_imputed) %in% c("cen_prob", "out_prob"), with = FALSE]
  
  return(data_weted)
}

# Apply the updated function
data_weted <- Wet_step(data_imputed, mod_pred, beta_m_new = beta_m, beta_m_0 = beta_m, LOD = LOD_value)







# Define the M-step
M_stepf_M <- function(data, LOD = LOD_value, beta_m) {
  # Define an internal function for optimization
  neg_log_likelihood <- function(beta_m, data, LOD) {
    # Extract weights
    weights <- data$w_norm
    # Calculate the mean for each sample using provided beta_m
    meanL <- as.matrix(data[, .(L1, L2, L3)]) %*% beta_m[c("L1", "L2", "L3")]
    means <- beta_m["A"] * data$A + meanL
    
    # Compute log-likelihood for each data point using the truncated normal distribution
    # Log transformation is done using the log() function
    log_likelihoods <- log(dtruncnorm(data$M, a = 0, b = LOD, mean = means, sd = beta_m["sd"]))
    
    # Incorporate weights into the likelihood calculation
    weighted_log_likelihood <- weights * log_likelihoods
    
    # Return the negative of the weighted log-likelihood sum
    return(-sum(weighted_log_likelihood))
  }
  
  # Subset the data where censored is TRUE
  censored_data <- data[data$censored == TRUE, ]
  optim_results <- optim(par = beta_m, fn = neg_log_likelihood, data = censored_data, LOD = LOD_value)
  
  # Results
  return(beta_m_new = optim_results$par)
}

M_stepf_M(data = data_weted, LOD = LOD_value, beta_m)








# EM Algorithm
EM_algorithm <- function(data, beta_m_0, S = 20, LOD = LOD_value, max_iter = 1000, tol = 1e-16) {
  iter <- 0
  diff <- Inf
  # I-Step
  data_imputed <- Impute_step(data, beta_m_0, S)
  mod_pred <- mod_update_glm(data_imputed)
  beta_m <- beta_m_0
  
  while (iter < max_iter && diff > tol) {
    # W-Step
    data_weted <- Wet_step(data_imputed, mod_pred, beta_m_new = beta_m, beta_m_0 = beta_m, LOD = LOD_value)
    
    # Maximization Step
    beta_m_new <- M_stepf_M(data = data_weted, LOD = LOD_value, beta_m)
    
    # Calculate convergence criteria
    diff <- sum((beta_m_new - beta_m)^2)
    cen_coef_old <- mod_pred$cen_coef
    out_coef_old <- mod_pred$out_coef

    # parameter updating
    beta_m <- beta_m_new
    mod_pred <- mod_update_glm(data_weted)
    
    diff <- diff + sum((mod_pred$cen_coef - cen_coef_old)^2) + sum((mod_pred$out_coef - out_coef_old)^2)
    iter <- iter + 1
    cat("Iteration:", iter, "Difference:", diff, "\n")
  }
  
  return(list(beta_m = beta_m, data_weted = data_weted))
}


# Define initial guesses for beta_m parameters
beta_m <- c(A = 0.5, L1 = 0.3, L2 = 0.2, L3 = 0.1, sd = 1)
LOD_value <- 1.0  # Set the LOD value
# Run the EM algorithm
dat_aug <- EM_algorithm(dataset, beta_m)

dat_aug <- dat_aug$data_weted
dim(dat_aug)
head(dat_aug)


gcomp <- function(data) {
  # Extract variables from the dataset
  y <- data$Y        # Assuming Y is the outcome variable
  m <- data$M        # Assuming M is the mediator
  a <- data$A        # Assuming A is the treatment
  l1 <- data$L1      # Assuming L1 is the first covariate
  l2 <- data$L2      # Assuming L2 is the second covariate
  l3 <- data$L3      # Assuming L3 is the third covariate
  weights <- data$w_norm  # Extract the weights
  
  # Combine L1, L2, L3 into a data frame for modeling convenience
  covariates <- data.frame(l1 = l1, l2 = l2, l3 = l3)
  
  # Fit a weighted linear model for y including m, a, and covariates
  lm_y <- lm(y ~ m + a + l1 + l2 + l3, weights = weights)
  
  # Predict potential outcomes under different treatments
  pred_y1 <- predict(lm_y, newdata = transform(covariates, a = 1, m = m))
  pred_y0 <- predict(lm_y, newdata = transform(covariates, a = 0, m = m))
  
  # Calculate the pseudo outcome
  pseudo <- pred_y1 - pred_y0
  
  # Fit a weighted linear model using the pseudo outcome
  lm_pseudo <- lm(pseudo ~ a + l1 + l2 + l3, weights = weights)
  
  # Predict the causal effect when a = 0
  pred_pseudo <- predict(lm_pseudo, newdata = transform(covariates, a = 0))
  
  # Calculate and return the estimate
  estimate <- mean(pred_pseudo)
  return(estimate)
}


estimate <- gcomp(dat_aug)
print(estimate)

dat_aug


# instantiate learners
mean_lrnr <- Lrnr_mean$new()
fglm_lrnr <- Lrnr_glm_fast$new()
lasso_lrnr <- Lrnr_glmnet$new(alpha = 1, nfolds = 3)
rf_lrnr <- Lrnr_ranger$new(num.trees = 200)

# create learner library and instantiate super learner ensemble
lrnr_lib <- Stack$new(mean_lrnr, fglm_lrnr, lasso_lrnr, rf_lrnr)
sl_lrnr <- Lrnr_sl$new(learners = lrnr_lib, metalearner = Lrnr_nnls$new())


# compute one-step estimate of the natural direct effect
nde_onestep <- medoutcon(
  W = dat_aug[, c("L1", "L2", "L3")],
  A = as.numeric(dat_aug$A),
  Z = NULL,
  M = dat_aug$M,
  Y = dat_aug$Y,
  obs_weights = dat_aug$w_norm,
  g_learners = lasso_lrnr,
  h_learners = lasso_lrnr,
  b_learners = lasso_lrnr,
  effect = "direct",
  estimator = "onestep",
  estimator_args = list(cv_folds = 5)
)
summary(nde_onestep)


# compute one-step estimate of the natural indirect effect
nie_onestep <- medoutcon(
  W = dat_aug[, c("L1", "L2", "L3")],
  A = as.numeric(dat_aug$A),
  Z = NULL,
  M = dat_aug$M,
  Y = dat_aug$Y,
  obs_weights = dat_aug$w_norm,
  g_learners = lasso_lrnr,
  h_learners = lasso_lrnr,
  b_learners = lasso_lrnr,
  effect = "indirect",
  estimator = "onestep",
  estimator_args = list(cv_folds = 5)
)
summary(nie_onestep)
