remotes::install_github("nhejazi/medoutcon")
devtools::install_github("tlverse/tlverse")
library(data.table) # Assuming 'data' is a data.table
library(dplyr)
library(sl3)
library(hal9001)
library(medoutcon)
library(truncnorm)
library(xgboost)
n_obs <- 5000

(data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = 0.3))
(LOD <- data$lod)
(LOD_value <- as.numeric(data$lod))

head(data$study_dat)
head(data$study_dat_full)

dataset <- data$study_dat
dataset


NaiveImpute <- function(data, LOD = LOD_value) {
  # Identify rows where M == "LOD"
  lod_rows <- data[data$M == "LOD", , drop = FALSE]
  
  # If there are LOD rows, replace M with LOD_value / 2
  if (nrow(lod_rows) > 0) {
    lod_rows$M <- as.numeric(LOD) / 2
  }
  
  # Combine original non-LOD data with imputed LOD rows
  data <- rbind(data[data$M != "LOD", ], lod_rows)
  # Ensure M is numeric
  data$M <- as.numeric(data$M)
  # Ensure outcome and weights are numeric
  data$censored <- as.numeric(data$censored)
  data$w_norm <- rep(1, nrow(data))
  data <- data %>% arrange(id)
  return(data)
}


NaiveImpute(dataset, LOD = LOD_value)

###############################
## Density of M below LOD
##
sample_truncated_proposal <- function(n, data, LOD, beta_m) {
  meanL <- as.matrix(data[c("L1", "L2", "L3")]) %*% beta_m[c("L1", "L2", "L3")]
  rtruncnorm(n, a = 0, b = LOD, mean = beta_m["A"] * data$A + meanL, sd = beta_m["sd"])
}

# Density Function to compute f(M)
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
LOD_value <- as.numeric(data$lod)  # Set the LOD value

########################################
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
  # Censoring Modeling (parametric model, does not include Y)
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
data_weted





########################################
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
EM_algorithm <- function(data, beta_m_0, S = 50, LOD = LOD_value, max_iter = 1000, tol = 1e-16) {
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
LOD_value <- as.numeric(data$lod)  # Set the LOD value
# Run the EM algorithm
dat_res <- EM_algorithm(dataset, beta_m)

dat_aug <- dat_res$data_weted
dat_res$beta_m
dim(dat_aug)
head(dat_aug)


gcomp.nde <- function(data) {
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
  lm_y <- glm(y ~ m + a + l1 + l2 + l3, weights = weights, family = binomial)
  
  # Predict potential outcomes under different treatments
  pred_y1 <- predict(lm_y, newdata = transform(covariates, a = 1, m = m), type = "response")
  pred_y0 <- predict(lm_y, newdata = transform(covariates, a = 0, m = m), type = "response")
  
  # Fit a weighted linear model using the pseudo outcome
  lm_y1 <- glm(pred_y1 ~ a + l1 + l2 + l3, weights = weights, family = binomial)
  lm_y0 <- glm(pred_y0 ~ a + l1 + l2 + l3, weights = weights, family = binomial)
  
  # Predict the causal effect when a = 0
  y1 <- predict(lm_y1, newdata = transform(covariates, a = 0), type = "response")
  y0 <- predict(lm_y0, newdata = transform(covariates, a = 0), type = "response")
  
  # Calculate and return the estimate
  nde.rd <- mean(y1 - y0); nde.rr <- mean(y1)/mean(y0); nde.or <- (mean(y1)/(1-mean(y1))) / (mean(y0)/(1-mean(y0)))
  return(list(nde.rd = nde.rd, nde.rr = nde.rr, nde.or = nde.or))
}


nde_gcom_glm <- gcomp.nde(dat_aug)
nde_gcom_glm$nde.rd

# Perform bootstrap
gcom_boots_ci <- function(data, n_bootstrap = 10, alpha = 0.05, type = c("NDE", "NIE")) {
  n <- nrow(data)
  boots_est <- numeric(n_bootstrap)
  for (i in seq_len(n_bootstrap)) {
    # Resample data with replacement
    bootstrap_sample <- data[sample(1:n, size = n, replace = TRUE), ]
    # Calculate the point estimate for the bootstrap sample
    if (type == "NDE") {
      bres <- gcomp.nde(bootstrap_sample)
      boots_est[i] <- bres$nde.rd  # Extract NDE estimate
    } else if (type == "NIE") {
      bres <- gcomp.nie(bootstrap_sample)
      boots_est[i] <- bres$nie.rd  # Extract NIE estimate
    }
  }
  
  # Calculate the 95% CI using percentiles
  ci_lower <- quantile(boots_est, probs = alpha / 2)
  ci_upper <- quantile(boots_est, probs = 1 - alpha / 2)
  
  list(
    point_estimate = mean(boots_est),  # Mean of bootstrap estimates
    ci = c(ci_lower, ci_upper)  # Confidence interval
  )
}

# Example usage with your dataset
(res1 <- gcom_boots_ci(dat_aug, type = "NDE"))

gcomp.nie <- function(data) {
  # Extract variables from the dataset
  y <- data$Y        # Outcome variable
  m <- data$M        # Mediator
  a <- data$A        # Treatment
  l1 <- data$L1      # Covariate 1
  l2 <- data$L2      # Covariate 2
  l3 <- data$L3      # Covariate 3
  weights <- data$w_norm  # Weights for weighted regression
  
  # Combine covariates into a data frame for modeling convenience
  covariates <- data.frame(l1 = l1, l2 = l2, l3 = l3)
  
  # Fit a weighted linear model for y including m, a, and covariates
  lm_y <- glm(y ~ m + a + l1 + l2 + l3, weights = weights, family = binomial)
  
  # Predict potential outcomes under different treatments
  pred_y1 <- predict(lm_y, newdata = transform(covariates, a = 1, m = m),  type = "response")
  # Fit a weighted linear model using the predict potential outcomes
  lm_y1 <- glm(pred_y1 ~ a + l1 + l2 + l3, weights = weights,  family = binomial)
  
  # Predict the causal effect when a = 0
  y1 <- predict(lm_y1, newdata = transform(covariates, a = 1),  type = "response")
  y0 <- predict(lm_y1, newdata = transform(covariates, a = 0),  type = "response")
  
  # Calculate and return the estimate
  nie.rd <- mean(y1 - y0); nie.rr <- mean(y1)/mean(y0); nie.or <- (mean(y1)/(1-mean(y1))) / (mean(y0)/(1-mean(y0)))
  # Return the estimate of the indirect effect
  return(list(nie.rd = nie.rd, nie.rr = nie.rr, nie.or = nie.or))
}
nie_gcom_glm <- gcomp.nie(dat_aug)
nie_gcom_glm
(res2 <- gcom_boots_ci(dat_aug, type = "NIE")$ci)

est_med_effs <- function(dat_aug) {
  # Load necessary libraries
  # library(sl3)
  # library(tlverse)
  
  # Instantiate learners
  fglm_lrnr <- Lrnr_glm_fast$new()
  lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
  lasso_lrnr <- Lrnr_glmnet$new(alpha = 1, nfolds = 3)
  lrn_earth <- Lrnr_earth$new()
  rf_lrnr <- Lrnr_ranger$new(num.trees = 200)
  # fast highly adaptive lasso (HAL) implementation
  lrn_hal <- Lrnr_hal9001$new(max_degree = 2, num_knots = c(3,2), nfolds = 5)
  lrn_xgb <- Lrnr_xgboost$new()
  
  # Create learner library and instantiate super learner ensemble
  lrnr_lib <- Stack$new(fglm_lrnr, lrn_ridge, lasso_lrnr, lrn_xgb)
  sl_lrnr <- Lrnr_sl$new(learners = lrnr_lib, metalearner = Lrnr_nnls$new())
  
  # Compute one-step and TMLE estimates for NDE and NIE
  nde_onestep_glm <- medoutcon(
    W = dat_aug[, c("L1", "L2", "L3")],
    A = as.numeric(dat_aug$A),
    Z = NULL,
    M = dat_aug$M,
    Y = dat_aug$Y,
    obs_weights = dat_aug$w_norm,
    g_learners = sl_lrnr, # PS  # interaction here lasso of PS is not good
    h_learners = sl_lrnr, # modefied PS with Mediator.  hard to model P(A | M, L)
    b_learners = lasso_lrnr, # Outcome regresssion
    effect = "direct",
    estimator = "onestep",
    estimator_args = list(cv_folds = 3)
  )
  
  nde_tmle_glm <- medoutcon(
    W = dat_aug[, c("L1", "L2", "L3")],
    A = as.numeric(dat_aug$A),
    Z = NULL,
    M = dat_aug$M,
    Y = dat_aug$Y,
    obs_weights = dat_aug$w_norm,
    g_learners = sl_lrnr,
    h_learners = sl_lrnr,
    b_learners = lasso_lrnr,
    effect = "direct",
    estimator = "tmle",
    estimator_args = list(cv_folds = 3)
  )
  
  nie_onestep_glm <- medoutcon(
    W = dat_aug[, c("L1", "L2", "L3")],
    A = as.numeric(dat_aug$A),
    Z = NULL,
    M = dat_aug$M,
    Y = dat_aug$Y,
    obs_weights = dat_aug$w_norm,
    g_learners = sl_lrnr,
    h_learners = sl_lrnr,
    b_learners = lasso_lrnr,
    effect = "indirect",
    estimator = "onestep",
    estimator_args = list(cv_folds = 3)
  )
  
  nie_tmle_glm <- medoutcon(
    W = dat_aug[, c("L1", "L2", "L3")],
    A = as.numeric(dat_aug$A),
    Z = NULL,
    M = dat_aug$M,
    Y = dat_aug$Y,
    obs_weights = dat_aug$w_norm,
    g_learners = sl_lrnr,
    h_learners = sl_lrnr,
    b_learners = lasso_lrnr,
    effect = "indirect",
    estimator = "tmle",
    estimator_args = list(cv_folds = 3)
  )
  
  nde_onestep <- summary(nde_onestep_glm)
  nde_tmle <- summary(nde_tmle_glm)
  nie_onestep <- summary(nie_onestep_glm)
  nie_tmle <- summary(nie_tmle_glm)
  
  # Return the summaries as a named list
  return(list(
    nde_onestep_est = nde_onestep$param_est,
    nde_onestep_ci = c(nde_onestep$lwr_ci, nde_onestep$upr_ci),
    nde_tmle_est = nde_tmle$param_est,
    nde_tmle_ci = c(nde_tmle$lwr_ci, nde_tmle$upr_ci),
    nie_onestep_est = nie_onestep$param_est,
    nie_onestep_ci = c(nie_onestep$lwr_ci, nie_onestep$upr_ci),
    nie_tmle_est = nie_tmle$param_est,
    nie_tmle_ci = c(nie_tmle$lwr_ci, nie_tmle$upr_ci)
  ))
}

(medeff <- est_med_effs(dat_aug))
str(medeff)
str(medeff$nde_onestep_est)



MC_sim <- function(r, n_obs, type = c("Trad", "EM-MLE"), censor_rate = 0.3, beta_m = c(A = 0.5, L1 = 0.3, L2 = 0.2, L3 = 0.1, sd = 1)) {
  type <- match.arg(type)
  output_file <- paste0("MED_FEB04_", type, "_n_", n_obs, ".txt")
  
  for (i in 1:r) {
    data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = censor_rate)
    simdata<- data$study_dat
    LOD_value <- as.numeric(data$lod)
    
    if (type == "Trad") {
      data_aug <- NaiveImpute(simdata, LOD_value)
    } else if (type == "EM-MLE") {
      em_result <- EM_algorithm(simdata, beta_m)
      data_aug <- em_result$data_weted
      beta_m <- em_result$beta_m  # Persisting update
    }
    
    # G-computation
    nde_gcom_glm <- gcomp.nde(data_aug)
    CI.nde <- gcom_boots_ci(data_aug, type = "NDE")$ci
    nie_gcom_glm <- gcomp.nie(data_aug)
    CI.nie <- gcom_boots_ci(data_aug, type = "NIE")$ci
    
    # Custom mediation estimates
    medeff <- est_med_effs(data_aug)
    
    # Save results
    results <- c(i, 
                 nde_gcom_glm$nde.rd, CI.nde, 
                 nie_gcom_glm$nie.rd, CI.nie, 
                 medeff$nde_onestep_est, medeff$nde_onestep_ci,
                 medeff$nde_tmle_est, medeff$nde_tmle_ci, 
                 medeff$nie_onestep_est, medeff$nie_onestep_ci, 
                 medeff$nie_tmle_est, medeff$nie_tmle_ci)
    
    write(results, file = output_file, ncolumns = length(results), append = TRUE)
  }
  
  return(beta_m)  # Returning the final beta_m if needed
}


# Example usage
r <- 1000   # Number of replicates
n_obs <- 1000  # Number of observations per replicate

sim.res <- MC_sim(r = r, type = "EM-MLE", n_obs = n_obs)
sim.res <- MC_sim(r = r, type = "Trad", n_obs = n_obs)










MC_sim <- function(r, n_obs, type = c("Trad", "EM-MLE"), censor_rate = 0.3, beta_m = c(A = 0.5, L1 = 0.3, L2 = 0.2, L3 = 0.1, sd = 1)) {
  type <- match.arg(type)
  output_file <- paste0("MED_FEB04_", type, "_n_", n_obs, ".txt")
  
  for (i in 1:r) {
    data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = censor_rate)
    simdata<- data$study_dat
    LOD_value <- as.numeric(data$lod)
    
    if (type == "Trad") {
      data_aug <- NaiveImpute(simdata, LOD_value)
    } else if (type == "EM-MLE") {
      em_result <- EM_algorithm(simdata, beta_m)
      data_aug <- em_result$data_weted
      beta_m <- em_result$beta_m  # Persisting update
    }
    
    # G-computation
    nde_gcom_glm <- gcomp.nde(data_aug)
    CI.nde <- gcom_boots_ci(data_aug, type = "NDE")$ci
    nie_gcom_glm <- gcomp.nie(data_aug)
    CI.nie <- gcom_boots_ci(data_aug, type = "NIE")$ci
    
    # Save results
    results <- c(i, 
                 nde_gcom_glm$nde.rd, CI.nde, 
                 nie_gcom_glm$nie.rd, CI.nie)
    
    write(results, file = output_file, ncolumns = length(results), append = TRUE)
  }
  
  return(beta_m)  # Returning the final beta_m if needed
}




