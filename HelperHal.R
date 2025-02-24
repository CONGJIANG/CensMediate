remotes::install_github("nhejazi/medoutcon")
devtools::install_github("tlverse/tlverse")
library(data.table) # Assuming 'data' is a data.table
library(dplyr)
library(sl3)
library(hal9001)
library(medoutcon)
library(truncnorm)
library(xgboost)
library(logistf)
library(zoo) 
n_obs <- 300

(data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = 0.1))
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
  data$C <- as.numeric(data$C)
  data$w_norm <- rep(1, nrow(data))
  data <- data %>% arrange(id)
  return(data)
}


NaiveImpute(dataset, LOD = LOD_value)

###############################
## Density of M for samples below LOD
##
sample_density_proposal <- function(n, data, LOD, beta_m) {
  # Linear predictor for L variables
  meanL <- as.matrix(data[c("L1", "L2", "L3")]) %*% beta_m[c("L1", "L2", "L3")]
  # Mean calculation for the normal distribution
  means <- beta_m["A"] * data$A + meanL + beta_m["b0"]
  # Initialize an empty vector to store valid samples
  valid_samples <- numeric(0)
  
  # Log-normal sampling loop until we have enough valid samples
  while (length(valid_samples) < n) {
    # Draw more samples to fill the request
    remaining <- n - length(valid_samples)
    # Draw samples from a log-normal distribution
    samples <- rlnorm(remaining, meanlog = means, sdlog = beta_m["sd"])
    # Append valid samples (less than LOD)
    valid_samples <- c(valid_samples, samples[samples <= LOD])
  }
  
  # Return only the first n valid samples
  return(valid_samples[1:n])
}
# Density Function to compute f(M)
compute_f_M <- function(data, beta_m) {
  # Calculate linear predictors for the mean of the normal distribution
  meanL <- as.matrix(data[, .(L1, L2, L3)]) %*% beta_m[c("L1", "L2", "L3")]
  # Total predicted mean
  means <- beta_m["A"] * data$A + meanL + beta_m["b0"]
  
  # Compute the log-normal density for each M
  # Convert the mean to a log scale for log-normal parameters
  f_M <- dlnorm(data$M, meanlog = means, sdlog = beta_m["sd"])
  return(f_M)
}
########################################
# Define the E-step
Impute_step <- function(data, beta_m, S, LOD) {
  # Duplicate and impute LOD rows
  imputed_rows_list <- lapply(1:nrow(data[data$M == "LOD", ]), function(i) {
    replicated_row <- data[data$M == "LOD", ][rep(i, S), , drop = FALSE]
    replicated_row$M <- sample_density_proposal(
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
  data$C <- as.integer(data$C)
  return(data[order(data$id), ])
}

beta_m <- c(b0 = -3.75, A = 1.75, L1 = 2 , L2 = 1.25, L3 = -0.35, sd = 0.3)
LOD_value <- as.numeric(data$lod)  # Set the LOD value
(data_impu <- Impute_step(data = dataset, beta_m, S = 2, LOD = LOD_value))

#############
# METHOD 2
mod_update_hal <- function(data, basis_list = NULL) {
  # Check and assign w_norm if NULL (for the initial value of the weights)
  if (is.null(data$w_norm)) {
    data$w_norm <- 1
  }
  
  x_names <- c("L1", "L2", "L3", "A", "M")
  y_out_names <- "Y"
  weights_names <- "w_norm"
  
  # fit enumerate_basis on full data to get set of basis functions
  if (is.null(basis_list)) {
    basis_list <- hal9001::enumerate_basis(
      as.matrix(data[, ..x_names]),
      max_degree = 2, smoothness_orders = 1
    )
  }
  # Function to fit HAL model
  out_mod <- hal9001::fit_hal(
    X = as.matrix(data[, ..x_names]),
    Y = as.numeric(data[, get(y_out_names)]),
    weights = as.numeric(data[, get(weights_names)]),
    family = "binomial",
    basis_list = basis_list
  )
  
  out_probY1 <- predict(out_mod, new_data = as.matrix(data[, ..x_names]), type = "response")
  out_probY0 <- 1 - out_probY1
  
  # Create a new column for predicted probabilities based on observed Y values
  out_prob <- ifelse(data$Y == 1, out_probY1, out_probY0)
  
  # Return predictions, coefficients, and terms
  return(list(
    out_prob = out_prob,
    basis_list = basis_list
  ))
}

system.time({
  cv_halmod_pred <- mod_update_hal(data = data_impu)
})
# Access results
cv_halmod_pred$out_prob # Predicted probabilities for the outcome model
cv_halmod_pred$basis_list # Shared basis functions

# mod_pred are from mod_update_glm funciton, which inlucdes both out_prob and cen_prob
Wet_step <- function(data_imputed, mod_pred, beta_m_new, beta_m_0) {
  # Update probabilities
  data_imputed$cen_prob <- pmin(pmax(mod_pred$cen_prob, 1e-40), 1 - 1e-40)
  data_imputed$out_prob <- pmin(pmax(mod_pred$out_prob, 1e-40), 1 - 1e-40) 
  
  # Filter censored data
  censored_data <- data_imputed[data_imputed$C == 0, ]
  
  # Compute new and initial weights for censored data
  cenwet_new <- compute_f_M(censored_data, beta_m_new)
  cenwet_int <- compute_f_M(censored_data, beta_m_0)
  
  censored_data$w <- ((1- censored_data$cen_prob) * censored_data$out_prob) * (cenwet_new / cenwet_int)
  
  # Normalize weights by ID
  censored_data$w_norm <- ave(censored_data$w, censored_data$id, FUN = function(x) x / sum(x))
  
  # For uncensored data
  uncensored_data <- data_imputed[data_imputed$C == 1, ]
  uncensored_data$w_norm <- 1;   uncensored_data$w <- 1
  data_imputed <- rbind(censored_data, uncensored_data)
  
  # Remove temporary columns
  data_weted <- data_imputed[, !names(data_imputed) %in% c("cen_prob", "out_prob"), with = FALSE]
  return(data_weted[order(data_weted$id),])
}
(data_weted <- Wet_step(data_impu, cv_halmod_pred, beta_m_new = beta_m, beta_m_0 = beta_m))


########################################
# Define the M-step, here is for M density
M_stepf_M <- function(data, beta_m) {
  # Define an internal function for optimization
  neg_log_likelihood <- function(beta_m, data) {
    # Extract weights
    weights <- data$w_norm
    # Calculate the mean for each sample using provided beta_m
    meanL <- as.matrix(data[, .(L1, L2, L3)]) %*% beta_m[c("L1", "L2", "L3")]
    means <- beta_m["A"] * data$A + meanL + beta_m["b0"]
    
    # Compute log-likelihood for each data point using the truncated normal distribution
    # Log transformation is done using the log() function
    log_likelihoods <- log(dlnorm(data$M, meanlog = means, sdlog = beta_m["sd"]))
    
    # Incorporate weights into the likelihood calculation
    weighted_log_likelihood <- weights * log_likelihoods
    
    # Return the negative of the weighted log-likelihood sum
    return(-sum(weighted_log_likelihood))
  }
  
  optim_results <- optim(par = beta_m, fn = neg_log_likelihood, data = data)
  
  # Results
  return(beta_m_new = optim_results$par)
}

# Example usage:
# QM <- -3 + 1.5 * A + 1.75 * L1 + 1.5 * L2 - 0.25 * L3,       0.25
beta_m_start <- c(b0 = -2.75, A = 1.75, L1 = 2, L2 = 1.25, L3 = -0.35, sd = 0.5)
M_stepf_M(data = data_weted, beta_m = beta_m_start)


# EM Algorithm
EM_algorithm <- function(data, beta_m_0, S = 2, LOD = LOD_value, max_iter = 50, tol = 1e-10, verbose = TRUE) {
  iter <- 0; diff <- Inf; rel_diff <- Inf
  
  # I-Step: Initial Imputation
  data_imputed <- Impute_step(data, beta_m_0, S, LOD)
  folds <- origami::make_folds(n = nrow(data_imputed), V = 5)
  cv_halmod_pred <- mod_update_halcv(data_imputed, folds)
  halbasis <- cv_halmod_pred$basis_list
  
  beta_m <- beta_m_0
  start_time <- Sys.time()  # Track total runtime
  
  while (iter < max_iter && diff > tol && rel_diff > tol) {
    iter_time <- Sys.time()
    
    # W-Step
    data_weted <- Wet_step(data_imputed, cv_halmod_pred, beta_m_new = beta_m, beta_m_0 = beta_m_0)
    
    # M-Step
    beta_m_new <- M_stepf_M(data = data_weted, beta_m)
    
    # Compute Convergence Criteria
    diff <- sum((beta_m_new - beta_m)^2)
    
    cen_coef_old <- cv_halmod_pred$cen_coef
    out_coef_old <- cv_halmod_pred$out_coef
    
    # Parameter Update
    beta_m <- beta_m_new
    folds <- origami::make_folds(n = nrow(data_weted), V = 5)
    cv_halmod_pred <- mod_update_halcv(data_weted, folds, basis_list = halbasis)
    
    # Update diff with coefficient changes
    #diff <- diff + sum((mod_pred$cen_coef - cen_coef_old)^2) + sum((mod_pred$out_coef - out_coef_old)^2)
    
    # Fixed denominator issue in rel_diff calculation
    #rel_diff <- diff / (sum(beta_m^2) + sum(cen_coef_old^2) + sum(out_coef_old^2) + 1e-8)  # Prevent division by zero
    rel_diff <- diff / (sum(beta_m^2) + 1e-8)  # Prevent division by zero
    iter_time_end <- Sys.time()
    iter_duration <- round(difftime(iter_time_end, iter_time, units = "secs"), 4)
    
    if (verbose) {
      cat("Iteration:", iter, 
          "| Absolute Diff:", diff, 
          "| Relative Diff:", rel_diff,
          "| Time:", iter_duration, "s\n")
    }
    
    iter <- iter + 1
  }
  
  total_runtime <- round(difftime(Sys.time(), start_time, units = "secs"), 4)
  
  # Stopping condition logging
  if (iter >= max_iter) {
    message("EM Algorithm stopped: Maximum iterations reached (", max_iter, ").")
  } else if (diff <= tol || rel_diff <= tol) {
    message("EM Algorithm converged successfully in ", iter, " iterations.")
  }
  
  cat("Total Runtime:", total_runtime, "seconds\n")
  
  return(list(beta_m = beta_m, data_weted = data_weted))
}


# Define initial guesses for beta_m parameters
beta_m_0 = c(b0 = -2.75, A = 1.75, L1 = 2, L2 = 1.25, L3 = -0.35, sd = 0.35)
LOD_value <- as.numeric(data$lod)  # Set the LOD value
# Run the EM algorithm
dat_res <- EM_algorithm(data = dataset, beta_m_0 = beta_m_0)

dat_aug <- dat_res$data_weted
dat_res$beta_m




gcomp.nde <- function(data) {
  # Extract variables from the dataset
  y <- data$Y        # Assuming Y is the outcome variable
  m <- data$M        # Assuming M is the mediator
  a <- data$A        # Assuming A is the treatment
  l1 <- data$L1      # Assuming L1 is the first covariate
  l2 <- data$L2      # Assuming L2 is the second covariate
  l3 <- data$L3      # Assuming L3 is the third covariate
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  weights <- data$w_norm  # Extract the weights
  
  # Combine L1, L2, L3 into a data frame for modeling convenience
  covariates <- data.frame(l1 = l1, l2 = l2, l3 = l3)
  
  # Fit a weighted linear model for y including m, a, and covariates
  lm_y <- glm(y ~ m + a + l1 + l2 + l3 + m*a, weights = weights, family = binomial)
  
  # Predict potential outcomes under different treatments
  pred_y1 <- predict(lm_y, newdata = data.frame(l1 = l1, l2 = l2, l3 = l3, a = 1, m = m), type = "response")
  pred_y0 <- predict(lm_y, newdata = data.frame(l1 = l1, l2 = l2, l3 = l3, a = 0, m = m), type = "response")
  
  # Fit a weighted linear model using the predicted outcome
  lm_y1 <- glm(pred_y1 ~ a + l1 + l2 + l3 + a*l1 + a*l2 + a*l3, 
               weights = weights, family = quasibinomial(link = "logit"))
  lm_y0 <- glm(pred_y0 ~ a + l1 + l2 + l3 + a*l1 + a*l2 + a*l3, 
               weights = weights, family = quasibinomial(link = "logit"))
  
  # Predict the causal effect when a = 0
  y10 <- predict(lm_y1, newdata = data.frame(l1 = l1, l2 = l2, l3 = l3, a = 0), type = "response")
  y00 <- predict(lm_y0, newdata = data.frame(l1 = l1, l2 = l2, l3 = l3, a = 0), type = "response")
  
  # Calculate and return the estimate
  nde.rd <- mean(y10 - y00); nde.rr <- mean(y10)/mean(y00); nde.or <- (mean(y10)/(1-mean(y10))) / (mean(y00)/(1-mean(y00)))
  return(list(nde.rd = nde.rd, nde.rr = nde.rr, nde.or = nde.or))
}


# Perform bootstrap
gcom_boots_ci <- function(data, n_bootstrap = 5000, alpha = 0.05, 
                          type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  type <- match.arg(type)       # Ensure valid type (NDE or NIE)
  measure <- match.arg(measure) # Ensure valid measure (RD, RR, OR)
  
  n <- nrow(data)
  boots_est <- numeric(n_bootstrap)
  
  for (i in seq_len(n_bootstrap)) {
    # Resample data with replacement
    bootstrap_sample <- data[sample(seq_len(n), size = n, replace = TRUE), ]
    
    # Compute bootstrap estimate based on type and measure
    if (type == "NDE") {
      bres <- gcomp.nde(bootstrap_sample)
      boots_est[i] <- switch(measure, 
                             RD = bres$nde.rd, 
                             RR = bres$nde.rr, 
                             OR = bres$nde.or)
    } else if (type == "NIE") {
      bres <- gcomp.nie(bootstrap_sample)
      boots_est[i] <- switch(measure, 
                             RD = bres$nie.rd, 
                             RR = bres$nie.rr, 
                             OR = bres$nie.or)
    }
  }
  
  # Remove any potential NA values from computation
  boots_est <- boots_est[!is.na(boots_est)]
  
  # Calculate the percentile-based confidence interval
  ci_lower <- quantile(boots_est, probs = alpha / 2, na.rm = TRUE)
  ci_upper <- quantile(boots_est, probs = 1 - alpha / 2, na.rm = TRUE)
  
  list(
    point_estimate = mean(boots_est, na.rm = TRUE),  # Mean of bootstrap estimates
    ci = c(lower = ci_lower, upper = ci_upper)  # Named confidence interval
  )
}


# Example usage with your dataset
# (res1 <- gcom_boots_ci(dat_aug, n_bootstrap = 200, type = "NDE", measure = "RD"))

gcomp.nie <- function(data) {
  # Extract variables from the dataset
  y <- data$Y        # Outcome variable
  m <- data$M        # Mediator
  a <- data$A        # Treatment
  l1 <- data$L1      # Covariate 1
  l2 <- data$L2      # Covariate 2
  l3 <- data$L3      # Covariate 3
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  weights <- data$w_norm  # Weights for weighted regression
  
  # Combine covariates into a data frame for modeling convenience
  covariates <- data.frame(l1 = l1, l2 = l2, l3 = l3)
  
  # Fit a weighted linear model for y including m, a, and covariates
  lm_y <- glm(y ~ m + a + l1 + l2 + l3 + m*a, weights = weights, family = binomial)
  
  # Predict potential outcomes under different treatments
  pred_y1 <- predict(lm_y, newdata = transform(covariates, a = 1, m = m),  type = "response")
  
  # Fit a weighted linear model using the predict potential outcomes
  lm_y1 <- glm(pred_y1 ~ a + l1 + l2 + l3, weights = weights, family = quasibinomial(link = "logit"))
  
  # Predict the causal effect when a = 0
  y11 <- predict(lm_y1, newdata = transform(covariates, a = 1),  type = "response")
  y10 <- predict(lm_y1, newdata = transform(covariates, a = 0),  type = "response")
  
  # Calculate and return the estimate
  nie.rd <- mean(y11 - y10); nie.rr <- mean(y11)/mean(y10); nie.or <- (mean(y11)/(1-mean(y11))) / (mean(y10)/(1-mean(y10)))
  # Return the estimate of the indirect effect
  return(list(nie.rd = nie.rd, nie.rr = nie.rr, nie.or = nie.or))
}


MC_sim <- function(r, n_obs, type = c("Trad", "EM-MLE"), censor_rate, S = 2, beta_m) {
  type <- match.arg(type)
  output_file <- paste0("/n/home09/c55jiang/mediation/FEB19_B5", type, "_n_", n_obs, "rate", censor_rate, ".txt")
  
  for (i in 1:r) {
    data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = censor_rate)
    simdata<- data$study_dat
    LOD_value <- as.numeric(data$lod)
    
    if (type == "Trad") {
      data_aug <- NaiveImpute(simdata, LOD_value)
    } else if (type == "EM-MLE") {
      em_result <- EM_algorithm(simdata, beta_m, S=S, LOD = LOD_value)
      data_aug <- em_result$data_weted
      beta_m <- em_result$beta_m  # Persisting update
    }
    
    # G-computation
    nde_gcom_glm <- gcomp.nde(data_aug)
    CI.nde <- gcom_boots_ci(data_aug, n_bootstrap = 200, type = "NDE", measure = "RD")$ci
    nie_gcom_glm <- gcomp.nie(data_aug)
    CI.nie <- gcom_boots_ci(data_aug, n_bootstrap = 200, type = "NIE", measure = "RD")$ci
    
    # Save results
    results <- c(i, 
                 nde_gcom_glm$nde.rd, CI.nde, 
                 nie_gcom_glm$nie.rd, CI.nie, beta_m)
    
    write(results, file = output_file, ncolumns = length(results), append = TRUE)
  }
  
  return(beta_m)  # Returning the final beta_m if needed
}


# Example usage
r <- 1000   # Number of replicates
n_obs <- 1000  # Number of observations per replicate
beta_m <- c(b0 = -2.75, A = 1.75, L1 = 2, L2 = 1.25, L3 = -0.35, sd = 0.35)


#sim.res <- MC_sim(r = r, n_obs = n_obs, type = "Trad", censor_rate = 0.75,  beta_m = beta_m)
sim.res <- MC_sim(r = r, n_obs = n_obs, type = "EM-MLE", censor_rate = 0.1,  S = 5, beta_m = beta_m)





