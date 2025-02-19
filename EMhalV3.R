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
n_obs <- 300

(data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = 0.5))
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
folds <- origami::make_folds(n = nrow(data_impu), V = 3)
fold <- folds[[1]]
cv_hal_joint <- function(fold, data_in,
                         x_names, y_cen_names, y_out_names, weights_names,
                         lambda_seq1 = exp(seq(-0.5, -20, length = 1000)),
                         lambda_seq2 = exp(seq(-0.5, -20, length = 1000)),
                         basis_list) {
  ## 0) set training and validation folds for cross-validation via origami
  train_data <- origami::training(data_in)
  valid_data <- origami::validation(data_in)
  
  # extract censoring status, outcome, adjustment set, weights for given split
  x_train <- as.matrix(train_data[, ..x_names])
  y_cen_train <- as.numeric(train_data[, get(y_cen_names)])
  y_out_train <- as.numeric(train_data[, get(y_out_names)])
  weights_train <- as.numeric(train_data[, get(weights_names)])
  x_valid <- as.matrix(valid_data[, ..x_names])
  y_cen_valid <- as.numeric(valid_data[, get(y_cen_names)])
  y_out_valid <- as.numeric(valid_data[, get(y_out_names)])
  weights_valid <- as.numeric(valid_data[, get(weights_names)])
  
  ## 1) fit HAL over sequence of lambda values
  cen_mod <- hal9001::fit_hal(
    X = x_train, Y = y_cen_train, weights = weights_train, family = "binomial",
    lambda = lambda_seq1,
    # max_degree = NULL,
    basis_list = basis_list,
    fit_control = list(cv_select = FALSE)
  )
  
  out_mod <- hal9001::fit_hal(
    X = x_train, Y = y_out_train, weights = weights_train, family = "binomial",
    lambda = lambda_seq2,
    # reduce_basis = 2 / sqrt(nrow(x_train)),
    basis_list = basis_list,
    fit_control = list(cv_select = FALSE)
  )
  # get coefficients
  coef_mat_cen <- cen_mod$coefs #     cen_coef = coef(cen_mod)
  coef_mat_out <- out_mod$coefs #     out_coef = coef(out_mod)
  
  ## 2) predictions on validation data for each value of lambda
  valid_x_basis <- hal9001::make_design_matrix(x_valid, basis_list)
  valid_x_basis_clean <- valid_x_basis[, as.numeric(names(cen_mod$copy_map))]
  pred_mat <- cbind(rep(1, nrow(x_valid)), valid_x_basis_clean)
  
  # preds_valid <- as.matrix(pred_mat %*% coef_mat)
  hat_val_cen <- as.matrix(pred_mat %*% coef_mat_cen)
  hat_valid_cen <- apply(hat_val_cen, 2, stats::plogis) # for binary censoring indicator
  hat_val_out <- as.matrix(pred_mat %*% coef_mat_out)
  hat_valid_out <- apply(hat_val_out, 2, stats::plogis) # for binary outcome
  # OR ???
  # hat_valid_cen <- predict(cen_mod, new_data = x_val, type = "response")
  # hat_valid_out <- predict(out_mod, new_data = x_val, type = "response")
  
  return(list(
    hat_valid_cen = hat_valid_cen,
    hat_valid_out = hat_valid_out
  ))
}


# Helper function: Define negative log-likelihood function
# NOTE: using sl3::loss_loglik_binomial() below instead
nll <- function(obs, probs) {
  obs * log(probs) + (1 - obs) * log(1 - probs)
}


fit_cv_hal_joint <- function(data_in, folds,
                             x_names, y_cen_names, y_out_names, weights_names,
                             lambda_seq1, lambda_seq2,
                             basis_list) {
  # fit a cross-validated joint HAL for the Censoring and Outcome Models
  # NOTE: use set of basis functions discovered in initial enumerate_basis
  cv_fit_hal_joint <- origami::cross_validate(
    cv_fun = cv_hal_joint,
    folds = folds,
    data = data_in,
    # arguments passed to cv_fun
    x_names = x_names,
    y_cen_names = y_cen_names,
    y_out_names = y_out_names,
    weights_names = weights_names,
    lambda_seq1 = lambda_seq1,
    lambda_seq2 = lambda_seq2,
    basis_list = basis_list,
    # back to arguments for cross_validate
    use_future = FALSE,
    .combine = FALSE
  )
  
  # Extract validation indices from folds to reorder predictions
  idx_folds <- do.call(c, lapply(folds, `[[`, "validation_set"))
  
  lapply(cv_fit_hal_joint$hat_valid_out, dim)  # Check dimensions
  lapply(cv_fit_hal_joint$hat_valid_cen, dim)
  
  str(cv_fit_hal_joint$hat_valid_out)
  sapply(cv_fit_hal_joint$hat_valid_cen, ncol)  # Verify column consistency
  
  
  # Combine predictions and reorder based on validation-set indices
  cen_cv <- do.call(rbind, cv_fit_hal_joint$hat_valid_cen)[order(idx_folds), ]
  out_cv <- do.call(rbind, cv_fit_hal_joint$hat_valid_out)[order(idx_folds), ]
  
  # Compute CV-NLL(O_i) for each lambda in lambda_seq1 for censoring model
  cv_nll_cen <- apply(cen_cv, 2, function(cen_hat) {
    loss_loglik_binomial(cen_hat, data_in[, get(y_cen_names)])
  })
  # Compute CV-NLL(O_i) for each lambda in lambda_seq2 for outcome model
  cv_nll_out <- apply(out_cv, 2, function(out_hat) {
    loss_loglik_binomial(out_hat, data_in[, get(y_out_names)])
  })
  
  # get each best idx, separately
  best_lambda_cen_idx <- which.min(colMeans(cv_nll_cen))
  best_lambda_out_idx <- which.min(colMeans(cv_nll_out))
  
  # create a grid of all (lambda_seq1, lambda_seq2) combinations
  # compute sum of empirical risk for each pair (lambda1, lambda2)
  cv_nll_sum <- outer(colMeans(cv_nll_cen), colMeans(cv_nll_out), "+")
  
  # find the indices of the best combination (minimizing the sum)
  best_lambda_idx <- which(cv_nll_sum == min(cv_nll_sum), arr.ind = TRUE)
  best_lambda_idx <- as_tibble(best_lambda_idx)
  colnames(best_lambda_idx) <- c("cen", "out")
  return(list(best_lambda_idx = best_lambda_idx, best_lambda_cen_idx = best_lambda_cen_idx,
              best_lambda_out_idx = best_lambda_out_idx))
}


mod_update_halcv <- function(data, folds, lambda_seq1 = exp(seq(-0.5, -20, length = 1000)),
                             lambda_seq2 = exp(seq(-0.5, -20, length = 1000)),
                             basis_list = NULL) {
  if (is.null(data$w_norm)) {
    data$w_norm <- 1
  }
  
  x_names <- c("L1", "L2", "L3", "A", "M")
  y_cen_names <- "C"
  y_out_names <- "Y"
  weights_names <- "w_norm"
  
  # fit enumerate_basis on full data to get set of basis functions
  if (is.null(basis_list)) {
    basis_list <- hal9001::enumerate_basis(
      as.matrix(data[, ..x_names]),
      max_degree = 2, smoothness_orders = 1
    )
  }
  
  best_lambda_idx <- fit_cv_hal_joint(
    data_in = data, folds,
    x_names, y_cen_names, y_out_names, weights_names,
    lambda_seq1, lambda_seq2,
    basis_list
  )
  
  best_lambda_cen_idx <- best_lambda_idx$best_lambda_idx$cen[1]  # Extract first row if multiple
  best_lambda_out_idx <- best_lambda_idx$best_lambda_idx$out[1]
  
  # Update HAL models with best lambda
  cen_mod <- fit_hal(
    X = as.matrix(data[, ..x_names]),
    Y = as.numeric(data[, get(y_cen_names)]),
    weights = as.numeric(data[, get(weights_names)]),
    family = "binomial",
    lambda = lambda_seq1[best_lambda_cen_idx],
    # max_degree = 3, reduce_basis = 1 / sqrt(nrow(X)), smoothness_orders = 0,
    basis_list = basis_list, fit_control = list(cv_select = FALSE)
  )
  
  out_mod <- fit_hal(
    X = as.matrix(data[, ..x_names]),
    Y = as.numeric(data[, get(y_out_names)]),
    weights = as.numeric(data[, get(weights_names)]),
    family = "binomial",
    lambda = lambda_seq2[best_lambda_out_idx],
    # max_degree = 3, reduce_basis = 1 / sqrt(nrow(X)), smoothness_orders = 0,
    basis_list = basis_list, fit_control = list(cv_select = FALSE)
  )
  
  
  out_probY1 <- predict(out_mod, new_data = as.matrix(data[, ..x_names]), type = "response")
  out_probY0 <- 1 - out_probY1
  
  # Create a new column for predicted probabilities based on observed Y values
  out_prob <- ifelse(data$Y == 1, out_probY1, out_probY0)
  
  return(list(
    cen_mod = cen_mod,
    out_mod = out_mod,
    cen_coef = coef(cen_mod),
    out_coef = coef(out_mod),
    cen_prob = predict(cen_mod, new_data = as.matrix(data[, ..x_names]), type = "response"),
    out_prob = out_prob,
    basis_list = basis_list,
    best_lambda_idx = best_lambda_idx
  ))
}

system.time({
  cv_halmod_pred <- mod_update_halcv(data = data_impu, folds)
})
# Access results
cv_halmod_pred$cen_prob # Predicted probabilities for the censoring model
cv_halmod_pred$out_prob # Predicted probabilities for the outcome model

cv_halmod_pred$cen_coef # Predicted probabilities for the censoring model
cv_halmod_pred$out_coef # Predicted probabilities for the outcome model

cv_halmod_pred$best_lambda # Optimal lambda selected
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
# beta_m_start <- c(b0 = -2.75, A = 1.75, L1 = 2, L2 = 1.25, L3 = -0.35, sd = 0.5)
#M_stepf_M(data = data_weted, beta_m = beta_m_start)

# EM Algorithm
EM_algorithm <- function(data, beta_m_0, S = 2, LOD = LOD_value, max_iter = 50, tol = 1e-20, verbose = TRUE) {
  iter <- 0; diff <- Inf; rel_diff <- Inf; beta_m <- beta_m_0 
  
  data_imputed <- Impute_step(data, beta_m_0, S, LOD)
  mod_pred <- mod_update_glm(data_imputed)
  
  start_time <- Sys.time()  # Track runtime
  
  while (iter < max_iter && diff > tol && rel_diff > tol) {
    iter_time <- Sys.time()
    
    # W-Step
    data_weted <- Wet_step(data_imputed, mod_pred, beta_m_new = beta_m, beta_m_0 = beta_m_0)
    
    # Maximization Step
    beta_m_new <- M_stepf_M(data = data_weted, beta_m)
    
    # Compute convergence criteria
    diff <- sum((beta_m_new - beta_m)^2)
    
    cen_coef_old <- mod_pred$cen_coef
    out_coef_old <- mod_pred$out_coef
    
    # Parameter update
    beta_m <- beta_m_new
    mod_pred <- mod_update_glm(data_weted)
    
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
gcom_boots_ci <- function(data, n_bootstrap = 2000, alpha = 0.05, 
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
  output_file <- paste0("/n/home09/c55jiang/mediation/MED_FEB17_", type, "_n_", n_obs, "rate", censor_rate, ".txt")
  
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
sim.res <- MC_sim(r = r, n_obs = n_obs, type = "EM-MLE", censor_rate = 0.5,  S = 5, beta_m = beta_m)





