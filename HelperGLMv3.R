remotes::install_github("nhejazi/medoutcon")
devtools::install_github("tlverse/tlverse")
library(data.table) # Assuming 'data' is a data.table
library(dplyr)
library(sl3)
library(hal9001)
library(medoutcon)
library(truncnorm)
library(xgboost)
library(future.apply)
library(dplyr) 
#n_obs <- 100000



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


#NaiveImpute(dataset, LOD = LOD_value)

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
  data$C <- as.numeric(data$C)
  return(data[order(data$id), ])
}

#(data_impu <- Impute_step(data = dataset, beta_m, S = 10, LOD = LOD_value))



mod_update_glm <- function(data) {
  # Check and assign w_norm if NULL
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  # Censoring Modeling (parametric model, does not include Y)
  #cen_mod <- glm(C ~ L1 + L2 + L3 + A + M + A*M, data = data,family = binomial, weights = w_norm)
  
  #cen_prob <- predict(cen_mod, type = "response")
  # Outcome Modeling
  out_mod <- glm(Y ~ L1 + L2 + L3 + A + M + A*M, data = data, family = binomial, weights = w_norm)
  out_probY1 <- predict(out_mod, type = "response")
  out_probY0 <- 1 - out_probY1
  
  # Create a new column for predicted probabilities based on observed Y values
  out_prob <- ifelse(data$Y == 1, out_probY1, out_probY0)
  
  # Return predictions and coefficients
  return(list(
    #cen_prob = cen_prob, 
    out_prob = out_prob,
    #cen_coef = cen_mod$coefficients,
    out_coef = out_mod$coefficients
  ))
}



# mod_pred are from mod_update_glm funciton, which inlucdes both out_prob and cen_prob
Wet_step <- function(data_imputed, mod_pred, beta_m_new, beta_m_0) {
  # Update probabilities
  #data_imputed$cen_prob <- pmin(pmax(mod_pred$cen_prob, 1e-40), 1 - 1e-40)
  data_imputed$out_prob <- mod_pred$out_prob
  
  # Filter censored data
  censored_data <- data_imputed[data_imputed$C == 0, ]
  
  # Compute new and initial weights for censored data
  cenwet_new <- compute_f_M(censored_data, beta_m_new)
  cenwet_int <- compute_f_M(censored_data, beta_m_0)
  
  censored_data$w <- censored_data$out_prob * (cenwet_new / cenwet_int)
  
  # Normalize weights by ID
  censored_data$w_norm <- ave(censored_data$w, censored_data$id, FUN = function(x) x / sum(x))
  
  # For uncensored data
  uncensored_data <- data_imputed[data_imputed$C == 1, ]
  uncensored_data$w_norm <- 1;   uncensored_data$w <- 1
  data_imputed <- rbind(censored_data, uncensored_data)
  
  # Remove temporary columns
  data_weted <- data_imputed[, !names(data_imputed) %in% c("out_prob"), with = FALSE]
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
    diff <- diff + sum((mod_pred$out_coef - out_coef_old)^2)
    
    # Fixed denominator issue in rel_diff calculation
    rel_diff <- diff / (sum(beta_m^2) + sum(out_coef_old^2) + 1e-8)  # Prevent division by zero
    # rel_diff <- diff / (sum(beta_m^2) + 1e-8)  # Prevent division by zero
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
library(future.apply)
gcom_boots_ci <- function(data, n_bootstrap = 1000, alpha = 0.05, 
                          type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  
  # Parallel bootstrap estimation using future_lapply()
  boots_est <- future_lapply(seq_len(n_bootstrap), function(i) {
    bootstrap_sample <- data[sample(nrow(data), size = nrow(data), replace = TRUE), ]
    
    # Compute bootstrap estimate based on type and measure
    bres <- if (type == "NDE") gcomp.nde(bootstrap_sample) else gcomp.nie(bootstrap_sample)
    
    # Extract the corresponding estimate
    if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
    if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
    if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
    
    return(NA)  # Handle unexpected measure values safely
  })
  
  # Convert to numeric vector and remove NAs
  boots_est <- na.omit(unlist(boots_est))
  
  # Check for empty bootstrap results
  if (length(boots_est) == 0) {
    warning("All bootstrap estimates are NA. Check data or computation methods.")
    return(NULL)
  }
  
  # Compute confidence intervals
  ci_lower <- quantile(boots_est, probs = alpha / 2, na.rm = TRUE)
  ci_upper <- quantile(boots_est, probs = 1 - alpha / 2, na.rm = TRUE)
  
  list(
    point_estimate = mean(boots_est, na.rm = TRUE),
    variance = var(boots_est, na.rm = TRUE),
    ci = c(lower = ci_lower, upper = ci_upper)
  )
}



#######
#m-out-of-n bootstrap confidence intervals
#######
library(future.apply)
# Set up the parallel plan (e.g., multisession for all systems or multicore for Unix-like systems)
plan(multisession)  # Use multisession for portability across all systems

# Function for m-out-of-n bootstrap confidence intervals
gcom_boots_ci_m_out_n <- function(data, m, n_bootstrap = 1000, alpha = 0.05, 
                                  type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  
  boots_est <- numeric(n_bootstrap)
  
  # Parallel bootstrap estimation using future_lapply
  boots_est <- future_lapply(seq_len(n_bootstrap), function(i) {
    bootstrap_sample <- data[sample(nrow(data), size = m, replace = TRUE), ]
    
    # Compute bootstrap estimate based on type and measure
    bres <- if (type == "NDE") gcomp.nde(bootstrap_sample) else gcomp.nie(bootstrap_sample)
    
    # Extract the corresponding estimate
    if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
    if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
    if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
    
    return(NA)  # Handle unexpected measure values safely
  }, future.seed = TRUE)  # Ensure reproducibility
  
  boots_est <- na.omit(unlist(boots_est))  # Flatten and remove NAs
  
  if (length(boots_est) == 0) {
    warning("All bootstrap estimates are NA. Check data or computation methods.")
    return(NULL)
  }
  
  # Scale variance by (n/m) for m-out-of-n bootstrap correction
  point_estimate <- mean(boots_est, na.rm = TRUE)
  scaled_variance <- var(boots_est, na.rm = TRUE) * (m / nrow(data))
  ci <- quantile(boots_est, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  
  list(
    point_estimate = point_estimate,
    variance = scaled_variance,
    ci = ci
  )
}

# Function to compute effect estimates
compute_eff <- function(data, indices, type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  sample_data <- data[indices, ]  # Resample data
  bres <- if (type == "NDE") gcomp.nde(sample_data) else gcomp.nie(sample_data)
  
  # Extract the corresponding estimate
  if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
  if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
  if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
}


# Function to choose c_i using double bootstrap ideas (parallelized with early stopping)
choose_ci <- function(data, p_cen = 0.5, alpha_seq = seq(0, 2, length.out = 100), B1 = 100, B2 = 50, nominal_coverage = 0.95,
                      type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  n <- nrow(data)  # Original sample size
  original_estimate <- compute_eff(data, 1:n, type = type, measure = measure)  # Estimate from original data
  
  # Initialize variables
  coverage_rates <- numeric(length(alpha_seq))
  optimal_c <- NA
  opt_m <- NA
  early_stop <- FALSE  # Flag to indicate early stopping
  
  # Loop over alpha_seq manually to allow early stopping
  for (i in seq_along(alpha_seq)) {
    if (early_stop) break  # Exit the loop if early stop is triggered
    
    alpha <- alpha_seq[i]
    
    # Compute c_i
    c_i <- (1 + alpha * (1 - p_cen)) / (1 + alpha)
    m_i <- floor(n^c_i)  # Compute m_i
    
    # Parallelize the first-stage bootstrap (B1 iterations)
    coverage_rate <- future_sapply(1:B1, function(b1) {
      first_stage_indices <- sample(1:n, n, replace = TRUE)  # Resample data
      first_stage_estimate <- compute_eff(data, first_stage_indices, type = type, measure = measure)
      
      # Second-stage (nested) bootstrap
      nested_res <- gcom_boots_ci_m_out_n(data[first_stage_indices, ], m_i, n_bootstrap = B2, type = type, measure = measure)
      nested_ci <- nested_res$ci
      
      # Check if the original estimate falls within the CI
      if (original_estimate >= nested_ci[1] && original_estimate <= nested_ci[2]) {
        return(1)  # Coverage success
      } else {
        return(0)  # Coverage failure
      }
    }, future.seed = TRUE)  # Ensure reproducibility
    
    # Calculate the coverage rate for this alpha
    coverage_rates[i] <- mean(coverage_rate)
    
    # Check if the coverage rate exceeds the nominal coverage
    if (coverage_rates[i] >= nominal_coverage) {
      optimal_c <- alpha
      opt_m <- m_i
      early_stop <- TRUE  # Trigger early stopping
      cat("Stopping early at alpha =", alpha, "with coverage rate =", coverage_rates[i], "\n")
    }
    
    cat("Iteration", i, "alpha =", alpha, "coverage rate =", coverage_rates[i], "\n")
  }
  
  # If no alpha achieved the nominal coverage, select the one with the highest coverage rate
  if (is.na(optimal_c)) {
    optimal_c <- alpha_seq[which.max(coverage_rates)]
    opt_m <- floor(n^optimal_c)
    cat("No alpha achieved the nominal coverage. Selecting alpha with the highest coverage rate:", optimal_c, "\n")
  }
  
  return(list(
    optimal_c = optimal_c,
    coverage_rates = data.frame(c = alpha_seq, coverage_rate = coverage_rates),
    opt_m = opt_m
  ))
}


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




sim_iteration <- function(i, n_obs, type = c("Trad", "EM-MLE", "Oracle"), censor_rate, S, beta_m) {
  type <- match.arg(type)
  output_file <- paste0("/n/home09/c55jiang/mediation/res/MAR11_moutn", type, "_n_", n_obs, "rate", censor_rate, ".txt")
  data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = censor_rate)
  simdata <- data$study_dat
  LOD_value <- as.numeric(data$lod)
  
  if (type == "Trad") {
    data_aug <- NaiveImpute(simdata, LOD_value)
  } else if (type == "EM-MLE") {
    em_result <- EM_algorithm(simdata, beta_m, S = S, LOD = LOD_value)
    data_aug <- em_result$data_weted
    beta_mres <- em_result$beta_m  
  } else if (type == "Oracle") {
    data_aug <- data$study_dat_full
  }
  
  # G-computation
  nde_gcom_glm <- gcomp.nde(data_aug)
  
  choose_ciNDE <- choose_ci(data_aug, p_cen = censor_rate, type = "NDE", measure = "RD")
  CI.nderes <- gcom_boots_ci_m_out_n(data_aug, m = choose_ciNDE$opt_m, n_bootstrap = 1000, type = "NDE", measure = "RD")
  CI.nde <- CI.nderes$ci
  var.nde <- CI.nderes$variance
  
  
  nie_gcom_glm <- gcomp.nie(data_aug)
  
  choose_ciNIE <- choose_ci(data_aug, p_cen = censor_rate, type = "NIE", measure = "RD")
  CI.nieres <- gcom_boots_ci_m_out_n(data_aug, m = choose_ciNIE$opt_m, n_bootstrap = 1000, type = "NIE", measure = "RD")
  CI.nie <- CI.nieres$ci
  var.nie <- CI.nieres$variance
  
  
  # Save results
  results <- c(i,nde_gcom_glm$nde.rd, CI.nde, var.nde, 
               nie_gcom_glm$nie.rd, CI.nie, var.nie,  beta_mres)
  
  write(results, file = output_file, ncolumns = length(results), append = TRUE)
  return(NULL)  # Returning NULL keeps memory usage low
}


parallelly::availableCores()
plan(list(
  tweak(multisession, workers = availableCores() %/% 10L),
  tweak(multisession, workers = I(10L))
))

r <- 1000   # Number of replicates
n_obs <- 300  # Number of observations per replicate
beta_m <- c(b0 = -2.75, A = 1.75, L1 = 2, L2 = 1.25, L3 = -0.35, sd = 0.35)
results <- future_lapply(
  1:r, 
  function(i) sim_iteration(i, 
                            n_obs = n_obs, 
                            type = "EM-MLE", 
                            censor_rate = 0.5, 
                            S = 5, 
                            beta_m = beta_m)
)

