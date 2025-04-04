remotes::install_github("nhejazi/medoutcon")
devtools::install_github("tlverse/tlverse")
library(data.table) # Assuming 'data' is a data.table
library(dplyr)
library(sl3)
library(hal9001)
library(medoutcon)
library(truncnorm)
library(xgboost)
library(zoo) 
n_obs <- 1000

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


library(R6)
Lrnr_density_hse <- R6Class(
  classname = "Lrnr_density_hse",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(mean_learner = NULL, ...) {
      if (is.null(mean_learner)) {
        mean_learner <- make_learner(Lrnr_glm_fast)
      }
      params <- list(mean_learner = mean_learner, ...)
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("density"),
    .train = function(task) {
      mean_learner <- self$params$mean_learner
      mean_fit <- mean_learner$train(task)
      
      # TODO: maybe these should be cv errors?
      mean_preds <- mean_fit$predict()
      errors <- task$Y - mean_preds
      dens_fit <- density(errors)
      fit_object <- list(mean_fit = mean_fit, dens_fit = dens_fit)
      return(fit_object)
    },
    .predict = function(task) {
      mean_fit <- self$fit_object$mean_fit
      dens_fit <- self$fit_object$dens_fit
      mean_preds <- mean_fit$predict(task)
      errors <- task$Y - mean_preds
      dens_preds <- approx(dens_fit$x, dens_fit$y, errors, rule = 2)$y
      # dens_preds[is.na(dens_preds)] <- 0
      return(list(dens_preds = dens_preds, errors = errors))
    },
    .required_packages = c()
  )
)



# Define function to compute estimated density f(M) using SL.
compute_f_M_HAL <- function(data) {
  if (is.null(data$w_norm)) {
    data$w_norm <- 1
  }
  # Convert new data to an sl3 task
  dens_task <- make_sl3_Task(
    data = data,
    covariates = c("L1", "L2", "L3", "A"),
    outcome = "M",
    weights = "w_norm"
  )
  
  mean_learner1 <- Lrnr_hal9001$new()   
  lrnr_density_hse <- Lrnr_density_hse$new(mean_learner = mean_learner1)
  fit_density_hse <- lrnr_density_hse$train(dens_task)
  preds_density_hse <- fit_density_hse$predict()
  
  return(list(f_M = preds_density_hse$dens_preds, errors = preds_density_hse$errors))
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

#system.time({
#  halmod_pred <- mod_update_hal(data = data_impu)
#})
# Access results
#halmod_pred$out_prob # Predicted probabilities for the outcome model
#halmod_pred$basis_list # Shared basis functions


Wet_step <- function(data_imputed, mod_pred, beta_m_0) {
  # Update probabilities
  #data_imputed$cen_prob <- pmin(pmax(mod_pred$cen_prob, 1e-40), 1 - 1e-40)
  data_imputed$out_prob <- mod_pred$out_prob
  
  Mdens_HAL <- compute_f_M_HAL(data_imputed)
  data_imputed$Mdens_new <- Mdens_HAL$f_M
  Mdens_error <- Mdens_HAL$errors
  
  # Filter censored data
  censored_data <- data_imputed[data_imputed$C == 0, ]
  cendens_new <- censored_data$Mdens_new
  cendens_int <- compute_f_M(censored_data, beta_m_0)
  
  censored_data$w <- censored_data$out_prob * (cendens_new / cendens_int)
  
  # Normalize weights by ID
  censored_data$w_norm <- ave(censored_data$w, censored_data$id, FUN = function(x) x / sum(x))
  
  # For uncensored data
  uncensored_data <- data_imputed[data_imputed$C == 1, ]
  uncensored_data$w_norm <- 1
  uncensored_data$w <- 1
  data_imputed <- rbind(censored_data, uncensored_data)
  
  # Remove temporary columns
  data_weted <- data_imputed[, !names(data_imputed) %in% c("out_prob"), with = FALSE]
  return(list(wet_dat = data_weted[order(data_weted$id),], dens_errors = Mdens_error))
}
# (data_weted <- Wet_step(data_imputed = data_impu, mod_pred = halmod_pred, beta_m_0 = beta_m))



# EM Algorithm
EM_algorithm <- function(data, beta_m_0, S = 2, LOD = LOD_value, max_iter = 100, tol = 1e-8, verbose = TRUE) {
  iter <- 0; diff <- Inf; rel_diff <- Inf
  
  # I-Step: Initial Imputation
  data_imputed <- Impute_step(data, beta_m_0, S, LOD)
  halmod_pred <- mod_update_hal(data_imputed)
  halbasis <- halmod_pred$basis_list
  wet_step <- Wet_step(data_imputed = data_imputed, mod_pred = halmod_pred, beta_m_0 = beta_m_0)
  data_weted <- wet_step$wet_dat
  
  start_time <- Sys.time()  # Track total runtime
  
  while (iter < max_iter && diff > tol && rel_diff > tol) {
    iter_time <- Sys.time()
    
    # W-Step
    wet_step <- Wet_step(data_imputed = data_weted, mod_pred = halmod_pred, beta_m_0 = beta_m_0)
    
    data_weted <- wet_step$wet_dat
    errors <- wet_step$dens_errors
    
    # M-Step
    errors_new <- compute_f_M_HAL(data = data_weted)$errors
    
    # Compute Convergence Criteria
    diff <- sum((errors_new - errors)^2)
    
    out_coef_old <- halmod_pred$out_coef
    
    # Parameter Update
    halmod_pred <- mod_update_hal(data_weted, basis_list = halbasis)
    
    # Update diff with coefficient changes
    #diff <- diff + sum((halmod_pred$out_coef - out_coef_old)^2)
    
    # Fixed denominator issue in rel_diff calculation
    # rel_diff <- diff / (sum(errors^2) + sum(out_coef_old^2) + 1e-8)  # Prevent division by zero
    rel_diff <- diff / (sum(errors^2) + 1e-8)  # Prevent division by zero
    
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
  
  return(list(data_weted = data_weted))
}


# Define initial guesses for beta_m parameters
beta_m_0 = c(b0 = -2.75, A = 1.75, L1 = 2, L2 = 1.25, L3 = -0.35, sd = 0.35)
LOD_value <- as.numeric(data$lod)  # Set the LOD value
# Run the EM algorithm
#dat_res <- EM_algorithm(data = dataset, beta_m_0 = beta_m_0)

#dat_aug <- dat_res$data_weted





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


library(parallel)
gcom_boots_ci_m_out_n <- function(data, n_bootstrap = 1000, alpha = 0.05, 
                                  type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  
  n <- nrow(data)
  m <- floor(n^0.9)  # Adjusted bootstrap sample size
  boots_est <- numeric(n_bootstrap)
  
  num_cores <- detectCores() - 1  # Use available CPU cores
  
  # Parallel bootstrap estimation
  boots_est <- mclapply(seq_len(n_bootstrap), function(i) {
    bootstrap_sample <- data[sample(n, size = m, replace = TRUE), ]
    
    # Compute bootstrap estimate based on type and measure
    bres <- if (type == "NDE") gcomp.nde(bootstrap_sample) else gcomp.nie(bootstrap_sample)
    
    # Extract the corresponding estimate
    if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
    if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
    if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
    
    return(NA)  # Handle unexpected measure values safely
  }, mc.cores = num_cores)
  
  boots_est <- na.omit(unlist(boots_est))  # Flatten and remove NAs
  
  if (length(boots_est) == 0) {
    warning("All bootstrap estimates are NA. Check data or computation methods.")
    return(NULL)
  }
  
  
  # Scale variance by (n/m) for m-out-of-n bootstrap correction
  point_estimate <- mean(boots_est, na.rm = TRUE)
  scaled_variance <- var(boots_est, na.rm = TRUE) * (m / n)
  ci_lower <- quantile(boots_est, probs = alpha / 2, na.rm = TRUE) 
  ci_upper <- quantile(boots_est, probs = 1 - alpha / 2, na.rm = TRUE) 
  
  list(
    point_estimate = point_estimate,
    variance = scaled_variance,
    ci = c(lower = ci_lower, upper = ci_upper)
  )
}



# Example usage with your dataset
#(res1 <- gcom_boots_ci(dat_aug, n_bootstrap = 200, type = "NIE", measure = "RD"))
#(res2 <- gcom_boots_ci_m_out_n(data = data_aug, n_bootstrap = 200, type = "NIE", measure = "RD"))


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
  output_file <- paste0("/n/home09/c55jiang/mediation/res/MAR11HAL_v3", type, "_n_", n_obs, "rate", censor_rate, ".txt")
  data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = censor_rate)
  simdata <- data$study_dat
  LOD_value <- as.numeric(data$lod)
  
  if (type == "Trad") {
    data_aug <- NaiveImpute(simdata, LOD_value)
  } else if (type == "EM-MLE") {
    em_result <- EM_algorithm(data = simdata,  beta_m_0 = beta_m, S=S, LOD = LOD_value)
    data_aug <- em_result$data_weted
    
  } else if (type == "Oracle") {
    data_aug <- data$study_dat_full
  }
  
  # G-computation
  nde_gcom_glm <- gcomp.nde(data_aug)
  CI.nderes <- gcom_boots_ci(data_aug, n_bootstrap = 1000, type = "NDE", measure = "RD")
  CI.nde <- CI.nderes$ci
  var.nde <- CI.nderes$variance
  nie_gcom_glm <- gcomp.nie(data_aug)
  CI.nieres <- gcom_boots_ci(data_aug, n_bootstrap = 1000, type = "NIE", measure = "RD")
  CI.nie <- CI.nieres$ci
  var.nie <- CI.nieres$variance
  
  # Save results
  results <- c(i,nde_gcom_glm$nde.rd, CI.nde, var.nde, 
               nie_gcom_glm$nie.rd, CI.nie, var.nie)
  
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

