library(data.table) # Assuming 'data' is a data.table
library(dplyr)
library(sl3)
library(hal9001)
library(medoutcon)
library(truncnorm)
library(origami)
library(future.apply)
n_obs <- 5000

(data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = 0.3))
(LOD <- data$lod)

head(data$study_dat)
head(data$study_dat_full)

dataset <- data$study_dat
dataset



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


#############
# METHOD 1
mod_update_hal <- function(data) {
  # Check and assign w_norm if NULL (for the initial value of the weights)
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  
  # Function to fit HAL model and extract details
  
  fit_hal_model <- function(X, Y, weights) {
    mod <- fit_hal(X = X, Y = Y, family = "binomial", weights = weights)
    summary_mod <- summary(mod)$table
    list(
      model = mod,
      coef = summary_mod$coef,
      terms = summary_mod$term,
      prob = predict(mod, new_data = X, type = "response")
    )
  }
  # separatly do the cv
  # Fit censoring model
  cen_res <- fit_hal_model(
    X = data[, c("L1", "L2", "L3", "A", "M")], 
    Y = data$censored, 
    weights = data$w_norm
  )
  
  # Fit outcome model
  out_res <- fit_hal_model(
    X = data[, c("L1", "L2", "L3", "A", "M")], 
    Y = data$Y, 
    weights = data$w_norm
  )
  
  # Return predictions, coefficients, and terms
  return(list(
    cen_prob = cen_res$prob, 
    out_prob = out_res$prob,
    cen_coef = cen_res$coef,
    out_coef = out_res$coef,
    cen_base = cen_res$terms,
    out_base = out_res$terms
  ))
}

system.time({ mod_pred <- mod_update_hal(data_imputed) })
mod_pred$cen_coef


mod_update_halcv <- function(data, lambda_seq = exp(seq(-0.5, -20, length = 1000)), 
                             basis_list = NULL) {
  # Check and assign w_norm if NULL (for the initial value of the weights)
  if (is.null(data$w_norm)) { data$w_norm <- 1 }
  
  # Joint cross-validated HAL fitting for censoring and outcome models
  fit_cv_hal_joint <- function(data, X, Y_cen, Y_out, weights, lambda_seq, basis_list) {
    # Define the HAL model for both outcomes
    fit_cv_hal <- function(lambda_seq, basis_list) {
      # Fit HAL for censoring model
      cen_mod <- fit_hal(
        X = X, Y = Y_cen, family = "binomial", 
        lambda = lambda_seq, weights = weights, 
        max_degree = 2, reduce_basis = 1 / sqrt(nrow(X)), smoothness_orders = 0,
        basis_list = basis_list
      )
      # Fit HAL for outcome model
      out_mod <- fit_hal(
        X = X, Y = Y_out, family = "binomial", 
        lambda = lambda_seq, weights = weights, 
        max_degree = 2, reduce_basis = 1 / sqrt(nrow(X)), smoothness_orders = 0,
        basis_list = basis_list
      )
      # Extract negative log-likelihoods for both models
      cen_nll <- -sum(log(predict(cen_mod, new_data = X, type = "response")[Y_cen == 1])) -
        sum(log(1 - predict(cen_mod, new_data = X, type = "response")[Y_cen == 0]))
      out_nll <- -sum(log(predict(out_mod, new_data = X, type = "response")[Y_out == 1])) -
        sum(log(1 - predict(out_mod, new_data = X, type = "response")[Y_out == 0]))
      list(
        cen_mod = cen_mod,
        out_mod = out_mod,
        joint_loss = cen_nll + out_nll,  # Joint negative log-likelihood
        cen_coef = coef(cen_mod),
        out_coef = coef(out_mod),
        cen_prob = predict(cen_mod, new_data = X, type = "response"),
        out_prob = predict(out_mod, new_data = X, type = "response"),
        cen_basis_list = cen_mod$basis_list,
        out_basis_list = out_mod$basis_list
      )
    }
    
    # Perform initial HAL fit to discover basis functions
    if (is.null(basis_list)) {
      initial_fit <- fit_cv_hal(lambda_seq = lambda_seq, basis_list = NULL)
      basis_list <- initial_fit$cen_basis_list  # Use censoring model's basis list as reference
    }
    
    # Evaluate models across lambda sequence using shared basis_list
    results <- lapply(lambda_seq, function(lambda) fit_cv_hal(lambda_seq = lambda_seq, basis_list = basis_list))
    joint_losses <- sapply(results, function(res) res$joint_loss)
    best_lambda_idx <- which.min(joint_losses)  # Find best lambda index
    best_result <- results[[best_lambda_idx]]  # Get best result
    best_lambda <- lambda_seq[best_lambda_idx]  # Get best lambda value
    
    return(list(
      cen_mod = best_result$cen_mod,
      out_mod = best_result$out_mod,
      cen_prob = best_result$cen_prob,
      out_prob = best_result$out_prob,
      cen_coef = best_result$cen_coef,
      out_coef = best_result$out_coef,
      best_lambda = best_lambda,
      basis_list = basis_list  # Return shared basis list
    ))
  }
  
  # Define covariates and outcomes
  X <- as.matrix(data[, c("L1", "L2", "L3", "A", "M")])
  Y_cen <- data$censored
  Y_out <- data$Y
  weights <- data$w_norm
  
  # Perform joint cross-validated HAL fitting
  joint_result <- fit_cv_hal_joint(data, X, Y_cen, Y_out, weights, lambda_seq, basis_list)
  
  # Return predictions, coefficients, and optimal lambda
  return(list(
    cen_prob = joint_result$cen_prob,
    out_prob = joint_result$out_prob,
    cen_coef = joint_result$cen_coef,
    out_coef = joint_result$out_coef,
    best_lambda = joint_result$best_lambda,
    basis_list = joint_result$basis_list
  ))
}






# Run the updated HAL fitting function
system.time({ cv_halmod_pred <- mod_update_halcv(data_imputed)})
# Access results
cv_halmod_pred$cen_prob     # Predicted probabilities for the censoring model
cv_halmod_pred$out_prob     # Predicted probabilities for the outcome model
cv_halmod_pred$best_lambda  # Optimal lambda selected
cv_halmod_pred$basis_list   # Shared basis functions


#############
# METHOD 2

data = data_imputed
fold = 1
k_folds = 10
lambda_seq1 = exp(seq(-0.5, -20, length = 20))
lambda_seq2 = exp(seq(-0.5, -20, length = 20))
basis_list = NULL


cv_hal_joint <- function(fold, data_in, x_names, y_cen_names, y_out_names, weights_names,
                         lambda_seq1 = exp(seq(-0.5, -20, length = 1e4)), 
                         lambda_seq2 = exp(seq(-0.5, -20, length = 1e4)), 
                         basis_list = NULL) {
  
  ## 0) cross-validation via origami
  train_data <- origami::training(data_in)
  valid_data <- origami::validation(data_in)
  x_train <- as.matrix(train_data[, ..x_names])
  y_cen_train <- as.numeric(train_data[, get(y_cen_names)])
  y_out_train <- as.numeric(train_data[, get(y_out_names)])
  x_valid <- as.matrix(valid_data[, ..x_names])
  y_cen_valid <- as.numeric(valid_data[, get(y_cen_names)])
  y_out_valid  <- as.numeric(valid_data[, get(y_out_names)])
  ## for weights
  weights_train <- as.numeric(train_data[, get(weights_names)])
  weights_valid <- as.numeric(valid_data[, get(weights_names)])
  
  
  ## 1) fit HAL over sequence of lambda values
  cen_mod <- hal9001::fit_hal(
    X = x_train, Y = y_cen_train, weights = weights_train, family = "binomial", 
    lambda = lambda_seq1,
    max_degree = NULL,
    return_x_basis = TRUE,
    basis_list = basis_list,fit_control = list(cv_select = FALSE)
  )
  
  out_mod <- hal9001::fit_hal(
    X = x_train, Y = y_out_train, weights = weights_train, family = "binomial", 
    lambda = lambda_seq2,
    max_degree = 2, reduce_basis = 2 / sqrt(nrow(X_train)), smoothness_orders = 0,
    basis_list = basis_list, fit_control = list(cv_select = FALSE),
    return_x_basis = TRUE
  )
  
  # Compute negative log-likelihood on validation set
  
  ## 2) predictions on validation data for each value of lambda
  valid_x_basis <- hal9001::make_design_matrix(x_valid, cen_mod$basis_list)
  valid_x_basis_clean <- valid_x_basis[, as.numeric(names(cen_mod$copy_map))]
  pred_mat <- cbind(rep(1, nrow(x_valid)), valid_x_basis_clean)
  
  ?? preds_valid <- as.matrix(pred_mat %*% coef_mat)
  
  
  cen_prob_val <- predict(cen_mod, new_data = X_val, type = "response")
  out_prob_val <- predict(out_mod, new_data = X_val, type = "response")
  
  return(list(
    joint_loss_val = cen_nll_val + out_nll_val,
    cen_mod = cen_mod,
    out_mod = out_mod,
    cen_coef = coef(cen_mod),
    out_coef = coef(out_mod)
  ))
}


fit_cv_hal_joint <- function(data_in, folds, 
                             x_names, y_cen_names, y_out_names, weights_names, 
                             lambda_seq1, lambda_seq2) {
  # fit enumerate_basis on full data to get set of basis functions
  if (is.null(basis_list)) {
    basis_list <- hal9001::enumerate_basis(as.matrix(data_in[, ..x_names]), max_degree = 3, smoothness_orders = 0)
  }
  
  # fit a cross-validated joint HAL for the Censoring and Outcome Models
  # NOTE: use set of basis functions discovered in initial enumerate_basis
  cv_fit_hal_joint <- origami::cross_validate(
    cv_fun = cv_hal_joint,
    folds = folds[-1],
    data = data_in,
    x_names, y_cen_names, y_out_names, weights_names,
    lambda_seq1 = lambda_seq1, 
    lambda_seq2 = lambda_seq2, 
    basis_list = basis_list,
    use_future = FALSE,
    .combine = FALSE
  )
  
  # Extract validation indices from folds
  idx_folds <- do.call(c, lapply(folds, `[[`, "validation_set"))
  
  # Combine predictions across folds for both Cen and Out model
  cen_cv <- do.call(rbind, c(list(gn_cv_fit_hal_initial$hat_valid),
                             cv_fit_hal_joint$hat_valid_cen))[idx_folds, ]
  
  # Combine predictions across folds for outcome model
  out_cv <- do.call(rbind, c(list(qn_cv_fit_hal_initial$hat_valid),
                             cv_fit_hal_joint$hat_valid_out))[idx_folds, ]
  
  # Define negative log-likelihood function
  nll <- function(obs, probs) {
    -sum(obs * log(probs) + (1 - obs) * log(1 - probs))
  }
  
  # Compute CV-NLL for each lammda1 for censoring model
  cv_nloglik_cen <- apply(cen_cv, 2, function(cen_hat) {
    nll(obs = data_in[, get(y_cen_names)], probs = cen_hat)
  })
  
  # Compute CV-NLL for each lammda2 for outcome model
  cv_nloglik_out <- apply(out_cv, 2, function(out_hat) {
    nll(obs = data_in[, get(y_out_names)], probs = out_hat)
  })
  
  # Create a grid of all (lambda_seq1, lambda_seq2) combinations
  cv_nloglik_sum <- outer(cv_nloglik_cen, cv_nloglik_out, "+")  # Compute sum for each pair
  
  # Find the indices of the best combination (minimizing the sum)
  best_lambda_idx <- which(cv_nloglik_sum == min(cv_nloglik_sum), arr.ind = TRUE)
  best_lambda1 <- best_lambda_idx[1, 1]  # Best lambda1 index for cen model
  best_lambda2 <- best_lambda_idx[1, 2]  # Best lambda2 index for out model
  
  return(list(best_lambda = c(best_lambda1, best_lambda2)))
}


mod_update_halcv <- function(data, lambda_seq1 = exp(seq(-0.5, -20, length = 20)), 
                             lambda_seq2 = exp(seq(-0.5, -20, length = 20)), 
                             basis_list = NULL, k_folds = 5, plot_diagnostics = TRUE) {
  
  if (is.null(data$w_norm)) { data$w_norm <- 1 }
  
  X <- as.matrix(data[, c("L1", "L2", "L3", "A", "M")])
  Y_cen <- data$censored
  Y_out <- data$Y
  weights <- data$w_norm
  
  cv_result <- fit_cv_hal_joint(data = data_in, folds, 
                                x_names, y_cen_names, y_out_names, weights_names, 
                                lambda_seq1, lambda_seq2)
  best_lambda <- cv_result$best_lambda
  # Plot diagnostics
  if (plot_diagnostics) {
    library(fields)
    image.plot(log(lambda_seq1), log(lambda_seq2), t(cv_result$mean_joint_losses),
               xlab = "log(Lambda1)", ylab = "log(Lambda2)", main = "Validation Loss (Joint Negative Log-Likelihood)")
    points(log(best_lambda[1]), log(best_lambda[2]), col = "red", pch = 19, cex = 1.5)
  }
  
 # Update HAL models with best lambda
  cen_mod <- fit_hal(
    X = as.matrix(data_in[, ..x_names]),
    Y = as.numeric(data_in[, get(y_cen_names)]),
    weights = as.numeric(data_in[, get(weights_names)]), 
    family = "binomial", 
    lambda = best_lambda[1], 
    # max_degree = 3, reduce_basis = 1 / sqrt(nrow(X)), smoothness_orders = 0,
    basis_list = basis_list, fit_control = list(cv_select = FALSE)
  )
  
  out_mod <- fit_hal(
    X = as.matrix(data_in[, ..x_names]),
    Y = as.numeric(data_in[, get(y_out_names)]),
    weights = as.numeric(data_in[, get(weights_names)]), 
    family = "binomial", 
    lambda = best_lambda[1], 
    # max_degree = 3, reduce_basis = 1 / sqrt(nrow(X)), smoothness_orders = 0,
    basis_list = basis_list, fit_control = list(cv_select = FALSE)
  )
  
  return(list(
    cen_mod = cen_mod,
    out_mod = out_mod,
    cen_coef = coef(cen_mod),
    out_coef = coef(out_mod),
    basis_list = cen_mod$basis_list
  ))
}

# Run the updated HAL fitting function
system.time({ cv_halmod_pred <- mod_update_halcv(data_imputed)})
# Access results
cv_halmod_pred$cen_prob     # Predicted probabilities for the censoring model
cv_halmod_pred$out_prob     # Predicted probabilities for the outcome model
cv_halmod_pred$best_lambda  # Optimal lambda selected
cv_halmod_pred$basis_list   # Shared basis functions


fit_cv_hal_joint <- function(data, X, Y_cen, Y_out, weights, lambda_seq1, lambda_seq2, basis_list, k_folds) {
  if (is.null(basis_list)) {
    basis_list <- hal9001::enumerate_basis(X, max_degree = 3, smoothness_orders = 0)
  }
  
  # K-fold Cross-Validation
  folds <- sample(rep(1:k_folds, length.out = nrow(X)))
  mean_joint_losses <- matrix(NA, length(lambda_seq1), length(lambda_seq2))
  
  for (fold in 1:k_folds) {
    train_idx <- which(folds != fold)
    val_idx <- which(folds == fold)
    
    X_train <- X[train_idx, ]
    Y_cen_train <- Y_cen[train_idx]
    Y_out_train <- Y_out[train_idx]
    weights_train <- weights[train_idx]
    
    X_val <- X[val_idx, ]
    Y_cen_val <- Y_cen[val_idx]
    Y_out_val <- Y_out[val_idx]
    
    for (i in seq_along(lambda_seq1)) {
      for (j in seq_along(lambda_seq2)) {
        loss <- fit_cv_hal(lambda_seq1[i], lambda_seq2[j], X_train, Y_cen_train, Y_out_train, weights_train, basis_list, X_val, Y_cen_val, Y_out_val)
        mean_joint_losses[i, j] <- ifelse(is.na(mean_joint_losses[i, j]), loss, mean_joint_losses[i, j] + loss)
      }
    }
  }
  
  mean_joint_losses <- mean_joint_losses / k_folds
  best_lambda_idx <- which(mean_joint_losses == min(mean_joint_losses), arr.ind = TRUE)
  
  best_lambda1 <- lambda_seq1[best_lambda_idx[1]]
  best_lambda2 <- lambda_seq2[best_lambda_idx[2]]
  
  return(list(
    best_lambda = c(best_lambda1, best_lambda2),
    lambda_seq1 = lambda_seq1,
    lambda_seq2 = lambda_seq2,
    mean_joint_losses = mean_joint_losses
  ))
}





# mod_pred are from mod_update_hal funciton, which inlucdes both out_prob and cen_prob
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
  data_imputed$w <- data_imputed$cen_prob * data_imputed$out_prob * 
    c(noncen, cenwet_new / cenwet_int)
  
  # Normalize weights by ID
  data_imputed$w_norm <- ave(data_imputed$w, data_imputed$id, FUN = function(x) x / sum(x))
  
  # Remove temporary columns
  data_imputed <- data_imputed[, !names(data_imputed) %in% c("cen_prob", "out_prob"), with = FALSE]
  
  return(data_imputed)
}

# Apply the updated function
data_weted <- Wet_step(data_imputed, cv_halmod_pred, beta_m_new = beta_m, beta_m_0 = beta_m, LOD = LOD_value)






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
  mod_pred <- mod_update_halcv(data_imputed)
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
    mod_pred <- mod_update_halcv(data_weted, basis_list = mod_pred$basis_list)
    
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
dat_aug <- EM_algorithm(dataset, beta_m)

dat_aug <- dat_aug$data_weted
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


nde_gcom <- gcomp.nde(dat_aug)
nde_gcom



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
  nid.rd <- mean(y1 - y0); nid.rr <- mean(y1)/mean(y0); nid.or <- (mean(y1)/(1-mean(y1))) / (mean(y0)/(1-mean(y0)))
  # Return the estimate of the indirect effect
  return(list(nid.rd = nid.rd, nid.rr = nid.rr, nid.or = nid.or))
}
nie_gcom <- gcomp.nie(dat_aug)
nie_gcom

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
  g_learners = lasso_lrnr,  #PS. interations hurts using super leaner
  h_learners = lasso_lrnr,  #xgboots 
  b_learners = lasso_lrnr,  #Outcome Model
  effect = "direct",
  estimator = "onestep",
  estimator_args = list(cv_folds = 5)
)
summary(nde_onestep)


# compute tmle estimate of the natural direct effect
nde_tmle <- medoutcon(
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
  estimator = "tmle",
  estimator_args = list(cv_folds = 5)
)
summary(nde_onestep)
summary(nde_tmle)


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


nie_tmle <- medoutcon(
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
  estimator = "tmle",
  estimator_args = list(cv_folds = 5)
)

nde_gcom
nie_gcom
summary(nde_onestep)
summary(nde_tmle)
summary(nie_onestep)
summary(nie_tmle)
