library(data.table) # Assuming 'data' is a data.table
n_obs <- 50
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
  return(rbind(data[data$M != "LOD", ], do.call(rbind, imputed_rows_list)))
}

(data_imputed <- Impute_step(data = dataset, beta_m, S = 10))





Wet_step <- function(data_imputed, beta_m_new, beta_m_0) {
  data_imputed$w_norm <- 1
  # Censoring Modeling (non parametric assumption, include Y? )
  cen_mod <- glm(censored ~ L1 + L2 + L3 + A + M, data = data_imputed, family = binomial, weights = w_norm)
  cen_prob <- predict(cen_mod, type = "response")
  # Outcome Modeling
  out_mod <- glm(Y ~ L1 + L2 + L3 + A + M, data = data_imputed, family = binomial, weights = w_norm)
  out_prob <- predict(out_mod, type = "response")
  
  data_imputed$cen_prob <- cen_prob; data_imputed$out_prob <- out_prob;
  censored_data <- data_imputed[data_imputed$censored == TRUE, ]
  cenwet_new <- compute_f_M(censored_data, LOD = LOD_value, beta_m_new)
  cenwet_int <- compute_f_M(censored_data, LOD = LOD_value, beta_m_0)
  noncen <- rep(1, nrow(data_imputed) - nrow(censored_data))
  data_imputed$w <- cen_prob * out_prob * c(noncen, cenwet_new/cenwet_int)
  # Normalize weights
  data_imputed$w_norm <- ave(data_imputed$w, data_imputed$id, FUN = function(x) x / sum(x))
  data_imputed <- data_imputed[, !names(data_imputed) %in% c("cen_prob", "out_prob"), with = FALSE]
  return(data_imputed) 
}

(data_weted <- Wet_step(data_imputed, beta_m_new = beta_m, beta_m_0 = beta_m))





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
EM_algorithm <- function(data, beta_m_0, S = 5, LOD = LOD_value, max_iter = 1000, tol = 1e-6) {
  iter <- 0
  diff <- Inf
  # I-Step
  data_imputed <- Impute_step(data, beta_m_0, S)
  while (iter < max_iter && diff > tol) {
    # W-Step
    data_weted <- Wet_step(data_imputed, beta_m_new = beta_m, beta_m_0 = beta_m_0)
    
    # Maximization Step
    beta_m_new <- M_stepf_M(data = data_weted, LOD = LOD_value, beta_m)
    
    # Calculate convergence criteria
    diff <- sum((beta_m_new - beta_m)^2)
    beta_m <- beta_m_new
    iter <- iter + 1
    cat("Iteration:", iter, "Difference:", diff, "\n")
  }
  
  return(list(beta_m = beta_m, data_weted = data_weted))
}


# Define initial guesses for beta_m parameters
beta_m <- c(A = 0.5, L1 = 0.3, L2 = 0.2, L3 = 0.1, sd = 1)
LOD_value <- 1.0  # Set the LOD value
# Run the EM algorithm
EM_algorithm(dataset, beta_m)

