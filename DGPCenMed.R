# Load necessary package
library(data.table)
# Observational Study with an assay lower limit (Lower Limit of Detection, LOD) for the mediator variable
#
# parameters:
#  - n_obs:  Number of study units to generate for the observational study.
#  - return_dgp: Logical indicating whether to return the data-generating
#                functions in addition to the simulated data.
#
# return:
#  A data frame of `n_obs` study units from an observationl study in which 
#  an assay lower limit for the mediator variable, i.e.,
# - study_dat: A data frame of `n_obs` study units with censoring for M1.
# - study_dat_full: A data frame of `n_obs` study units without censoring for M1.

Study_dgp_Mcensor <- function(n_obs = 5000, random = FALSE, shape = 3, scale = 1, lod = NULL, censor_rate = 0.7, return_dgp = FALSE) {
  # Generate baseline covariates (e.g., demographic characteristics)
  L1 <- rbinom(n_obs, 1, 0.7)
  L2 <- rbinom(n_obs, 1, 0.5)
  L3 <- rbinom(n_obs, 1, 0.25)
  
  # Generate infection status as a function of baseline covariates, depending on randomization
  gA <- function(L1, L2, L3, random) {
    if (random) {
      return(rep(0.5, length(L1)))  # Fixed probability if randomized
    } else {
      return(plogis(-1.25 + 0.25 * L1 + 1.25 * L2 + 0.5 * L3 - 2.5 * L1 * L3))
    }
  }
  gA_obs <- gA(L1, L2, L3, random)
  A <- rbinom(n_obs, 1, gA_obs)
  
  # Generate post-infection biomarker (RNA quantity) from a log-normal distribution
  QM <- function(A, L1, L2, L3) {
    QM <- A * (L1 + 2 * L2 - L3)
    return(list(QM = QM))
  }
  QM_obs <- QM(A, L1, L2, L3)
  M_full <- round(rweibull(n_obs, shape = shape, scale = scale) * (1 + QM_obs$QM), 3)
  
  # Apply fixed LOD for censoring
  if (is.null(lod)) {
    lod <- quantile(M_full, probs = censor_rate)  # Set LOD as a fixed quantile of simulated data if not provided
  }
  censored <- M_full < lod                         # Identify which values are censored
  M <- ifelse(censored,"LOD" , M_full)             # Apply censoring %lod/2
  
  # Generate the outcome of interest
  QY <- function(A, M, L1, L2, L3) {
    0.5 + 2 * A + 0.5 * M + 0.25 * L1 - 0.5 * L2 -0.5 * L3
  }
  QY_obs <- QY(A, M_full, L1, L2, L3)
  Yprob <- plogis(QY_obs)  # Convert to probabilities using the logistic function
  Y <- rbinom(n_obs, 1, Yprob) 
  
  # Create data frame and return
  study_dat <- data.table::setDT(
    data.frame(L1, L2, L3, A, M, Y, censored = censored)
  )
  study_dat_full <- data.table::setDT(
    data.frame(L1, L2, L3, A, M = M_full, Y)
  )
  # Collect DGP functions to optionally return
  dgp_funs <- list(gA = gA, QM = QM, QY = QY)
  
  # Return data frame and optionally DGP functions
  if (return_dgp) {
    out <- list(
      "study_dat" = study_dat,
      "study_dat_full" = study_dat_full,
      "dgp_funs" = dgp_funs
    )
  } else {
    out <- list(
      "study_dat" = study_dat,
      "study_dat_full" = study_dat_full
    )
  }
  return(out)
}



n_obs <- 5000
(data <- Study_dgp_Mcensor(n_obs = n_obs, lod = 1.0))
sum(data$study_dat$censored)/n_obs *100

head(data$study_dat)
head(data$study_dat_full)


summary(data$study_dat$M1)
summary(data$study_dat_full$M1)



# Approximate Average Treatment Effect (ATE), Natural Direct and Indirect Effects
#
# parameters:
#  - gen_data: A function that generates data from a data-generating process.
#  - n_truth: Number of observations to generate for an "asymptotic" sample
#              in which the plug-in estimator of the ATE will be used.
#
# return:
#   A list containing:
#   - psi_est: Estimate of the ATE based on the plug-in estimator.
#   var_eif: Variance bound via the efficient influence function (EIF).
#   nde: Natural direct effect (NDE).
#   nie: Natural indirect effect (NIE).

get_truth <- function(gen_data, random = FALSE, n_truth = 1e7) {
  # Compute large data set from data-generating mechanism
  dgp <- gen_data(n_obs = n_truth, random =random, return_dgp = TRUE)
  wow_so_much_data <- dgp$study_dat_full
  
  # Extract DGP helper functions
  gA <- dgp$dgp_funs$gA
  QM <- dgp$dgp_funs$QM
  QY <- dgp$dgp_funs$QY
  
  # Compute exposure (infection) mechanism in "asymptotic" sample
  gA_mech <- with(wow_so_much_data, gA(L1, L2, L3, random))
  
  # compute post-infection labs and biomarkers in "asymptotic" sample
  QM_Anat_mech <- with(wow_so_much_data, QM(A, L1, L2, L3))
  
  # compute post-infection labs and biomarkers in "asymptotic" sample for
  # counterfactual data where all units are exposed (A = 1)
  QM_Ais1_mech <- with(wow_so_much_data, QM(1, L1, L2, L3))
  M1_Ais1 <- rnorm(n_truth, QM_Ais1_mech$QM, 1)
  
  # compute post-infection labs and biomarkers in "asymptotic" sample for
  # counterfactual data where all units are unexposed (A = 0)
  QM_Ais0_mech <- with(wow_so_much_data, QM(0, L1, L2, L3))
  M1_Ais0 <- rnorm(n_truth, QM_Ais0_mech$QM, 1)
  
  # Compute outcome of interest in "asymptotic" sample
  QY_mech <- with(wow_so_much_data, QY(A, M, L1, L2, L3))
  QY_Ais1_mech <- with(wow_so_much_data, QY(1, M1_Ais1, L1, L2, L3))
  QY_Ais0_mech <- with(wow_so_much_data, QY(0, M1_Ais0, L1, L2, L3))
  

  # Compute average treatment effect (ATE) via plug-in in "asymptotic" sample
  ate_approx <- mean(QY_Ais1_mech - QY_Ais0_mech)
  
  # Compute direct and indirect effects
  QY_corss_mech <- with(wow_so_much_data, QY(1, M1_Ais0, L1, L2, L3))
  nie <- mean(QY_Ais1_mech - QY_corss_mech)  # Natural Indirect Effect
  nde <- mean(QY_corss_mech - QY_Ais0_mech)  # Natural Direct Effect

  
  # Compute the variance bound via the EIF
  eif_Ais1 <- (wow_so_much_data$A == 1) / gA_mech * (wow_so_much_data$Y - QY_mech) + QY_Ais1_mech
  eif_Ais0 <- (wow_so_much_data$A == 0) / (1 - gA_mech) * (wow_so_much_data$Y - QY_mech) + QY_Ais0_mech
  eif_ate <- (eif_Ais1 - eif_Ais0) - ate_approx
  
  # EIF for 
  return(list(
    psi_est = ate_approx,
    var_eif = var(eif_ate),
    nde = nde,
    nie = nie
  ))
}

get_truth(Study_dgp_Mcensor)



n <- 100
Z <- rbinom(n, size = 1, prob = 0.5)

# Step 2: Generate X from Weibull distribution with shape = 0.75 and scale depending on Z
shape_X <- 0.75
scale_X <- 0.25 + 0.25 * Z
X <- rweibull(n, shape = shape_X, scale = scale_X)
hist(X)
summary(X)
# Step 3: Generate the continuous outcome Y
e <- rnorm(n, mean = 0, sd = 1)  # Standard normal error term
Y <- 1 + 0.5 * X + 0.25 * Z + e

# Step 4: Generate censoring times C from an exponential distribution with specified rate
C <- rexp(n, rate = 2)
hist(C)
# Step 5: Construct W = min(X, C) and Δ = I(X <= C)
W <- pmin(X, C)
Delta <- as.numeric(X <= C)




# Function for the EM algorithm
semiparametric_em <- function(L, A, M, Y, init_params, max_iter = 100, tol = 1e-6) {
  # L: covariates (matrix of n rows, p covariates)
  # A: treatment/exposure (vector of n)
  # M: observed mediator (vector of n)
  # Y: outcome (vector of n)
  # init_params: list of initial functions and probabilities, f^0(.), P^0(.)
  # max_iter: maximum number of iterations
  # tol: convergence tolerance
  
  n <- length(A)  # Sample size
  current_params <- init_params  # Initialize parameters
  diff <- tol + 1  # Initialize difference for convergence check
  iter <- 0  # Initialize iteration counter
  
  # E-step: function to compute posterior density
  e_step <- function(L, A, M, Y, params) {
    p_t <- numeric(n)  # Vector for posterior density
    
    for (i in 1:n) {
      num <- params$f_M_star_L(L[i,]) * params$f_A_M_star_L(A[i], params$M_star, L[i,]) *
        params$f_Y_M_star_A_L(Y[i], params$M_star, A[i], L[i,]) * 
        params$f_M_M_star_A_L(M[i], params$M_star, A[i], L[i,])
      
      denom <- integrate(function(M_star) {
        params$f_M_star_L(L[i,]) * params$f_A_M_star_L(A[i], M_star, L[i,]) *
          params$f_Y_M_star_A_L(Y[i], M_star, A[i], L[i,]) *
          params$f_M_M_star_A_L(M[i], M_star, A[i], L[i,])
      }, -Inf, Inf)$value
      
      p_t[i] <- num / denom  # Posterior density
    }
    
    return(p_t)
  }
  
  # M-step: function to update the parameters using the posterior density
  m_step <- function(p_t, L, A, M, Y) {
    new_params <- list()
    
    # Update the function f(M*)
    new_params$f_M_star <- mean(p_t)
    
    # Update f(M* | L)
    new_params$f_M_star_L <- function(L) {
      mean(p_t[L == L])
    }
    
    # Update f(M* | A, L)
    new_params$f_M_star_A_L <- function(A, L) {
      mean(p_t[A == A & L == L])
    }
    
    # Update f(M* | Y, A, L)
    new_params$f_M_star_Y_A_L <- function(Y, A, L) {
      mean(p_t[Y == Y & A == A & L == L])
    }
    
    # Update f(Y | M*, A, L)
    new_params$f_Y_M_star_A_L <- function(Y, M_star, A, L) {
      numerator <- mean(Y[M_star == M_star & A == A & L == L])
      denominator <- mean(M_star[M_star == M_star & A == A & L == L])
      numerator / denominator
    }
    
    # Update the probability P(A | M*, L)
    new_params$P_A_M_star_L <- function(A, M_star, L) {
      numerator <- mean(A[M_star == M_star & L == L])
      denominator <- mean(M_star[M_star == M_star & L == L])
      numerator / denominator
    }
    
    return(new_params)
  }
  
  # Main EM loop
  while (iter < max_iter & diff > tol) {
    # E-step
    p_t <- e_step(L, A, M, Y, current_params)
    
    # M-step
    new_params <- m_step(p_t, L, A, M, Y)
    
    # Check convergence
    diff <- sum(abs(unlist(new_params) - unlist(current_params)))
    current_params <- new_params
    iter <- iter + 1
    
    cat("Iteration:", iter, "Difference:", diff, "\n")
  }
  
  # Return final parameters
  return(current_params)
}

# Example usage:
# Assuming L, A, M, Y are the data vectors/matrix and init_params is a list with initial parameter estimates
# Replace the initial functions (init_params) with your own starting estimates

init_params <- list(
  f_M_star_L = function(L) { rnorm(1, mean = 0, sd = 1) },  # Example function
  f_A_M_star_L = function(A, M_star, L) { dnorm(A, mean = M_star, sd = 1) },
  f_Y_M_star_A_L = function(Y, M_star, A, L) { dnorm(Y, mean = M_star + A, sd = 1) },
  f_M_M_star_A_L = function(M, M_star, A, L) { dnorm(M, mean = M_star + A, sd = 1) }
)

# Run the EM algorithm
final_params <- semiparametric_em(L, A, M, Y, init_params)



set.seed(123)  # For reproducibility

# Sample size
n <- 100

# Generate data
L <- rnorm(n, mean = 50, sd = 10)  # Covariate (e.g., age)
A <- rbinom(n, size = 1, prob = 0.5)  # Treatment (e.g., smoking status)
M <- 0.5 * L + 0.3 * A + rnorm(n, mean = 0, sd = 1)  # Mediator (e.g., blood pressure)
Y <- 0.7 * M + 0.4 * A + 0.2 * L + rnorm(n, mean = 0, sd = 1)  # Outcome (e.g., health score)

# Visualize data
head(data.frame(L, A, M, Y))


# Initial parameter functions (guesses)
init_params <- list(
  f_M_star_L = function(L) { rnorm(1, mean = mean(L), sd = 1) },  # f(M* | L)
  f_A_M_star_L = function(A, M_star, L) { dnorm(A, mean = M_star, sd = 1) },  # f(A | M*, L)
  f_Y_M_star_A_L = function(Y, M_star, A, L) { dnorm(Y, mean = M_star + A, sd = 1) },  # f(Y | M*, A, L)
  f_M_M_star_A_L = function(M, M_star, A, L) { dnorm(M, mean = M_star + A, sd = 1) }  # f(M | M*, A, L)
)


# Run the EM algorithm
final_params <- semiparametric_em(L, A, M, Y, init_params)

# Display final parameter estimates
print(final_params)





