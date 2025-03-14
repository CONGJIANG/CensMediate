# Clear the entire environment
rm(list = ls())
# Load necessary package
library(data.table)
# Observational Study with an assay lower limit (Lower Limit of Detection, LOD) for the mediator variable
#
# parameters:
#  - n_obs:  Number of study units to generate for the observational study.
#  - return_dgp: Logical indicating whether to return the data-generating
#                functions in addition to the simulated data.
#  - random: Randomization trial or not. FALSE means observational study
#
# return:
#  A data frame of `n_obs` study units from an observationl study in which 
#  an assay lower limit for the mediator variable, i.e.,
# - study_dat: A data frame of `n_obs` study units with censoring for M1.
# - study_dat_full: A data frame of `n_obs` study units without censoring for M1.

Study_dgp_Mcensor <- function(n_obs = 5000, random = FALSE, lod = NULL, censor_rate = 0.3, return_dgp = TRUE) {
  # Generate baseline covariates (e.g., demographic characteristics)
  L1 <- rbinom(n_obs, 1, 0.7)
  L2 <- rbinom(n_obs, 1, 0.5)
  L3 <- rbinom(n_obs, 1, 0.25)
  
  # Generate exposure status as a function of baseline covariates, depending on randomization
  gA <- function(L1, L2, L3, random) {
    if (random) {
      return(rep(0.5, length(L1)))  # Fixed probability if randomized
    } else {
      return(plogis(-1 + 0.5 * L1 + 1.25 * L2 + 0.75 * L3 - 1.25 * L1 * L3)) # interactions
    }
  }
  gA_obs <- gA(L1, L2, L3, random)
  A <- rbinom(n_obs, 1, gA_obs)
  # summary(A); hist(A)
  
  # Generate post-infection biomarker (RNA quantity) from a log-normal distribution
  QM <- function(A, L1, L2, L3) {
    QM <- -3 + 1.5 * A + 1.75 * L1 + 1.5 * L2 - 0.25 * L3
    return(list(QM = QM))
  }
  
  QM_obs <- QM(A, L1, L2, L3)
  M_full <- rlnorm(n_obs, QM_obs$QM, 0.25)
  summary(M_full); hist(M_full)
  # different distribution, e.g., round(rweibull(n_obs, shape = shape, scale = scale) * (1 + QM_obs$QM), 3)
  
  # Apply fixed LOD for censoring
  if (is.null(lod)) {
    lod <- quantile(M_full, probs = censor_rate)  # Set LOD as a fixed quantile of simulated data if not provided
  }
  C <- M_full >= lod  # Identify uncensored values (C = 1 if M_full >= lod)
  M <- M_full
  M[!C] <- "LOD"  # Replace values below LOD with "LOD"
  
  # Generate the outcome of interest
  QY <- function(A, M, L1, L2, L3) {
    -1 + 2.5 * A + 1.75 * M + 0.5*A*M - 2.25 * L1 - 1.75 * L2 - 1.5 * L3 # additive model lasso is ok. 
  }
  QY_obs <- QY(A, M_full, L1, L2, L3)
  Yprob <- plogis(QY_obs)  # Convert to probabilities using the logistic function
  Y <- rbinom(n_obs, 1, Yprob) 
  summary(Yprob); hist(Y)
  hist(Yprob)
  
  # Create data frame and return
  study_dat <- data.table::setDT(
    data.frame(id = seq(1, n_obs), L1, L2, L3, A, M, Y, C = C)
  )
  study_dat_full <- data.table::setDT(
    data.frame(id = seq(1, n_obs), L1, L2, L3, A, M = M_full, Y)
  )
  # Collect DGP functions to optionally return
  dgp_funs <- list(gA = gA, QM = QM, QY = QY)
  
  # Return data frame and optionally DGP functions
  if (return_dgp) {
    out <- list(
      "study_dat" = study_dat,
      "study_dat_full" = study_dat_full,
      "dgp_funs" = dgp_funs,
      "lod" = lod
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
(data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = 0.3))
data$lod

head(data$study_dat)
head(data$study_dat_full)

summary(data$study_dat$M)
summary(data$study_dat_full$M)



# Approximate Average Treatment Effect (ATE), Natural Direct and Indirect Effects (NDE/NIE)
#
# parameters:
#  - gen_data: A function that generates data from a data-generating process.
#  - n_truth: Number of observations to generate for an "asymptotic" sample
#              in which the plug-in estimator of the ATE will be used.
#  - random: Randomization trial or not. FALSE means observational study
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
  M1_Ais1 <- rlnorm(n_truth, QM_Ais1_mech$QM, 0.25)
  
  # compute post-infection labs and biomarkers in "asymptotic" sample for
  # counterfactual data where all units are unexposed (A = 0)
  QM_Ais0_mech <- with(wow_so_much_data, QM(0, L1, L2, L3))
  M1_Ais0 <- rlnorm(n_truth, QM_Ais0_mech$QM, 0.25)
  
  # Compute outcome of interest in "asymptotic" sample
  Y_mech <- with(wow_so_much_data, QY(A, M, L1, L2, L3))
  QY_mech <- rbinom(n_truth, 1, plogis(Y_mech) )
  
  Y_Ais1_mech <- with(wow_so_much_data, QY(1, M1_Ais1, L1, L2, L3))
  Y_Ais0_mech <- with(wow_so_much_data, QY(0, M1_Ais0, L1, L2, L3))
  
  QY_Ais1_mech <- rbinom(n_truth, 1, plogis(Y_Ais1_mech) )
  QY_Ais0_mech <- rbinom(n_truth, 1, plogis(Y_Ais0_mech) )
  
  
  # Compute average treatment effect (ATE) via plug-in in "asymptotic" sample
  ate_rd <- mean(QY_Ais1_mech) - mean(QY_Ais0_mech)
  ate_rr <- mean(QY_Ais1_mech)/mean(QY_Ais0_mech)
  ate_or <- (mean(QY_Ais1_mech)/(1 - mean(QY_Ais1_mech)))/(mean(QY_Ais0_mech)/(1 - mean(QY_Ais0_mech)))
  
  # Compute direct and indirect effects
  Y_corss_mech <- with(wow_so_much_data, QY(1, M1_Ais0, L1, L2, L3))
  QY_corss_mech <- rbinom(n_truth, 1, plogis(Y_corss_mech) )
  
  nie_rd <- mean(QY_Ais1_mech) - mean(QY_corss_mech)  # Natural Indirect Effect RD
  nie_rr <- mean(QY_Ais1_mech)/mean(QY_corss_mech)  # Natural Indirect Effect RR
  nie_or <- (mean(QY_Ais1_mech)/(1 - mean(QY_Ais1_mech)))/(mean(QY_corss_mech)/(1-mean(QY_corss_mech)) ) # Natural Indirect Effect OR
  
  nde_rd <- mean(QY_corss_mech) - mean(QY_Ais0_mech)  # Natural Direct Effect RD
  nde_rr <- mean(QY_corss_mech)/mean(QY_Ais0_mech)  # Natural Direct Effect RD
  nde_or <- (mean(QY_corss_mech)/(1-mean(QY_corss_mech)) )/(mean(QY_Ais0_mech)/(1 - mean(QY_Ais0_mech)))  # Natural Direct Effect RD
  
  
  # Compute the variance bound via the EIF
  # EIF for ATE
  eif_Ais1 <- (wow_so_much_data$A == 1) / gA_mech * (wow_so_much_data$Y - QY_mech) + QY_Ais1_mech
  eif_Ais0 <- (wow_so_much_data$A == 0) / (1 - gA_mech) * (wow_so_much_data$Y - QY_mech) + QY_Ais0_mech
  eif_ate <- (eif_Ais1 - eif_Ais0) - ate_rd
  
  # EIF for NIE and NDE
  return(list(
    psi_rd = ate_rd,
    psi_rr = ate_rr,
    psi_or = ate_rr,
    var_eif = var(eif_ate),
    nde_rd = nde_rd,
    nde_rr = nde_rr,
    nde_or = nde_or,
    nie_rd = nie_rd,
    nie_rr = nie_rr,
    nie_or = nie_or
  ))
}

#(truth <- get_truth(gen_data = Study_dgp_Mcensor))


# $nde_rd [1] 0.4131817; $nie_rd [1] 0.3686884
# $nde_rd [1] 0.4135181  $nie_rd [1] 0.3683857
# $nde_rd [1]  0.4131935  $nie_rd [1] 0.3689446
