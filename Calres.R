# Define output file path
output_file <- "MEDJAN29EM-MLEn_1000.txt"

# Read data
res <- read.table(output_file, header = FALSE, sep = " ", stringsAsFactors = FALSE)
ncol(res)
head(res)
# Assign appropriate column names
colnames(res) <- c("i", 
                   "nde_gcom_glm", "CI_nde_L", "CI_nde_U",
                   "nie_gcom_glm", "CI_nie_L", "CI_nie_U",
                   "nde_onestep", "CI_onestep_nde_L", "CI_onestep_nde_U",
                   "nde_tmle", "CI_tmle_nde_L", "CI_tmle_nde_U",
                   "nie_onestep", "CI_onestep_nie_L", "CI_onestep_nie_U",
                   "nie_tmle", "CI_tmle_nie_L", "CI_tmle_nie_U")

# Define true parameter values
nde_true <- 0.433
nie_true <- 0.083

# Define method columns and their corresponding true values
method_true_values <- list(
  "nde_gcom_glm" = nde_true, "nie_gcom_glm" = nie_true,
  "nde_onestep" = nde_true,  "nde_tmle" = nde_true,
  "nie_onestep" = nie_true,  "nie_tmle" = nie_true
)

# Define CI bounds for each method
ci_bounds_pairs_list <- list(
  "nde_gcom_glm" = c("CI_nde_L", "CI_nde_U"),
  "nie_gcom_glm" = c("CI_nie_L", "CI_nie_U"),
  "nde_onestep"  = c("CI_onestep_nde_L", "CI_onestep_nde_U"),
  "nde_tmle"     = c("CI_tmle_nde_L", "CI_tmle_nde_U"),
  "nie_onestep"  = c("CI_onestep_nie_L", "CI_onestep_nie_U"),
  "nie_tmle"     = c("CI_tmle_nie_L", "CI_tmle_nie_U")
)

# Function to compute bias, RMSE, and CI coverage
compute_metrics <- function(data, method_true_values, ci_bounds_pairs_list) {
  results_df <- data.frame(Method = character(), Mean_Bias = numeric(), Median_Bias = numeric(),
                           MSE = numeric(), MC_SE = numeric(), CI_Coverage = numeric(), 
                           stringsAsFactors = FALSE)
  
  for (method_name in names(method_true_values)) {
    estimates <- data[[method_name]]
    theta_true <- method_true_values[[method_name]]
    
    # Calculate bias, RMSE, and standard error
    bias <- estimates - theta_true
    mean_bias <- mean(bias, na.rm = TRUE)
    median_bias <- median(bias, na.rm = TRUE)
    mse <- mean(bias^2, na.rm = TRUE)
    mc_standard_error <- sd(estimates, na.rm = TRUE) / sqrt(sum(!is.na(estimates)))
    
    # Compute CI coverage
    ci_bounds <- ci_bounds_pairs_list[[method_name]]
    lower_bounds <- data[[ci_bounds[1]]]
    upper_bounds <- data[[ci_bounds[2]]]
    
    valid_indices <- !is.na(lower_bounds) & !is.na(upper_bounds)
    ci_covered <- sum(theta_true >= lower_bounds[valid_indices] & theta_true <= upper_bounds[valid_indices])
    ci_coverage <- ci_covered / sum(valid_indices, na.rm = TRUE)
    
    # Append results
    results_df <- rbind(results_df, data.frame(Method = method_name, Mean_Bias = mean_bias, 
                                               Median_Bias = median_bias, MSE = mse, 
                                               MC_SE = mc_standard_error, CI_Coverage = ci_coverage))
  }
  
  return(results_df)
}

# Compute and display results
results <- compute_metrics(res, method_true_values, ci_bounds_pairs_list)
print(results)
