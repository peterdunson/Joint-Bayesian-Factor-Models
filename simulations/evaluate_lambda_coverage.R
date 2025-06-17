# evaluate_lambda_coverage.R
# --------------------------------------
# Evaluate 95% CI coverage and CI symmetry for Lambda (factor loadings)
# For both TEB-FAR and Baseline Stan fits
# --------------------------------------

library(rstan)

# -------------- USER SETTINGS --------------
scenario <- 1
n_train <- 1000

sim_dir      <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations"
tebfar_dir   <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/teb_far"
baseline_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model"

sim_file           <- sprintf("%s/sim_scen%d_%d.rds", sim_dir, scenario, n_train)
tebfar_fit_file    <- sprintf("%s/stan_tebfar_fit_scen%d.rds", sim_dir, scenario)
baseline_fit_file  <- sprintf("%s/mgps_fit_scen%d.rds", sim_dir, scenario)

# -------------- Load Data & True Parameters --------------
sim         <- readRDS(sim_file)
Lambda_true <- sim$Lambda  # [p+1, K] or [p, K], check your simulation code

# -------------- Helper Function --------------
evaluate_ci <- function(posterior_array, Lambda_true, method_name = "TEB-FAR") {
   # posterior_array: [iterations, rows, cols]
   # Lambda_true: [rows, cols]
   lower <- apply(posterior_array, c(2, 3), quantile, probs = 0.025)
   upper <- apply(posterior_array, c(2, 3), quantile, probs = 0.975)
   mean_est <- apply(posterior_array, c(2, 3), mean)
   
   coverage <- (Lambda_true >= lower) & (Lambda_true <= upper)
   coverage_rate <- mean(coverage)
   
   # CI symmetry: (upper - mean) vs (mean - lower)
   upper_diff <- upper - mean_est
   lower_diff <- mean_est - lower
   symmetry <- abs(upper_diff - lower_diff) / (upper_diff + lower_diff + 1e-10)  # Avoid 0/0
   
   cat(sprintf(
      "%s Lambda 95%% coverage: %.2f%%\n", method_name, 100 * coverage_rate
   ))
   
   cat(sprintf(
      "%s average CI symmetry (0 = perfect): %.3f\n", method_name, mean(symmetry)
   ))
   
   hist(symmetry, breaks = 30, main = sprintf("%s Lambda CI Symmetry", method_name), xlab = "Symmetry measure")
   invisible(list(
      coverage = coverage,
      coverage_rate = coverage_rate,
      symmetry = symmetry
   ))
}

# -------------- TEB-FAR --------------------
if (file.exists(tebfar_fit_file)) {
   tebfar_fit <- readRDS(tebfar_fit_file)
   tebfar_post <- rstan::extract(tebfar_fit)
   # posterior$Lambda: [iterations, p+1, K]
   # Make sure Lambda_true and posterior$Lambda match in shape
   res_tebfar <- evaluate_ci(tebfar_post$Lambda, Lambda_true, method_name = "TEB-FAR")
}

# -------------- Baseline -------------------
if (file.exists(baseline_fit_file)) {
   baseline_fit <- readRDS(baseline_fit_file)
   baseline_post <- rstan::extract(baseline_fit)
   res_baseline <- evaluate_ci(baseline_post$Lambda, Lambda_true, method_name = "Baseline")
}
