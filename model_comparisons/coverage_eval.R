# coverage_eval.R
# ----------------------------------------------------------
# Evaluates 95% credible interval coverage for factor loadings
# for both Baseline (MGPS) and TEB-FAR models
# ----------------------------------------------------------

library(rstan)

# --- USER SETTINGS ---
scenario <- 1
n_train <- 1000

sim_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations"
tebfar_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/teb_far"
baseline_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model"

sim_file <- sprintf("%s/sim_scen%d_%d.rds", sim_dir, scenario, n_train)
sim <- readRDS(sim_file)
Lambda_true <- sim$Lambda

# Function to compute coverage
coverage_rate <- function(Lambda_post, Lambda_true) {
   n_samples <- dim(Lambda_post)[1]
   n_rows <- dim(Lambda_post)[2]
   n_factors <- dim(Lambda_post)[3]
   coverage_matrix <- matrix(NA, nrow = n_rows, ncol = n_factors)
   for (i in 1:n_rows) {
      for (j in 1:n_factors) {
         samples <- Lambda_post[, i, j]
         ci <- quantile(samples, probs = c(0.025, 0.975))
         coverage_matrix[i, j] <- (Lambda_true[i, j] >= ci[1]) & (Lambda_true[i, j] <= ci[2])
      }
   }
   mean(coverage_matrix)
}

# ------------------------------
# --- Evaluate TEB-FAR Model ---
# ------------------------------
setwd(tebfar_dir)
tebfar_fit_file <- sprintf("%s/stan_tebfar_fit_scen%d.rds", sim_dir, scenario)
if (file.exists(tebfar_fit_file)) {
   tebfar_fit <- readRDS(tebfar_fit_file)
   tebfar_post <- rstan::extract(tebfar_fit)
   # If Lambda is [samples, rows, factors] in your Stan output
   # Sometimes it is [samples, rows, factors]; sometimes [samples, factors, rows]
   # If wrong shape, use aperm() to rearrange
   if (length(dim(tebfar_post$Lambda)) == 3 && dim(tebfar_post$Lambda)[2] == dim(Lambda_true)[1]) {
      Lambda_tebfar_post <- tebfar_post$Lambda
   } else if (length(dim(tebfar_post$Lambda)) == 3 && dim(tebfar_post$Lambda)[3] == dim(Lambda_true)[1]) {
      Lambda_tebfar_post <- aperm(tebfar_post$Lambda, c(1, 3, 2))
   } else {
      stop("Unexpected Lambda dimensions for TEB-FAR!")
   }
   tebfar_coverage <- coverage_rate(Lambda_tebfar_post, Lambda_true)
   cat(sprintf("TEB-FAR Lambda 95%% coverage: %.2f%%\n", 100 * tebfar_coverage))
} else {
   cat("TEB-FAR fit file not found.\n")
}

# --------------------------------------
# --- Evaluate Baseline Factor Model  ---
# --------------------------------------
setwd(baseline_dir)
baseline_fit_file <- sprintf("%s/mgps_fit_scen%d.rds", sim_dir, scenario)
if (file.exists(baseline_fit_file)) {
   baseline_fit <- readRDS(baseline_fit_file)
   baseline_post <- rstan::extract(baseline_fit)
   # Adjust for dimension order if needed
   if (length(dim(baseline_post$Lambda)) == 3 && dim(baseline_post$Lambda)[2] == dim(Lambda_true)[1]) {
      Lambda_baseline_post <- baseline_post$Lambda
   } else if (length(dim(baseline_post$Lambda)) == 3 && dim(baseline_post$Lambda)[3] == dim(Lambda_true)[1]) {
      Lambda_baseline_post <- aperm(baseline_post$Lambda, c(1, 3, 2))
   } else {
      stop("Unexpected Lambda dimensions for Baseline!")
   }
   baseline_coverage <- coverage_rate(Lambda_baseline_post, Lambda_true)
   cat(sprintf("Baseline Lambda 95%% coverage: %.2f%%\n", 100 * baseline_coverage))
} else {
   cat("Baseline fit file not found.\n")
}

