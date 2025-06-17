# model_comparisons/post_process_results.R

# -----------------------------
# Comprehensive Postprocessing Script
# For Bayesian Factor Model Simulation Studies
# -----------------------------

library(tidyverse)
library(rstan)
library(bayesplot)
library(pROC)

# --- User Settings ---

scenario <- 1
sim_dir <- "../simulations"
fit_dir <- sim_dir  # or customize if fits are elsewhere

# ---- 1. Load Simulation Data ----
sim_file <- sprintf("%s/sim_scen%d_1000.rds", sim_dir, scenario)
sim <- readRDS(sim_file)
Y <- sim$Y      # n x (p+1)
true_Lambda <- sim$Lambda
true_Omega <- sim$Omega

# --- 2. Load Model Fits (Update paths/names as needed) ---
tebfar_fit <- readRDS(file.path(fit_dir, sprintf("stan_tebfar_fit_scen%d_5.rds", scenario))) #UPDATE
baseline_fit <- readRDS(file.path(fit_dir, sprintf("mgps_fit_scen%d_5.rds", scenario))) #UPDATE
# Add others as needed

# --- 3. Predictive Metrics: Test MSE, AUC (if binary y) ---
# Assume you have y_test, y_pred_tebfar, y_pred_baseline, etc., saved in previous runs
# If not, compute them here using your workflow

# (Example) Suppose you saved these as .rds (else adapt this section)
results <- readRDS(sprintf("%s/model_comparison_summary_scen%d.rds", fit_dir, scenario))
# results should be a data.frame with columns: Model, y_test, y_pred, Test_MSE

# Plot MSE by model
ggplot(results, aes(x=Model, y=Test_MSE, fill=Model)) +
   geom_col(width=0.6) +
   labs(title="Test-set MSE by Model", y="Test MSE", x="Model") +
   theme_minimal()

# --- 4. Parameter Recovery & Interval Coverage ---

# Helper: compute 95% posterior intervals and check coverage
coverage_eval <- function(fit, true_Lambda, param_prefix="Lambda") {
   est <- rstan::extract(fit, pars=param_prefix)[[1]] # dim: draws x row x col
   coverage <- matrix(NA, nrow=nrow(true_Lambda), ncol=ncol(true_Lambda))
   for (i in 1:nrow(true_Lambda)) {
      for (j in 1:ncol(true_Lambda)) {
         draws <- est[, i, j]
         ci <- quantile(draws, c(0.025, 0.975))
         coverage[i, j] <- (true_Lambda[i, j] >= ci[1]) & (true_Lambda[i, j] <= ci[2])
      }
   }
   prop_covered <- mean(coverage, na.rm=TRUE)
   list(coverage=coverage, prop_covered=prop_covered)
}

# TEB-FAR coverage
tebfar_cov <- coverage_eval(tebfar_fit, true_Lambda, "Lambda")
cat("TEB-FAR Lambda 95% interval coverage:", tebfar_cov$prop_covered, "\n")
# Baseline coverage
baseline_cov <- coverage_eval(baseline_fit, true_Lambda, "Lambda")
cat("Baseline Lambda 95% interval coverage:", baseline_cov$prop_covered, "\n")

# --- 5. Trace, Rhat, and ESS Diagnostics ---
# Stan's summary already reports Rhat and ESS, but visualize with bayesplot
params_to_diag <- c("psi[1]", "psi[2]", "tau[1]", "tau[2]") # adjust as needed
mcmc_trace(as.array(tebfar_fit), pars=params_to_diag) + ggtitle("TEB-FAR Trace")
mcmc_rhat(rhat(tebfar_fit)) + ggtitle("TEB-FAR Rhat")
mcmc_ess_bulk(bulkESS(tebfar_fit)) + ggtitle("TEB-FAR ESS (bulk)")

# --- 6. Parameter Recovery Plots ---
# Plot true vs posterior mean for Lambda
extract_mean <- function(fit, param_prefix="Lambda") {
   est <- rstan::extract(fit, pars=param_prefix)[[1]]
   apply(est, c(2,3), mean)
}

mean_Lambda_tebfar <- extract_mean(tebfar_fit, "Lambda")
mean_Lambda_baseline <- extract_mean(baseline_fit, "Lambda")

# Combine into one long data.frame for ggplot
param_df <- tibble(
   True = as.vector(true_Lambda),
   TEBFAR = as.vector(mean_Lambda_tebfar),
   Baseline = as.vector(mean_Lambda_baseline)
)
ggplot(param_df, aes(x=True, y=TEBFAR)) +
   geom_point(alpha=0.5) +
   geom_abline(slope=1, intercept=0, color="red") +
   labs(title="Lambda Parameter Recovery (TEB-FAR)", x="True", y="Posterior Mean")

ggplot(param_df, aes(x=True, y=Baseline)) +
   geom_point(alpha=0.5) +
   geom_abline(slope=1, intercept=0, color="blue") +
   labs(title="Lambda Parameter Recovery (Baseline)", x="True", y="Posterior Mean")

# --- 7. Save All Summaries ---
summary_df <- tibble(
   Model = c("TEB-FAR", "Baseline"),
   Test_MSE = c(results$Test_MSE[results$Model=="TEB-FAR"],
                results$Test_MSE[results$Model=="Baseline"]),
   Lambda_Coverage = c(tebfar_cov$prop_covered, baseline_cov$prop_covered)
)
write_csv(summary_df, file=sprintf("%s/model_summary_scen%d.csv", fit_dir, scenario))

# --- 8. (Optional) Print Key Info for Quick Reporting ---
print(summary_df)
