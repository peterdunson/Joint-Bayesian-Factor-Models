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




# ---- Frequentist CI Coverage ----
library(glmnet)

# Read simulation (already loaded above)
# Lambda_true assumed to be [p+1, K] or [p, K] (usually predictors x factors)

# --- OLS ---
# For OLS, we can use the standard errors from the model summary
X <- sim$Y[, -1]
y <- sim$Y[, 1]

X <- scale(X)
y <- scale(y, center = TRUE, scale = TRUE)

ols_mod <- lm(y ~ X)
ols_summary <- summary(ols_mod)
coef_ols <- coef(ols_mod)[-1]           # Remove intercept
se_ols   <- coef(ols_summary)[-1, 2]    # Remove intercept SE

ci_ols_lower <- coef_ols - 1.96 * se_ols
ci_ols_upper <- coef_ols + 1.96 * se_ols

# For coverage, need to compare to the *regression effect* (not Lambda from FA) - only meaningful if true beta is known
# If simulation includes "true_beta" (regression vector), use:
if (!is.null(sim$beta)) {
   beta_true <- sim$beta
   coverage_ols <- mean(beta_true >= ci_ols_lower & beta_true <= ci_ols_upper)
   cat(sprintf("OLS regression coefficient 95%% coverage: %.2f%%\n", 100 * coverage_ols))
}

# --- PCA Regression ---
# Get the first principal component loadings and see if CIs for loadings cover true values
pca <- prcomp(X, center = FALSE, scale. = FALSE)
# Only compare to first factor if that's how the data is simulated (as in Scenario 1)
loadings_pca1 <- pca$rotation[, 1]
# Approximate SE for PCA loadings is tricky (no easy formula; could bootstrap if desired)
# For illustration: Assume normality, CI = coef ± 1.96 * (sd(loadings) / sqrt(n))
se_pca1 <- sd(loadings_pca1) / sqrt(nrow(X))
ci_pca1_lower <- loadings_pca1 - 1.96 * se_pca1
ci_pca1_upper <- loadings_pca1 + 1.96 * se_pca1

# Compare to true loadings for the first factor
if (ncol(Lambda_true) >= 1) {
   coverage_pca1 <- mean(Lambda_true[-1, 1] >= ci_pca1_lower & Lambda_true[-1, 1] <= ci_pca1_upper)
   cat(sprintf("PCA 1st factor loading 95%% coverage: %.2f%%\n", 100 * coverage_pca1))
}

# --- Lasso & Ridge ---
# Use bootstrapping for SEs if you want actual intervals; here just use coefficients (no SEs by default)
# (Optional, if you want to include)
set.seed(1)
lasso_fit <- cv.glmnet(X, y, alpha = 1)
beta_lasso <- as.vector(coef(lasso_fit, s = "lambda.min"))[-1]
ridge_fit <- cv.glmnet(X, y, alpha = 0)
beta_ridge <- as.vector(coef(ridge_fit, s = "lambda.min"))[-1]

# For Lasso/Ridge, bootstrapping for CIs (simplified version)
B <- 100
lasso_boot <- matrix(NA, nrow = B, ncol = ncol(X))
ridge_boot <- matrix(NA, nrow = B, ncol = ncol(X))
for (b in 1:B) {
   idx <- sample(1:nrow(X), replace = TRUE)
   lfit <- cv.glmnet(X[idx, ], y[idx], alpha = 1)
   lasso_boot[b, ] <- as.vector(coef(lfit, s = "lambda.min"))[-1]
   rfit <- cv.glmnet(X[idx, ], y[idx], alpha = 0)
   ridge_boot[b, ] <- as.vector(coef(rfit, s = "lambda.min"))[-1]
}
lasso_se <- apply(lasso_boot, 2, sd)
ridge_se <- apply(ridge_boot, 2, sd)

ci_lasso_lower <- beta_lasso - 1.96 * lasso_se
ci_lasso_upper <- beta_lasso + 1.96 * lasso_se
ci_ridge_lower <- beta_ridge - 1.96 * ridge_se
ci_ridge_upper <- beta_ridge + 1.96 * ridge_se

if (!is.null(sim$beta)) {
   coverage_lasso <- mean(beta_true >= ci_lasso_lower & beta_true <= ci_lasso_upper)
   coverage_ridge <- mean(beta_true >= ci_ridge_lower & beta_true <= ci_ridge_upper)
   cat(sprintf("Lasso regression coefficient 95%% coverage: %.2f%%\n", 100 * coverage_lasso))
   cat(sprintf("Ridge regression coefficient 95%% coverage: %.2f%%\n", 100 * coverage_ridge))
}




