# ---------------------------------------------------------------
# Evaluate All Factor Models on Simulated Data
# ---------------------------------------------------------------

library(rstan)
library(ggplot2)
library(MatrixCorrelation)
library(tidyverse)
library(pheatmap)

set.seed(42)
scenario <- 1
n_train <- 1000
K <- 5
sim_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations"

# ---- File names ----
sim_file <- sprintf("%s/sim_scen1_1000.rds", sim_dir)
mgps_file <- sprintf("%s/mgps_fit_scen1_5.rds", sim_dir)
ssl_file <- sprintf("%s/ssl_factor_fit_scen1_5.rds", sim_dir)
horseshoe_file <- sprintf("%s/stan_horseshoe_fit_scen1_5.rds", sim_dir)
tebfar_file <- sprintf("%s/stan_tebfar_fit_scen1_5.rds", sim_dir)
vimsfa_cavi_file <- sprintf("%s/vimsfa_cavi_fit_scen1_5.rds", sim_dir)
vimsfa_svi_file  <- sprintf("%s/vimsfa_svi_fit_scen1_5.rds", sim_dir)

# ---- Load simulation data ----
sim <- readRDS(sim_file)
yX <- sim$Y      # n x (p+1): first column is y
Lambda_true <- sim$Lambda # (p+1) x K

X <- yX[, -1, drop=FALSE]
y <- yX[, 1]

# ---- Standardize ----
X <- scale(X)
y <- scale(y, center = TRUE, scale = TRUE)

# ---- Train/test split (75/25) ----
n <- nrow(X)
test_frac <- 0.25
test_idx <- sample(1:n, size = floor(n * test_frac))
train_idx <- setdiff(1:n, test_idx)

X_train <- X[train_idx, , drop=FALSE]
y_train <- y[train_idx]
X_test  <- X[test_idx, , drop=FALSE]
y_test  <- y[test_idx]

# ---- Predict from factors helper ----
predict_y_from_factors <- function(Xtest, Lambda, psi, y_col = 1) {
   p <- ncol(Xtest)
   K <- ncol(Lambda)
   Lambda_x <- Lambda[2:(p+1), , drop=FALSE]
   Lambda_y <- Lambda[y_col, , drop=FALSE]
   psi_x <- psi[2:(p+1)]
   psi_y <- psi[y_col]
   Psi_x_inv <- diag(1 / psi_x)
   Lt_PsiInv <- t(Lambda_x) %*% Psi_x_inv
   V_eta <- solve(diag(K) + Lt_PsiInv %*% Lambda_x)
   preds <- numeric(nrow(Xtest))
   for (i in 1:nrow(Xtest)) {
      M_eta <- V_eta %*% Lt_PsiInv %*% Xtest[i, ]
      preds[i] <- as.numeric(Lambda_y %*% M_eta)
   }
   preds
}

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

# ---- Evaluate/posterior mean for Stan models ----
evaluate_stan_fit <- function(fit_file, method_name, Lambda_true = NULL) {
   if (!file.exists(fit_file)) {
      cat("Fit file missing:", fit_file, "\n")
      return(list(mse = NA))
   }
   fit <- readRDS(fit_file)
   post <- rstan::extract(fit)
   Lambda <- apply(post$Lambda, c(2,3), mean)
   psi <- apply(post$psi, 2, mean)
   # Prediction
   pred <- predict_y_from_factors(X_test, Lambda, psi)
   mse <- mean((y_test - pred)^2)
   # Heatmap
   print(pheatmap::pheatmap(Lambda, cluster_rows = FALSE, cluster_cols = FALSE,
                            main = sprintf("%s Factor Loadings", method_name)))
   # Histogram
   hist(as.numeric(Lambda), breaks = 50, main = sprintf("%s Posterior Mean Loadings", method_name))
   # Sparsity
   sparse_frac <- mean(abs(Lambda) < 0.05)
   cat(sprintf("[%s] Test-set MSE: %.4f | Fraction near-zero loadings (<0.05): %.2f\n",
               method_name, mse, sparse_frac))
   # --- Coverage and CI symmetry ---
   if (!is.null(Lambda_true)) {
      # Try permutations if dims mismatched
      post_lambda_dim <- dim(post$Lambda)
      true_dim <- dim(Lambda_true)
      if (length(post_lambda_dim) == 3 && all(post_lambda_dim[2:3] == true_dim)) {
         evaluate_ci(post$Lambda, Lambda_true, method_name = method_name)
      } else if (length(post_lambda_dim) == 3 && post_lambda_dim[3] == true_dim[1]) {
         Lambda_perm <- aperm(post$Lambda, c(1, 3, 2))
         if (all(dim(Lambda_perm)[2:3] == true_dim)) {
            evaluate_ci(Lambda_perm, Lambda_true, method_name = method_name)
         } else {
            cat("Dimension mismatch for Lambda in CI eval\n")
         }
      } else {
         cat("Dimension mismatch for Lambda in CI eval\n")
      }
   }
   invisible(list(mse = mse, Lambda = Lambda, sparse_frac = sparse_frac))
}



# ---- Evaluate all ----
mgps_out      <- evaluate_stan_fit(mgps_file, "MGPS", Lambda_true = Lambda_true)
ssl_out       <- evaluate_stan_fit(ssl_file, "SSL", Lambda_true = Lambda_true)
horseshoe_out <- evaluate_stan_fit(horseshoe_file, "Horseshoe", Lambda_true = Lambda_true)
tebfar_out    <- evaluate_stan_fit(tebfar_file, "TEB-FAR", Lambda_true = Lambda_true)

# ---- Print summary of test-set MSEs (only for Stan models) ----
cat("Test-set MSEs:\n")
cat("MGPS:      ", mgps_out$mse, "\n")
cat("SSL:       ", ssl_out$mse, "\n")
cat("Horseshoe: ", horseshoe_out$mse, "\n")
cat("TEB-FAR:   ", tebfar_out$mse, "\n")

cat("\nFraction of near-zero loadings (< 0.05):\n")
cat("MGPS:      ", mgps_out$sparse_frac, "\n")
cat("SSL:       ", ssl_out$sparse_frac, "\n")
cat("Horseshoe: ", horseshoe_out$sparse_frac, "\n")
cat("TEB-FAR:   ", tebfar_out$sparse_frac, "\n")
