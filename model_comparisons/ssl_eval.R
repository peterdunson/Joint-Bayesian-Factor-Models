# ---------------------------------------------------------------
# Evaluate Spike-and-Slab (SSL) Factor Model (Stan) on Simulated Data
# ---------------------------------------------------------------

library(rstan)
library(pheatmap)
library(tidyverse)

# --- USER SETTINGS ---
scenario <- 1
n_train  <- 1000   # use 200 or 1000 as appropriate
K        <- 5      # number of factors (should match model fit)
sim_dir  <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations"

ssl_fit_file <- sprintf("%s/ssl_factor_fit_scen%d_%d.rds", sim_dir, scenario, K)
sim_file     <- sprintf("%s/sim_scen%d_%d.rds", sim_dir, scenario, n_train)

# --- Load simulation data and Stan fit ---
if (!file.exists(sim_file)) stop("Simulation file not found: ", sim_file)
if (!file.exists(ssl_fit_file)) stop("Stan fit not found: ", ssl_fit_file)

sim <- readRDS(sim_file)
yX <- sim$Y      # n x (p+1): first column is y
Lambda_true <- sim$Lambda # (p+1) x K, or check dimensions

# --- Standardize predictors and outcome ---
X <- yX[, -1, drop=FALSE]
y <- yX[, 1]
X <- scale(X)
y <- scale(y, center = TRUE, scale = TRUE)

# --- Train/test split ---
set.seed(42)
n <- nrow(X)
test_frac <- 0.25
test_idx <- sample(1:n, size = floor(n * test_frac))
train_idx <- setdiff(1:n, test_idx)
X_train <- X[train_idx, , drop=FALSE]
y_train <- y[train_idx]
X_test  <- X[test_idx, , drop=FALSE]
y_test  <- y[test_idx]

# --- Extract Stan posterior ---
ssl_fit  <- readRDS(ssl_fit_file)
ssl_post <- rstan::extract(ssl_fit)
Lambda_ssl_post <- ssl_post$Lambda # [iterations, p+1, K]
psi_ssl_post    <- ssl_post$psi    # [iterations, p+1]

# --- Posterior mean loadings & residuals ---
Lambda_ssl <- apply(Lambda_ssl_post, c(2, 3), mean)  # (p+1) x K
psi_ssl    <- apply(psi_ssl_post, 2, mean)           # (p+1)

# ---- Visualization ----
pheatmap(Lambda_ssl,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Posterior Mean Loadings (SSL Factor Model)",
         color = colorRampPalette(c("blue", "white", "red"))(50))
hist(as.numeric(Lambda_ssl), breaks = 50, main = "Histogram of Posterior Mean Loadings")
sparse_frac <- mean(abs(Lambda_ssl) < 0.05)
cat("Fraction of near-zero loadings (< 0.05):", round(sparse_frac, 3), "\n")

# ---- Prediction function (same as horseshoe/TEB-FAR) ----
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
ssl_pred <- predict_y_from_factors(X_test, Lambda_ssl, psi_ssl)
ssl_mse <- mean((y_test - ssl_pred)^2)
cat("SSL Factor Model Test-set MSE (Scenario =", scenario, ", n_train =", n_train, "):\n")
cat("MSE:", ssl_mse, "\n")

# ---- Coverage and CI symmetry ----
evaluate_ci <- function(posterior_array, Lambda_true, method_name = "SSL") {
   lower <- apply(posterior_array, c(2, 3), quantile, probs = 0.025)
   upper <- apply(posterior_array, c(2, 3), quantile, probs = 0.975)
   mean_est <- apply(posterior_array, c(2, 3), mean)
   coverage <- (Lambda_true >= lower) & (Lambda_true <= upper)
   coverage_rate <- mean(coverage)
   upper_diff <- upper - mean_est
   lower_diff <- mean_est - lower
   symmetry <- abs(upper_diff - lower_diff) / (upper_diff + lower_diff + 1e-10)
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

# Ensure posterior dims match true dims
if (all(dim(Lambda_ssl_post)[2:3] == dim(Lambda_true))) {
   res_ssl <- evaluate_ci(Lambda_ssl_post, Lambda_true, method_name = "SSL")
} else if (length(dim(Lambda_ssl_post)) == 3 &&
           dim(Lambda_ssl_post)[3] == dim(Lambda_true)[1]) {
   Lambda_perm <- aperm(Lambda_ssl_post, c(1, 3, 2))
   if (all(dim(Lambda_perm)[2:3] == dim(Lambda_true))) {
      res_ssl <- evaluate_ci(Lambda_perm, Lambda_true, method_name = "SSL")
   } else {
      cat("Dimension mismatch for SSL: check K and simulation!\n")
   }
} else {
   cat("Dimension mismatch for SSL: check K and simulation!\n")
}



