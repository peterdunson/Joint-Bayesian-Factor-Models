# ---------------------------------------------------------------
# Visualize Stan Horseshoe Factor Model Loadings Matrix
# ---------------------------------------------------------------

library(rstan)
library(pheatmap)

# ---- Load Stan fit ----
scenario <- 1  # or 2, 3, etc. as appropriate
fit <- readRDS(sprintf("stan_horseshoe_fit_scen%d_5.rds", scenario))

# ---- Extract posterior draws ----
post <- rstan::extract(fit)
# Check parameter names and dimensions if unsure:
# str(post)

# ---- Posterior mean of factor loadings ----
# This assumes: post$Lambda is [iterations, P, K]
Lambda_post <- apply(post$Lambda, c(2, 3), mean)

# ---- Print and basic inspect ----
print(dim(Lambda_post))         # Should be [P, K]
print(round(Lambda_post, 2)[1:5,])  # First 5 rows

# ---- Visualize as heatmap (recommended) ----
pheatmap(Lambda_post,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Posterior Mean Loadings (Horseshoe Factor Model)",
         color = colorRampPalette(c("blue", "white", "red"))(50))

# ---- Optional: Image plot (base R) ----
image(t(Lambda_post), main = "Posterior Mean Factor Loadings (Horseshoe)",
      xlab = "Variables", ylab = "Factors", axes = FALSE, col = heat.colors(30))
box()

# ---- Optional: Print histogram of loadings ----
hist(as.numeric(Lambda_post), breaks = 50, main = "Histogram of Posterior Mean Loadings")

# ---- Optional: Show fraction of near-zero loadings (sparsity) ----
sparse_frac <- mean(abs(Lambda_post) < 0.05)
cat("Fraction of near-zero loadings (< 0.05):", round(sparse_frac, 3), "\n")




# ---- Calculate Mean Squared Error (MSE) ----

# ---------------------------------------------------------------
# Evaluate Horseshoe Factor Model (Stan) on Simulated Data
# ---------------------------------------------------------------

library(rstan)
library(tidyverse)

# --- USER SETTINGS ---
scenario <- 1
n_train  <- 1000   # use 200 or 1000 as appropriate

sim_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations"
horseshoe_fit_file <- sprintf("%s/stan_horseshoe_fit_scen%d_5.rds", sim_dir, scenario)

# --- Load simulation data ---
sim_file <- sprintf("%s/sim_scen%d_%d.rds", sim_dir, scenario, n_train)
sim <- readRDS(sim_file)
yX <- sim$Y      # n x (p+1): first column is y
X <- yX[, -1, drop=FALSE]
y <- yX[, 1]

# --- Standardize predictors and outcome (centering and scaling) ---
X <- scale(X)
y <- scale(y, center = TRUE, scale = TRUE)

# --- Train/test split (75/25), set seed for reproducibility ---
set.seed(42)
n <- nrow(X)
test_frac <- 0.25
test_idx <- sample(1:n, size = floor(n * test_frac))
train_idx <- setdiff(1:n, test_idx)

X_train <- X[train_idx, , drop=FALSE]
y_train <- y[train_idx]
X_test  <- X[test_idx, , drop=FALSE]
y_test  <- y[test_idx]

# --- Load Stan fit and extract posterior means ---
if (file.exists(horseshoe_fit_file)) {
   horseshoe_fit <- readRDS(horseshoe_fit_file)
   horseshoe_post <- rstan::extract(horseshoe_fit)
   Lambda_horseshoe <- apply(horseshoe_post$Lambda, c(2,3), mean)  # (p+1) x K
   psi_horseshoe <- apply(horseshoe_post$psi, 2, mean)             # (p+1)
   
   # --- Prediction function (same as TEB-FAR etc.) ---
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
   
   # --- Predict and compute MSE on test set ---
   horseshoe_pred <- predict_y_from_factors(X_test, Lambda_horseshoe, psi_horseshoe)
   horseshoe_mse <- mean((y_test - horseshoe_pred)^2)
   
   cat("Horseshoe Factor Model Test-set MSE (Scenario =", scenario, ", n_train =", n_train, "):\n")
   cat("MSE:", horseshoe_mse, "\n")
} else {
   cat("Stan fit not found:", horseshoe_fit_file, "\n")
}




# ---- Calculate Coverage ----


# coverage_eval_horseshoe.R
# ----------------------------------------------------------
# Evaluates 95% credible interval coverage for factor loadings
# for the Horseshoe Factor Model (Stan)
# ----------------------------------------------------------

library(rstan)

# --- USER SETTINGS ---
scenario <- 1
n_train <- 1000

sim_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations"
horseshoe_fit_file <- sprintf("%s/stan_horseshoe_fit_scen%d_5.rds", sim_dir, scenario)

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

# ----------------------------------------------
# --- Evaluate Horseshoe Factor Model (Stan) ---
# ----------------------------------------------
if (file.exists(horseshoe_fit_file)) {
   horseshoe_fit <- readRDS(horseshoe_fit_file)
   horseshoe_post <- rstan::extract(horseshoe_fit)
   # Adjust for dimension order if needed
   if (length(dim(horseshoe_post$Lambda)) == 3 && dim(horseshoe_post$Lambda)[2] == dim(Lambda_true)[1]) {
      Lambda_horseshoe_post <- horseshoe_post$Lambda
   } else if (length(dim(horseshoe_post$Lambda)) == 3 && dim(horseshoe_post$Lambda)[3] == dim(Lambda_true)[1]) {
      Lambda_horseshoe_post <- aperm(horseshoe_post$Lambda, c(1, 3, 2))
   } else {
      stop("Unexpected Lambda dimensions for Horseshoe!")
   }
   horseshoe_coverage <- coverage_rate(Lambda_horseshoe_post, Lambda_true)
   cat(sprintf("Horseshoe Factor Model Lambda 95%% coverage: %.2f%%\n", 100 * horseshoe_coverage))
} else {
   cat("Horseshoe fit file not found.\n")
}


#ci symmetry check



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

# -------------- HORSESHOE FACTOR MODEL --------------------
if (file.exists(horseshoe_fit_file)) {
   horseshoe_fit <- readRDS(horseshoe_fit_file)
   horseshoe_post <- rstan::extract(horseshoe_fit)
   # Ensure posterior Lambda dims match truth
   if (all(dim(horseshoe_post$Lambda)[2:3] == dim(Lambda_true))) {
      res_horseshoe <- evaluate_ci(horseshoe_post$Lambda, Lambda_true, method_name = "Horseshoe")
   } else if (length(dim(horseshoe_post$Lambda)) == 3 &&
              dim(horseshoe_post$Lambda)[3] == dim(Lambda_true)[1]) {
      # Sometimes Stan output is [iterations, K, p], need to aperm
      Lambda_perm <- aperm(horseshoe_post$Lambda, c(1, 3, 2))
      if (all(dim(Lambda_perm)[2:3] == dim(Lambda_true))) {
         res_horseshoe <- evaluate_ci(Lambda_perm, Lambda_true, method_name = "Horseshoe")
      } else {
         cat("Dimension mismatch for Horseshoe: check K and simulation!\n")
         cat("Posterior Lambda dim: ", paste(dim(horseshoe_post$Lambda), collapse = " x "), "\n")
         cat("True Lambda dim: ", paste(dim(Lambda_true), collapse = " x "), "\n")
      }
   } else {
      cat("Dimension mismatch for Horseshoe: check K and simulation!\n")
      cat("Posterior Lambda dim: ", paste(dim(horseshoe_post$Lambda), collapse = " x "), "\n")
      cat("True Lambda dim: ", paste(dim(Lambda_true), collapse = " x "), "\n")
   }
}



