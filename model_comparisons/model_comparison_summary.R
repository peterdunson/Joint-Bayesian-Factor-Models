# model_comparison_summary.R
# Compare TEB-FAR, Baseline Factor Model, Lasso, Ridge, OLS, PCA+Regression

library(rstan)
library(glmnet)
library(tidyverse)

# --- USER SETTINGS ---
scenario <- 1        # 1, 2, or 3
n_train  <- 1000     # 200 or 1000

sim_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations"
tebfar_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/teb_far"
baseline_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model"

sim_file <- sprintf("%s/sim_scen%d_%d.rds", sim_dir, scenario, n_train)
sim <- readRDS(sim_file)
yX <- sim$Y      # n x (p+1): first column is y
X <- yX[, -1, drop=FALSE]
y <- yX[, 1]

# --- Standardize predictors and outcome (centering and scaling) ---
X <- scale(X)
y <- scale(y, center = TRUE, scale = TRUE)

# --- Train/test split (75/25) ---
set.seed(42)
n <- nrow(X)
test_frac <- 0.25
test_idx <- sample(1:n, size = floor(n * test_frac))
train_idx <- setdiff(1:n, test_idx)

X_train <- X[train_idx, , drop=FALSE]
y_train <- y[train_idx]
X_test  <- X[test_idx, , drop=FALSE]
y_test  <- y[test_idx]

# ---------------------------------------------------------------
# ---- TEB-FAR (Stan) ----
# ---------------------------------------------------------------
setwd(tebfar_dir)
tebfar_fit_file <- sprintf("%s/stan_tebfar_fit_scen%d.rds", sim_dir, scenario)
if (file.exists(tebfar_fit_file)) {
   tebfar_fit <- readRDS(tebfar_fit_file)
   tebfar_post <- rstan::extract(tebfar_fit)
   Lambda_tebfar <- apply(tebfar_post$Lambda, c(2,3), mean)  # (p+1) x K
   psi_tebfar <- apply(tebfar_post$psi, 2, mean)             # (p+1)
   # Prediction
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
   tebfar_pred <- predict_y_from_factors(X_test, Lambda_tebfar, psi_tebfar)
   tebfar_mse <- mean((y_test - tebfar_pred)^2)
} else {
   tebfar_mse <- NA
}

# ---------------------------------------------------------------
# ---- Baseline Factor Model (Stan) ----
# ---------------------------------------------------------------
setwd(baseline_dir)
baseline_fit_file <- sprintf("%s/mgps_fit_scen%d.rds", sim_dir, scenario)
if (file.exists(baseline_fit_file)) {
   baseline_fit <- readRDS(baseline_fit_file)
   baseline_post <- rstan::extract(baseline_fit)
   Lambda_baseline <- apply(baseline_post$Lambda, c(2,3), mean)
   psi_baseline <- apply(baseline_post$psi, 2, mean)
   baseline_pred <- predict_y_from_factors(X_test, Lambda_baseline, psi_baseline)
   baseline_mse <- mean((y_test - baseline_pred)^2)
} else {
   baseline_mse <- NA
}

# ---------------------------------------------------------------
# ---- Lasso ----
# ---------------------------------------------------------------
lasso_mod <- cv.glmnet(X_train, y_train, alpha = 1)
lasso_pred <- predict(lasso_mod, newx = X_test, s = "lambda.min")
lasso_mse <- mean((y_test - lasso_pred)^2)

# ---------------------------------------------------------------
# ---- Ridge ----
# ---------------------------------------------------------------
ridge_mod <- cv.glmnet(X_train, y_train, alpha = 0)
ridge_pred <- predict(ridge_mod, newx = X_test, s = "lambda.min")
ridge_mse <- mean((y_test - ridge_pred)^2)

# ---------------------------------------------------------------
# ---- OLS ----
# ---------------------------------------------------------------
# Use a data frame with the correct column names for training and testing
X_train_df <- as.data.frame(X_train)
X_test_df <- as.data.frame(X_test)
colnames(X_test_df) <- colnames(X_train_df)  # Ensure matching names

ols_mod <- lm(y_train ~ ., data = X_train_df)
ols_pred <- predict(ols_mod, newdata = X_test_df)
ols_mse <- mean((y_test - ols_pred)^2)

# ---------------------------------------------------------------
# ---- PCA + Regression ----
# ---------------------------------------------------------------
# PCA on training set only, then regress on top M PCs
pca <- prcomp(X_train, center = FALSE, scale. = FALSE)
expl_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
M <- which(expl_var >= 0.95)[1]
X_train_pca <- pca$x[, 1:M, drop = FALSE]
X_test_pca <- predict(pca, newdata = X_test)[, 1:M, drop = FALSE]

# Convert to data frames with consistent column names
X_train_pca_df <- as.data.frame(X_train_pca)
X_test_pca_df <- as.data.frame(X_test_pca)
colnames(X_test_pca_df) <- colnames(X_train_pca_df)

pca_reg_mod <- lm(y_train ~ ., data = X_train_pca_df)
pca_pred <- predict(pca_reg_mod, newdata = X_test_pca_df)
pca_mse <- mean((y_test - pca_pred)^2)
# ---------------------------------------------------------------
# ---- Results ----
# ---------------------------------------------------------------
cat("Test-set MSEs (Scenario =", scenario, ", n_train =", n_train, ")\n")
cat("TEB-FAR:     ", tebfar_mse, "\n")
cat("Baseline:    ", baseline_mse, "\n")
cat("Lasso:       ", lasso_mse, "\n")
cat("Ridge:       ", ridge_mse, "\n")
cat("OLS:         ", ols_mse, "\n")
cat("PCA+Regress: ", pca_mse, "\n")
