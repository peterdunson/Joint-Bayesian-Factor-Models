# -------------------------------------------------------------
# Evaluate Frequentist Methods (+ VI-MSFA) on Simulated Data
# -------------------------------------------------------------

library(glmnet)
library(VIMSFA)
library(pheatmap)

# --- LOAD SIM DATA ---
sim_path <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds"
sim <- readRDS(sim_path)
X_full <- sim$Y[, -1, drop=FALSE]   # predictors (n x p)
y_full <- sim$Y[, 1]                # outcome (length n)
n <- nrow(X_full)
p <- ncol(X_full)
J <- 5   # num factors

# --- STANDARDIZE DATA (mean 0, variance 1) ---
X <- scale(X_full, center=TRUE, scale=TRUE)
y <- as.vector(scale(y_full, center=TRUE, scale=TRUE))

# --- TRAIN/TEST SPLIT (75/25) ---
set.seed(42)
test_frac <- 0.25
test_idx <- sample(1:n, size = floor(n * test_frac))
train_idx <- setdiff(1:n, test_idx)
X_train <- X[train_idx, , drop=FALSE]
y_train <- y[train_idx]
X_test  <- X[test_idx, , drop=FALSE]
y_test  <- y[test_idx]

# -------------- LASSO --------------
lasso_mod <- cv.glmnet(X_train, y_train, alpha = 1)
lasso_pred <- predict(lasso_mod, newx = X_test, s = "lambda.min")
lasso_mse <- mean((y_test - lasso_pred)^2)
cat(sprintf("[Lasso] Test-set MSE: %.4f\n", lasso_mse))

# -------------- RIDGE --------------
ridge_mod <- cv.glmnet(X_train, y_train, alpha = 0)
ridge_pred <- predict(ridge_mod, newx = X_test, s = "lambda.min")
ridge_mse <- mean((y_test - ridge_pred)^2)
cat(sprintf("[Ridge] Test-set MSE: %.4f\n", ridge_mse))

# -------------- OLS --------------
X_train_df <- as.data.frame(X_train)
X_test_df  <- as.data.frame(X_test)
colnames(X_test_df) <- colnames(X_train_df)
ols_mod <- lm(y_train ~ ., data = X_train_df)
ols_pred <- predict(ols_mod, newdata = X_test_df)
ols_mse <- mean((y_test - ols_pred)^2)
cat(sprintf("[OLS] Test-set MSE: %.4f\n", ols_mse))

# -------------- PCA + REGRESSION --------------
pca <- prcomp(X_train, center = FALSE, scale. = FALSE)
expl_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
M <- which(expl_var >= 0.95)[1]
X_train_pca <- pca$x[, 1:M, drop = FALSE]
X_test_pca  <- predict(pca, newdata = X_test)[, 1:M, drop = FALSE]
X_train_pca_df <- as.data.frame(X_train_pca)
X_test_pca_df  <- as.data.frame(X_test_pca)
colnames(X_test_pca_df) <- colnames(X_train_pca_df)
pca_reg_mod <- lm(y_train ~ ., data = X_train_pca_df)
pca_pred <- predict(pca_reg_mod, newdata = X_test_pca_df)
pca_mse <- mean((y_test - pca_pred)^2)
cat(sprintf("[PCA+Regress] Test-set MSE: %.4f\n", pca_mse))

# -------------- VI-MSFA (CAVI) --------------
cavi_est <- cavi_fa(X_train, J, scale=FALSE)
compute_factor_scores <- function(X_new, lambda, psi) {
   Psi_inv <- diag(1 / psi)
   A <- t(lambda) %*% Psi_inv %*% lambda + diag(ncol(lambda))
   B <- t(lambda) %*% Psi_inv
   t(apply(X_new, 1, function(x) as.numeric(solve(A, B %*% x))))
}
F_train_cavi <- compute_factor_scores(X_train, cavi_est$mean_lambda, cavi_est$mean_psi)
F_test_cavi  <- compute_factor_scores(X_test,  cavi_est$mean_lambda, cavi_est$mean_psi)
F_train_df <- as.data.frame(F_train_cavi)
colnames(F_train_df) <- paste0("F", 1:ncol(F_train_cavi))
F_test_df <- as.data.frame(F_test_cavi)
colnames(F_test_df) <- paste0("F", 1:ncol(F_test_cavi))
reg_cavi <- lm(y_train ~ ., data = F_train_df)
y_pred_cavi <- predict(reg_cavi, newdata = F_test_df)
mse_cavi <- mean((y_test - y_pred_cavi)^2)
cat(sprintf("[VI-MSFA CAVI] Test-set MSE: %.4f\n", mse_cavi))

# -------------- VI-MSFA (SVI) --------------
svi_est <- svi_fa(X_train, J, scale=FALSE)
F_train_svi <- compute_factor_scores(X_train, svi_est$mean_lambda, svi_est$mean_psi)
F_test_svi  <- compute_factor_scores(X_test,  svi_est$mean_lambda, svi_est$mean_psi)
F_train_df_svi <- as.data.frame(F_train_svi)
colnames(F_train_df_svi) <- paste0("F", 1:ncol(F_train_svi))
F_test_df_svi <- as.data.frame(F_test_svi)
colnames(F_test_df_svi) <- paste0("F", 1:ncol(F_test_svi))
reg_svi <- lm(y_train ~ ., data = F_train_df_svi)
y_pred_svi <- predict(reg_svi, newdata = F_test_df_svi)
mse_svi <- mean((y_test - y_pred_svi)^2)
cat(sprintf("[VI-MSFA SVI] Test-set MSE: %.4f\n", mse_svi))

# -------------- SUMMARY --------------
cat("\nSummary of Test-set MSEs:\n")
cat(sprintf("Lasso:       %.4f\n", lasso_mse))
cat(sprintf("Ridge:       %.4f\n", ridge_mse))
cat(sprintf("OLS:         %.4f\n", ols_mse))
cat(sprintf("PCA+Regress: %.4f\n", pca_mse))
cat(sprintf("VI-MSFA CAVI:%.4f\n", mse_cavi))
cat(sprintf("VI-MSFA SVI: %.4f\n", mse_svi))


# -------------- VI-MSFA CAVI: Sparsity and Visualizations --------------
Lambda_cavi <- cavi_est$mean_lambda
frac_nonzero_cavi <- mean(abs(Lambda_cavi) >= 0.05)
cat(sprintf("[VI-MSFA CAVI] Fraction of non-zero loadings (>=0.05): %.2f\n", frac_nonzero_cavi))
pheatmap(Lambda_cavi, cluster_rows=FALSE, cluster_cols=FALSE,
         main="VI-MSFA (CAVI) Mean Factor Loadings")
hist(as.numeric(Lambda_cavi), breaks=50, main="VI-MSFA (CAVI) Mean Loadings", xlab="Loading Value")

# -------------- VI-MSFA SVI: Sparsity and Visualizations --------------
Lambda_svi <- svi_est$mean_lambda
frac_nonzero_svi <- mean(abs(Lambda_svi) >= 0.05)
cat(sprintf("[VI-MSFA SVI] Fraction of non-zero loadings (>=0.05): %.2f\n", frac_nonzero_svi))
pheatmap(Lambda_svi, cluster_rows=FALSE, cluster_cols=FALSE,
         main="VI-MSFA (SVI) Mean Factor Loadings")
hist(as.numeric(Lambda_svi), breaks=50, main="VI-MSFA (SVI) Mean Loadings", xlab="Loading Value")


