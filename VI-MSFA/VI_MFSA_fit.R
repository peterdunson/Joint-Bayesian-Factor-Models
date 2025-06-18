# -------------------------------------------------------------
# Evaluate VI-MSFA (CAVI/SVI) on Simulated Data: Test-set MSE,
# Factor Loadings Heatmap, Histogram, and Sparsity
# -------------------------------------------------------------

# --- LIBRARIES ---
library(VIMSFA)             # For VI-MSFA (CAVI/SVI)
library(pheatmap)           # For heatmaps

# ---- LOAD SIM DATA ----
sim_path <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen1_1000.rds"
sim <- readRDS(sim_path)
X_full <- sim$Y[, -1, drop=FALSE]   # predictors (n x p)
y_full <- sim$Y[, 1]                # outcome (length n)
n <- nrow(X_full)
p <- ncol(X_full)
J <- 5   # num factors (should match K in your sim)

# ---- STANDARDIZE DATA (mean 0, variance 1) ----
X <- scale(X_full, center=TRUE, scale=TRUE)
y <- as.vector(scale(y_full, center=TRUE, scale=TRUE))

# ---- TRAIN/TEST SPLIT (75/25) ----
set.seed(42)
test_frac <- 0.25
test_idx <- sample(1:n, size = floor(n * test_frac))
train_idx <- setdiff(1:n, test_idx)
X_train <- X[train_idx, , drop=FALSE]
y_train <- y[train_idx]
X_test  <- X[test_idx, , drop=FALSE]
y_test  <- y[test_idx]

# ---- FACTOR SCORE PROJECTION FUNCTION ----
compute_factor_scores <- function(X_new, lambda, psi) {
   Psi_inv <- diag(1 / psi)
   A <- t(lambda) %*% Psi_inv %*% lambda + diag(ncol(lambda))
   B <- t(lambda) %*% Psi_inv
   t(apply(X_new, 1, function(x) as.numeric(solve(A, B %*% x))))
}

# ---- FIT VI-MSFA (CAVI) ON TRAIN SET ----
cavi_est <- cavi_fa(X_train, J, scale=FALSE)
F_train_cavi <- compute_factor_scores(X_train, cavi_est$mean_lambda, cavi_est$mean_psi)
F_test_cavi  <- compute_factor_scores(X_test,  cavi_est$mean_lambda, cavi_est$mean_psi)

# ---- FIT VI-MSFA (SVI) ON TRAIN SET ----
svi_est <- svi_fa(X_train, J, scale=FALSE)
F_train_svi <- compute_factor_scores(X_train, svi_est$mean_lambda, svi_est$mean_psi)
F_test_svi  <- compute_factor_scores(X_test,  svi_est$mean_lambda, svi_est$mean_psi)

# ---- REGRESSION FOR TEST-SET MSE ----
get_mse_from_factors <- function(F_train, F_test, y_train, y_test, method="") {
   F_train_df <- as.data.frame(F_train)
   colnames(F_train_df) <- paste0("F", 1:ncol(F_train))
   F_test_df <- as.data.frame(F_test)
   colnames(F_test_df) <- paste0("F", 1:ncol(F_test))
   reg <- lm(y_train ~ ., data = F_train_df)
   y_pred <- predict(reg, newdata = F_test_df)
   mse <- mean((y_test - y_pred)^2)
   cat(sprintf("[%s] Test-set MSE: %.4f\n", method, mse))
   mse
}

mse_cavi <- get_mse_from_factors(F_train_cavi, F_test_cavi, y_train, y_test, method="VI-MSFA CAVI")
mse_svi  <- get_mse_from_factors(F_train_svi,  F_test_svi,  y_train, y_test, method="VI-MSFA SVI")

# ---- VISUALIZE FACTOR LOADINGS (HEATMAP, HISTOGRAM) ----
# CAVI
Lambda_cavi <- cavi_est$mean_lambda
pheatmap(Lambda_cavi, cluster_rows=FALSE, cluster_cols=FALSE,
         main="VI-MSFA (CAVI) Mean Factor Loadings")
hist(as.numeric(Lambda_cavi), breaks=50, main="VI-MSFA (CAVI) Mean Loadings", xlab="Loading Value")
sparse_frac_cavi <- mean(abs(Lambda_cavi) < 0.05)
cat(sprintf("[VI-MSFA CAVI] Fraction near-zero loadings (<0.05): %.2f\n", sparse_frac_cavi))

# SVI
Lambda_svi <- svi_est$mean_lambda
pheatmap(Lambda_svi, cluster_rows=FALSE, cluster_cols=FALSE,
         main="VI-MSFA (SVI) Mean Factor Loadings")
hist(as.numeric(Lambda_svi), breaks=50, main="VI-MSFA (SVI) Mean Loadings", xlab="Loading Value")
sparse_frac_svi <- mean(abs(Lambda_svi) < 0.05)
cat(sprintf("[VI-MSFA SVI] Fraction near-zero loadings (<0.05): %.2f\n", sparse_frac_svi))
