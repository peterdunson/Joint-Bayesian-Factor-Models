# ---------------------------------------------------------------
# Evaluate All Factor Models on Simulated Data
# ---------------------------------------------------------------

library(rstan)
library(tidyverse)
library(abind)

set.seed(42)
scenario <- 2
n_train <- 1000
K <- 5
sim_dir <- "simulations"

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

evaluate_ci <- function(posterior_array, y_true, method_name = "TEB-FAR") {
   # posterior_array: [iterations, rows, cols]
   # Lambda_true: [rows, cols]
   lower <- apply(posterior_array, c(2), quantile, probs = 0.025)
   upper <- apply(posterior_array, c(2), quantile, probs = 0.975)
   mean_est <- apply(posterior_array, c(2), mean)
   
   coverage <- (y_true >= lower) & (y_true <= upper)
   coverage_rate <- mean(coverage)
   
   # CI symmetry: (upper - mean) vs (mean - lower)
   upper_diff <- upper - mean_est
   lower_diff <- mean_est - lower
   symmetry <- abs(upper_diff - lower_diff) / (upper_diff + lower_diff + 1e-10)  # Avoid 0/0
   
   cat(sprintf(
      "%s Y 95%% coverage: %.2f%%\n", method_name, 100 * coverage_rate
   ))
   
   cat(sprintf(
      "%s average CI symmetry (0 = perfect): %.3f\n", method_name, mean(symmetry)
   ))
   
   cat(sprintf(
      "%s mean squared prediction error in y: %.3f\n", method_name, mean((y_true - mean_est)^2)
   ))
}

fit <- readRDS(tebfar_file)
post <- rstan::extract(fit)
parray <- do.call(rbind, map(1:length(post$lp__), function(x){
   Lambda <- apply(post$Lambda, c(2,3), mean)
   psi <- apply(post$psi, 2, mean)
   
   y_vals = predict_y_from_factors(X_test, post$Lambda[x,,], 
                                   post$psi[x,])
   return(y_vals)
}))
evaluate_ci(parray, y_test, "TEB-FAR")

fit <- readRDS(mgps_file)
post <- rstan::extract(fit)
parray <- do.call(rbind, map(1:length(post$lp__), function(x){
   Lambda <- apply(post$Lambda, c(2,3), mean)
   psi <- apply(post$psi, 2, mean)
   
   y_vals = predict_y_from_factors(X_test, post$Lambda[x,,], 
                                   post$psi[x,])
   return(y_vals)
}))
evaluate_ci(parray, y_test, "MGSP")
