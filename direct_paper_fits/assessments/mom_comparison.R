# ----- Load sim2 data -----
sim2 <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000_1.rds")
X <- scale(sim2$Y, center = TRUE, scale = TRUE)

# ---- K=1 Method of Moments (MoM) ----
S <- cov(X)
e <- eigen(S)
d <- e$values
V <- e$vectors

sigma2_hat  <- mean(d[-1])
lambda_mag  <- sqrt(max(d[1] - sigma2_hat, 0))
lambda1_hat <- lambda_mag * V[, 1]
Sigma_hat   <- sigma2_hat * diag(ncol(X))
Sigma_inv   <- solve(Sigma_hat)
denom       <- as.numeric(crossprod(lambda1_hat, Sigma_inv %*% lambda1_hat) + 1)
w_vec       <- (Sigma_inv %*% lambda1_hat) / denom

eta_hat_1 <- as.numeric(X %*% w_vec)  # n-vector
F1_fit    <- outer(eta_hat_1, lambda1_hat)  # n x p
resid1_MoM <- X - F1_fit

# ---- Off-diagonal correlation histograms ----
get_offdiag <- function(M) M[lower.tri(M)]
cor_before <- cor(X)
cor_after_MoM <- cor(resid1_MoM)
r_before <- get_offdiag(cor_before)
r_after_MoM <- get_offdiag(cor_after_MoM)

cat("Mean absolute off-diagonal correlation (MoM):\n")
cat("   Before K=1:", round(mean(abs(r_before)), 4), "\n")
cat("   After  K=1:", round(mean(abs(r_after_MoM)), 4), "\n")

# ---- MoM DSC with permutation null ----
set.seed(42)
res_dsc_before <- dsc_with_permutation_null_obs(X, B = 500)
res_dsc_after  <- dsc_with_permutation_null_obs(resid1_MoM, B = 500)

cat("\n[MoM] DSC (original data):\n")
cat("Observed DSC:", round(res_dsc_before$dsc_obs, 3), "\n")
cat("Permutation null mean:", round(mean(res_dsc_before$dsc_null), 3), "\n")
cat("Null 95% interval:", round(quantile(res_dsc_before$dsc_null, c(0.025, 0.975)), 3), "\n")

cat("\n[MoM] DSC (after removing K=1):\n")
cat("Observed DSC:", round(res_dsc_after$dsc_obs, 3), "\n")
cat("Permutation null mean:", round(mean(res_dsc_after$dsc_null), 3), "\n")
cat("Null 95% interval:", round(quantile(res_dsc_after$dsc_null, c(0.025, 0.975)), 3), "\n")

# ---- Plot off-diagonal histograms before/after MoM ----
par(mfrow = c(1,2))
hist(r_before, breaks = 40, col = rgb(0.2,0.4,1,0.5), main = "MoM: Before", xlab = "Correlation", border = "white")
hist(r_after_MoM, breaks = 40, col = rgb(1,0.4,0.4,0.5), add = TRUE, border = "white")
legend("topright", legend = c("Before", "After MoM K=1"), fill = c(rgb(0.2,0.4,1,0.5), rgb(1,0.4,0.4,0.5)), border = NA)

# ---- Plot DSC nulls for before/after MoM ----
par(mfrow = c(1,2))
hist(res_dsc_before$dsc_null, breaks = 40, main = "MoM DSC Null (Before)", xlab = "DSC", col = "lightblue", border = "white")
abline(v = res_dsc_before$dsc_obs, col = "darkblue", lwd = 2)
legend("topright", legend = "Observed DSC", col = "darkblue", lwd = 2, bty = "n")

hist(res_dsc_after$dsc_null, breaks = 40, main = "MoM DSC Null (After K=1)", xlab = "DSC", col = "pink", border = "white")
abline(v = res_dsc_after$dsc_obs, col = "red", lwd = 2)
legend("topright", legend = "Observed DSC", col = "red", lwd = 2, bty = "n")

# ---- Compare to Bayesian fits: Overlay histograms ----
model_names <- c("MGSP", "HS", "SSL", "MASS", "Robust")
fit_list <- list(fit_MGSP, fit_HS, fit_SSL, fit_MASS, fit_RobustMGSP)

for (i in seq_along(model_names)) {
   Lambda_hat <- fit_list[[i]]$Lambda_hat
   eta_hat <- X %*% Lambda_hat %*% solve(t(Lambda_hat) %*% Lambda_hat)
   resid_bayes <- X - eta_hat %*% t(Lambda_hat)
   r_after_bayes <- get_offdiag(cor(resid_bayes))
   cat("\nMean abs off-diagonal correlation (", model_names[i], "):\n", sep = "")
   cat("   After K=1:", round(mean(abs(r_after_bayes)), 4), "\n")
   hist(r_before, breaks = 40, col = rgb(0.2,0.4,1,0.4), main = paste(model_names[i], "Before/After"), xlab = "Correlation", border = "white")
   hist(r_after_bayes, breaks = 40, col = rgb(1,0.7,0.2,0.5), add = TRUE, border = "white")
   legend("topright", legend = c("Before", "After K=1"), fill = c(rgb(0.2,0.4,1,0.4), rgb(1,0.7,0.2,0.5)), border = NA)
}
par(mfrow = c(1,1))
