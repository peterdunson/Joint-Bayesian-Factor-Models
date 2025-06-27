# --- Load simulation data ---
sim2 <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Y_all <- sim2$Y           # n x (p+1), including y as first column
X     <- sweep(Y_all, 2, colMeans(Y_all), FUN = "-")  # center all columns

# --- 1. K=1 MoM estimator (top factor) ---
S   <- cov(X)
e   <- eigen(S)
d   <- e$values
V   <- e$vectors

sigma2_hat  <- mean(d[-1])
lambda_mag  <- sqrt(max(d[1] - sigma2_hat, 0))
lambda1_hat <- lambda_mag * V[, 1]
Sigma_hat   <- sigma2_hat * diag(ncol(X))
Sigma_inv   <- solve(Sigma_hat)
denom       <- as.numeric(crossprod(lambda1_hat, Sigma_inv %*% lambda1_hat) + 1)
w_vec       <- (Sigma_inv %*% lambda1_hat) / denom

eta_hat_1   <- as.numeric(X %*% w_vec)   # n-vector

# --- 2. Remove fitted factor: project out rank-1 structure ---
F1_fit   <- outer(eta_hat_1, lambda1_hat)  # n x p+1
resid1   <- X - F1_fit

# --- 3. Compare correlation distributions: before vs after ---
cor_before <- cor(X)
cor_after  <- cor(resid1)

offdiag <- function(M) M[lower.tri(M)]
hist(offdiag(cor_before), breaks = 40, col = "#AAAAFF", main = "Off-diagonal correlations (Before/After K=1 removal)", xlab = "Correlation", freq = FALSE)
hist(offdiag(cor_after), breaks = 40, col = "#FFAAAA", add = TRUE, freq = FALSE)
legend("topright", legend = c("Before", "After K=1 MoM"), fill = c("#AAAAFF", "#FFAAAA"))

# --- 4. Summarize reduction in correlation magnitude ---
cat(sprintf("Mean absolute off-diag correlation before: %.3f\n", mean(abs(offdiag(cor_before)))))
cat(sprintf("Mean absolute off-diag correlation after : %.3f\n", mean(abs(offdiag(cor_after)))))





