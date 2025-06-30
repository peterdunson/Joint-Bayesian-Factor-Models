# --- Load and standardize data ---
X <- scale(readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")$Y)

# --- K=1 MoM estimator (use as BLUP) ---
S   <- cov(X)
e   <- eigen(S)
d   <- e$values
V   <- e$vectors

sigma2_hat  <- mean(d[-1])
lambda_hat  <- sqrt(max(d[1] - sigma2_hat, 0)) * V[, 1]
w_vec       <- lambda_hat / (sum(lambda_hat^2) + sigma2_hat)

# --- BLUP for all subjects ---
eta_blup <- as.numeric(X %*% w_vec)

# --- Project out BLUP factor (decorrelate) ---
F_fit   <- outer(eta_blup, lambda_hat)   # n x p
X_resid <- X - F_fit

# --- Off-diagonal correlations ---
get_offdiag <- function(M) M[lower.tri(M)]
r_before <- get_offdiag(cor(X))
r_after  <- get_offdiag(cor(X_resid))

# --- Plot overlayed histograms ---
breaks <- seq(-1, 1, length.out = 41)
hist(r_before, breaks = breaks, col = rgb(0.2, 0.4, 1, 0.5),
     main = "Off-diagonal correlations\nBefore and After BLUP Factor Removal",
     xlab = "Correlation", ylab = "Frequency", border = "white")
hist(r_after, breaks = breaks, col = rgb(1, 0.4, 0.4, 0.5), add = TRUE, border = "white")
legend("topright", legend = c("Before", "After BLUP removal"),
       fill = c(rgb(0.2, 0.4, 1, 0.5), rgb(1, 0.4, 0.4, 0.5)), border = NA)

cat("Mean |cor| before:", mean(abs(r_before)), "\n")
cat("Mean |cor| after :", mean(abs(r_after)), "\n")
