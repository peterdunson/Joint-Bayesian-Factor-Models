# --- Prepare matrix of phthalate variables (only adults, complete cases) ---
# Remove SEQN and RIDAGEYR, keep only the 19 variables
X_raw <- as.matrix(phthalates_adults[, phthalate_vars])

# Standardize each variable (center and scale)
X <- scale(X_raw)

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
F1_fit   <- outer(eta_hat_1, lambda1_hat)  # n x p
resid1   <- X - F1_fit

# --- 3. Compare correlation distributions: before vs after ---
cor_before <- cor(X)
cor_after  <- cor(resid1)






# --- Compute off-diagonal correlations ---
get_offdiag <- function(M) M[lower.tri(M)]
cor_before <- cor(X)
cor_after  <- cor(resid1)
r_before <- get_offdiag(cor_before)
r_after  <- get_offdiag(cor_after)

# --- Plot overlayed histograms ---
breaks <- seq(-1, 1, length.out = 41)  # 40 equal-width bins

hist(r_before,
     breaks = breaks,
     col = rgb(0.2, 0.4, 1, 0.5),   # blue, semi-transparent
     main = "Off-diagonal correlations before/after K=1 removal",
     xlab = "Correlation",
     ylab = "Frequency",
     border = "white")

hist(r_after,
     breaks = breaks,
     col = rgb(1, 0.4, 0.4, 0.5),   # red, semi-transparent
     add = TRUE,
     border = "white")

legend("topright",
       legend = c("Before", "After K=1 MoM"),
       fill = c(rgb(0.2, 0.4, 1, 0.5), rgb(1, 0.4, 0.4, 0.5)),
       border = NA)



# --- Print summary statistics and samples of off-diagonal correlations ---

cat("Summary of off-diagonal correlations BEFORE K=1 removal:\n")
print(summary(r_before))

cat("\nSummary of off-diagonal correlations AFTER K=1 removal:\n")
print(summary(r_after))

cat("\nFirst 10 off-diagonal correlations BEFORE:\n")
print(round(head(r_before, 10), 4))

cat("\nFirst 10 off-diagonal correlations AFTER:\n")
print(round(head(r_after, 10), 4))







set.seed(42)
res_dsc_before <- dsc_with_permutation_null_obs(X, B = 1000)

cat("\nDSC (original data):\n")
cat("Observed DSC:", round(res_dsc_before$dsc_obs, 3), "\n")
cat("Permutation null mean:", round(mean(res_dsc_before$dsc_null), 3), "\n")
cat("Null 95% interval:", round(quantile(res_dsc_before$dsc_null, c(0.025, 0.975)), 3), "\n")




set.seed(42)
res_dsc_after <- dsc_with_permutation_null_obs(resid1, B = 1000)

cat("\nDSC (after removing K=1):\n")
cat("Observed DSC:", round(res_dsc_after$dsc_obs, 3), "\n")
cat("Permutation null mean:", round(mean(res_dsc_after$dsc_null), 3), "\n")
cat("Null 95% interval:", round(quantile(res_dsc_after$dsc_null, c(0.025, 0.975)), 3), "\n")







par(mfrow = c(1,2))
hist(res_dsc_before$dsc_null,
     breaks = 40,
     main = "DSC Null (Before)",
     xlab = "DSC",
     col = "lightblue", border = "white")
abline(v = res_dsc_before$dsc_obs, col = "darkblue", lwd = 2)
legend("topright", legend = "Observed DSC", col = "darkblue", lwd = 2, bty = "n")

hist(res_dsc_after$dsc_null,
     breaks = 40,
     main = "DSC Null (After K=1)",
     xlab = "DSC",
     col = "pink", border = "white")
abline(v = res_dsc_after$dsc_obs, col = "red", lwd = 2)
legend("topright", legend = "Observed DSC", col = "red", lwd = 2, bty = "n")











##FIT:

#k=5

# 1. Center and scale data
Y <- scale(dat, center = TRUE, scale = TRUE)    # n x p

# 2. Load Stan fit and extract Lambda_hat (p x K)
fit_obj <- readRDS("fit_Joint_NHANES2017_phthalates_scale_all_randominit.rds")
Lambda_hat <- fit_obj$Lambda_hat
if (nrow(Lambda_hat) != ncol(Y)) Lambda_hat <- t(Lambda_hat)

# 3. Project out the fitted K=5 factors (Bayesian fit)
eta_hat <- Y %*% Lambda_hat %*% solve(t(Lambda_hat) %*% Lambda_hat)   # n x K
Y_hat   <- eta_hat %*% t(Lambda_hat)                                 # n x p
resid5  <- Y - Y_hat

# 4. Compute off-diagonal correlations
get_offdiag <- function(M) M[lower.tri(M)]
cor_before <- cor(Y)
cor_after  <- cor(resid5)
r_before <- get_offdiag(cor_before)
r_after  <- get_offdiag(cor_after)

# 5. Plot overlayed histograms
breaks <- seq(-1, 1, length.out = 41)
hist(r_before,
     breaks = breaks,
     col = rgb(0.2, 0.4, 1, 0.5),
     main = "Off-diagonal correlations\nbefore/after K=5 Bayesian fit",
     xlab = "Correlation",
     ylab = "Frequency",
     border = "white")
hist(r_after,
     breaks = breaks,
     col = rgb(1, 0.4, 0.4, 0.5),
     add = TRUE,
     border = "white")
legend("topright",
       legend = c("Before", "After K=5 Bayesian fit"),
       fill = c(rgb(0.2, 0.4, 1, 0.5), rgb(1, 0.4, 0.4, 0.5)),
       border = NA)

# 6. Print summary statistics
cat("Summary of off-diagonal correlations BEFORE K=5 Bayesian fit:\n")
print(summary(r_before))

cat("\nSummary of off-diagonal correlations AFTER K=5 Bayesian fit:\n")
print(summary(r_after))

cat("\nMean abs off-diagonal correlation BEFORE: ", mean(abs(r_before)), "\n")
cat("Mean abs off-diagonal correlation AFTER:  ", mean(abs(r_after)), "\n")





# --- K=1 ---

# 1. Load K=1 Stan fit
fit_obj1 <- readRDS("fit_Joint_NHANES1718_1.rds")
Lambda_hat1 <- fit_obj1$Lambda_hat
if (nrow(Lambda_hat1) != ncol(Y)) Lambda_hat1 <- t(Lambda_hat1)

# 2. Project out the fitted K=1 factor
eta_hat1 <- Y %*% Lambda_hat1 %*% solve(t(Lambda_hat1) %*% Lambda_hat1)   # n x 1
Y_hat1   <- eta_hat1 %*% t(Lambda_hat1)                                  # n x p
resid1   <- Y - Y_hat1

# 3. Compute off-diagonal correlations
cor_after1 <- cor(resid1)
r_after1 <- get_offdiag(cor_after1)

# 4. Plot overlayed histograms (before, after K=1, after K=5)
breaks <- seq(-1, 1, length.out = 41)
hist(r_before,
     breaks = breaks,
     col = rgb(0.2, 0.4, 1, 0.4),
     main = "Off-diagonal correlations\nBefore / K=1 / K=5",
     xlab = "Correlation",
     ylab = "Frequency",
     border = "white")
hist(r_after1,
     breaks = breaks,
     col = rgb(1, 0.7, 0.2, 0.5),    # orange
     add = TRUE,
     border = "white")
hist(r_after,
     breaks = breaks,
     col = rgb(1, 0.4, 0.4, 0.5),
     add = TRUE,
     border = "white")
legend("topright",
       legend = c("Before", "After K=1", "After K=5"),
       fill = c(rgb(0.2, 0.4, 1, 0.4), rgb(1, 0.7, 0.2, 0.5), rgb(1, 0.4, 0.4, 0.5)),
       border = NA)

# 5. Print summary statistics for K=1 as well
cat("\nSummary of off-diagonal correlations AFTER K=1 Bayesian fit:\n")
print(summary(r_after1))

cat("\nMean abs off-diagonal correlation AFTER K=1:  ", mean(abs(r_after1)), "\n")








