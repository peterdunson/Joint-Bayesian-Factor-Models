# Load the fit object and simulation data
fit <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment/fit_Joint_scen2_scale_all_randominit.rds")
sim <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")

# Center and scale the original data the same way as in model fitting
Y <- scale(sim$Y, center = TRUE, scale = TRUE)

# Extract Lambda_hat (p x K), and get n, p from Y
Lambda_hat <- fit$Lambda_hat
n <- nrow(Y)
p <- ncol(Y)
K <- ncol(Lambda_hat)

# Estimate factor scores
eta_hat <- Y %*% Lambda_hat %*% solve(t(Lambda_hat) %*% Lambda_hat)  # n x K

# Reconstruct fitted values
Y_hat <- eta_hat %*% t(Lambda_hat)  # n x p

# Compute residuals
resid_mat <- Y - Y_hat

# Compute correlation matrices
R_obs <- cor(Y)
R_resid <- cor(resid_mat)

# Get all off-diagonal correlations
off_diag_corrs_obs <- R_obs[lower.tri(R_obs)]
off_diag_corrs_resid <- R_resid[lower.tri(R_resid)]

# === Fisher z-transform ===
fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
z_obs <- fisher_z(off_diag_corrs_obs)
z_resid <- fisher_z(off_diag_corrs_resid)

# Plot histograms
par(mfrow=c(1,2))
hist(z_obs, breaks=30, main="Fisher z: Observed", xlab="z(correlation)")
hist(z_resid, breaks=30, main="Fisher z: Residual", xlab="z(correlation)")
par(mfrow=c(1,1))

# Print summary statistics
cat("Mean abs Fisher-z obs:", mean(abs(z_obs)), "\n")
cat("Mean abs Fisher-z resid:", mean(abs(z_resid)), "\n")

# --- DSC distance: mean, SD, skewness, kurtosis (fourth moment) ---
dsc_skew_kurt <- function(corrs, mu2 = 0, sd2 = NULL, sk2 = 0, ku2 = 0, n = NULL) {
   mu1 <- mean(corrs)
   sd1 <- sd(corrs)
   # Avoid division by zero in skewness and kurtosis if sd1 is tiny
   if (sd1 < 1e-12) {
      sk1 <- 0
      ku1 <- 0
   } else {
      sk1 <- mean((corrs - mu1)^3) / sd1^3
      ku1 <- mean((corrs - mu1)^4) / sd1^4
   }
   if (is.null(sd2)) {
      if (is.null(n)) stop("Need n (sample size) if sd2 is not given")
      sd2 <- 1 / sqrt(n - 3)
   }
   sk_diff <- (sign(sk1) * abs(sk1)^(1/3)) - (sign(sk2) * abs(sk2)^(1/3))
   ku_diff <- (sign(ku1) * abs(ku1)^(1/4)) - (sign(ku2) * abs(ku2)^(1/4))
   D <- sqrt((mu1 - mu2)^2 + (sd1 - sd2)^2 + sk_diff^2 + ku_diff^2)
   return(list(
      DSC = D, mu1 = mu1, sd1 = sd1, sk1 = sk1, ku1 = ku1,
      mu2 = mu2, sd2 = sd2, sk2 = sk2, ku2 = ku2
   ))
}

# --- Compute and print DSC with skewness and kurtosis ---
dist_sk_kurt_obs <- dsc_skew_kurt(z_obs, n = n)
dist_sk_kurt_resid <- dsc_skew_kurt(z_resid, n = n)

cat(sprintf(
   "DSC (mean, sd, skew, kurtosis) for observed [Fisher-z]: %.4f (mean=%.4f, sd=%.4f, skew=%.4f, kurtosis=%.4f)\n",
   dist_sk_kurt_obs$DSC, dist_sk_kurt_obs$mu1, dist_sk_kurt_obs$sd1, dist_sk_kurt_obs$sk1, dist_sk_kurt_obs$ku1
))
cat(sprintf(
   "DSC (mean, sd, skew, kurtosis) for residuals [Fisher-z]: %.4f (mean=%.4f, sd=%.4f, skew=%.4f, kurtosis=%.4f)\n",
   dist_sk_kurt_resid$DSC, dist_sk_kurt_resid$mu1, dist_sk_kurt_resid$sd1, dist_sk_kurt_resid$sk1, dist_sk_kurt_resid$ku1
))
