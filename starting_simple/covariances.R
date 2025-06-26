# Load Simulation 2
sim <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Y <- as.matrix(sim$Y)
n <- nrow(Y)
p <- ncol(Y)

B <- 100  # Number of bootstraps
set.seed(42)

# To store lambda estimates per bootstrap
lambda_mat <- matrix(NA, nrow=B, ncol=p)

compat_fraction <- numeric(B)
all_cpr <- numeric()

for (b in 1:B) {
   idx <- sample(1:n, replace=TRUE)
   Yb <- Y[idx, ]
   Cb <- cov(Yb)
   
   # Compatibility: for all triplets, is the sign product positive?
   triplet_compatible <- c()
   for (trip in combn(p, 3, simplify=FALSE)) {
      cpr <- Cb[trip[1], trip[2]]
      cps <- Cb[trip[1], trip[3]]
      crs <- Cb[trip[2], trip[3]]
      triplet_compatible <- c(triplet_compatible, (cpr * cps * crs >= 0))
   }
   compat_fraction[b] <- mean(triplet_compatible)
   
   # Estimate lambdas via method of moments, pairwise: c_{pr} = lambda_p * lambda_r
   # 1. Fix lambda_1 > 0 by convention
   # 2. For r > 1, lambda_r = c_{1r} / lambda_1
   # 3. lambda_1^2 = c_{11} - psi_1 (estimate psi_1 as c_{11} - max(abs(c_{1, -1})))
   # Here, just solve system for lambda_1 using all off-diagonals
   # We'll use absolute value for numerical stability.
   lambda1_vals <- sqrt(abs(Cb[1, -1]))
   lambda_1 <- mean(lambda1_vals) # take mean of sqrt(abs(c_{1r})) for r != 1
   # Now estimate rest
   lambda_vec <- numeric(p)
   lambda_vec[1] <- lambda_1
   for (j in 2:p) {
      lambda_vec[j] <- ifelse(lambda_1 == 0, 0, Cb[1, j] / lambda_1)
   }
   lambda_mat[b, ] <- lambda_vec
   
   # Store all off-diagonal c_{pr}
   all_cpr <- c(all_cpr, Cb[lower.tri(Cb)])
}

# Plot: histogram of fraction of compatible triplets per bootstrap
hist(compat_fraction, breaks=30, main="Fraction of compatible triplets per bootstrap", xlab="Fraction compatible")

cat("Mean fraction of compatible triplets per bootstrap:", mean(compat_fraction), "\n")

# Plot: histogram of all off-diagonal c_{pr} (all bootstraps)
hist(all_cpr, breaks=30, main="Distribution of Off-Diagonal Covariances (All Bootstraps)", xlab="Off-diagonal Covariance")

# Plot: lambda estimates distribution for variable 1 and variable 2
hist(lambda_mat[,1], breaks=30, main="Lambda_1 estimates (method of moments)", xlab="Estimated lambda_1")
hist(lambda_mat[,2], breaks=30, main="Lambda_2 estimates (method of moments)", xlab="Estimated lambda_2")




# Compute the sample correlation matrix
R_obs <- cor(Y)

# Get all off-diagonal correlations
off_diag_corrs <- R_obs[lower.tri(R_obs)]

# Plot the histogram
hist(off_diag_corrs, breaks=30, main="Distribution of Pairwise Correlations", xlab="Correlation")

# Optional: summary stats
summary(off_diag_corrs)

