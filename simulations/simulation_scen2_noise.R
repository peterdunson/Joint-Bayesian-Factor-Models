# simulation_scen2_noise.R
# Simulate *pure noise* data for Scenario 2 for TEB-FAR/Bayesian infinite factor regression
set.seed(123)   # For reproducibility

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations")

# ----- Parameters -----
n_train <- 1000   # Training samples
p       <- 20     # Predictors (not counting outcome)
K       <- 1      # Latent factors (ignored in noise sim, for compatibility)

# ----- Lambda: random, unused (keep structure for compatibility) -----
nonzero_each_col <- 10
Lambda <- matrix(0, nrow = p + 1, ncol = K)
for (i in seq_len(K)) {
   idx <- sample.int(p + 1, nonzero_each_col)
   Lambda[idx, i] <- rexp(nonzero_each_col, rate = 1)
   Lambda[, i] <- Lambda[, i] / sqrt(sum(Lambda[, i]^2)) * (i / K)
}

# ----- Noise Covariance -----
Sigma <- diag(rep(5, p + 1))

# ----- Omega: just Sigma on correlation scale -----
Omega <- Sigma
for (i in seq_len(nrow(Omega))) {
   for (j in seq_len(ncol(Omega))) {
      Omega[i, j] <- Omega[i, j] / sqrt(Omega[i, i] * Omega[j, j])
   }
}

# ----- Simulate data: just noise -----
library(mvtnorm)
eta <- matrix(NA, nrow = n_train, ncol = K)  # No true factors
yX  <- rmvnorm(n = n_train, sigma = Sigma)   # Just noise

# ----- Save for Stan/TEB-FAR -----
saveRDS(
   list(
      Y      = yX,      # Data matrix [n × (p+1)], col 1 is y
      Omega  = Omega,   # True correlation matrix (just noise)
      Lambda = Lambda,  # Loadings (random, unused here)
      eta    = eta      # No true factors
   ),
   file = "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000_noise.rds"
)

