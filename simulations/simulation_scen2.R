# simulation_scen2.R
# Simulate data for Scenario 2 for TEB-FAR/Bayesian infinite factor regression
set.seed(123)   # For reproducibility

# ----- Parameters -----
n_train <- 1000   # Training samples
p       <- 20     # Predictors (not counting outcome)
K       <- 5      # Latent factors

# ----- Lambda: each col gets nonzeros in random rows (including y) -----
nonzero_each_col <- 10
Lambda <- matrix(0, nrow = p + 1, ncol = K)
for (i in seq_len(K)) {
   idx <- sample.int(p + 1, nonzero_each_col)
   Lambda[idx, i] <- rexp(nonzero_each_col, rate = 1)
   # normalize to length i/K
   Lambda[, i] <- Lambda[, i] / sqrt(sum(Lambda[, i]^2)) * (i / K)
}

# ----- Noise Covariance -----
Sigma <- diag(rep(0.2, p + 1))

# ----- Covariance Matrix (for checking) -----
Omega <- Lambda %*% t(Lambda) + Sigma
# Normalize to correlation scale
for (i in seq_len(nrow(Omega))) {
   for (j in seq_len(ncol(Omega))) {
      Omega[i, j] <- Omega[i, j] / sqrt(Omega[i, i] * Omega[j, j])
   }
}

# ----- Simulate latent factors and data -----
# 1) latent factors η: n_train × K
eta <- matrix(rnorm(n_train * K), nrow = n_train, ncol = K)

# 2) noise ε: n_train × (p+1)
library(mvtnorm)
epsilon <- rmvnorm(n = n_train, sigma = Sigma)

# 3) observed data yX: n_train × (p+1)
#    first column of yX corresponds to the 'y' variable
yX <- eta %*% t(Lambda) + epsilon

# ----- Save for Stan/TEB-FAR -----
saveRDS(
   list(
      Y      = yX,      # Data matrix [n × (p+1)], col 1 is y
      Omega  = Omega,   # True correlation matrix
      Lambda = Lambda,  # True loadings
      eta    = eta      # True latent factors [n × K]
   ),
   file = "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds"
)

