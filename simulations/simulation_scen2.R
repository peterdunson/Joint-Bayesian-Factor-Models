# simulation_scen2.R
# Simulate data for Scenario 2 for TEB-FAR/Bayesian infinite factor regression
set.seed(123)   # For reproducibility

# ----- Parameters -----
n_train <- 1000   # Training samples
p <- 20          # Predictors (not counting outcome)
K <- 5         # Latent factors

# ----- Lambda: each col gets nonzeros in random rows (including y) -----
nonzero_each_col <- 10
Lambda <- matrix(0, nrow = p + 1, ncol = K)
for (i in 1:K) {
   idx <- sample(1:(p+1), nonzero_each_col)
   Lambda[idx, i] <- rexp(nonzero_each_col, rate=1)
   Lambda[, i] <- Lambda[, i] / sqrt(sum(Lambda[, i]^2)) * (i/K)
}

# ----- Noise Covariance -----
Sigma <- diag(rep(0.2, p + 1))

# ----- Covariance Matrix -----
Omega <- Lambda %*% t(Lambda) + Sigma
# Normalize for numerical stability
for (i in 1:nrow(Omega)) {
   for (j in 1:ncol(Omega)) {
      Omega[i, j] <- Omega[i, j] / sqrt(Omega[i, i] * Omega[j, j])
   }
}

# ----- Simulate data -----
library(mvtnorm)
yX <- rmvnorm(n = n_train, sigma = Omega)  # [n x (p+1)], col 1 is y

# ----- Save for Stan/TEB-FAR -----
saveRDS(
   list(
      Y = yX,            # Data matrix [n x (p+1)], first col is y
      Omega = Omega,     # True covariance
      Lambda = Lambda    # True loadings
   ),
   file = "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds"
)
