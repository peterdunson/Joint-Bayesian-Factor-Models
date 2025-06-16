# simulation_data_tebfar.R
# Simulate data for TEB-FAR/Bayesian infinite factor regression models
set.seed(123)   # For reproducibility

# ----- Simulation Parameters -----
n_train <- 200   # Number of training samples (adjust as needed)
n_test  <- 1000  # Number of test samples (for completeness)
p <- 20          # Number of predictors (excluding outcome)
K <- 10          # Number of latent factors (matches factor model truncation)

# ----- Construct Lambda -----
Lambda <- matrix(0, nrow = p, ncol = K)
nonzero_each_col <- 10
for (i in 1:K) {
   idx <- sample(1:p, nonzero_each_col)
   Lambda[idx, i] <- rexp(nonzero_each_col, rate=1)
   Lambda[,i] <- Lambda[,i] / sqrt(sum(Lambda[,i]^2)) * (i/K)
}
# Set true regression effect in the last factor only
lambda_y <- c(1, rep(0, K-1))
Lambda <- rbind(lambda_y, Lambda)    # First row is the outcome loading

# ----- Noise Covariance -----
Sigma <- diag(rep(0.2, p + 1))

# ----- Covariance Matrix -----
Omega <- Lambda %*% t(Lambda) + Sigma

# Normalize Omega for numerical stability (as in original script)
for (i in 1:nrow(Omega)) {
   for (j in 1:ncol(Omega)) {
      Omega[i,j] <- Omega[i,j] / sqrt(Omega[i,i] * Omega[j,j])
   }
}

# ----- Simulate Data -----
library(mvtnorm)
yX <- rmvnorm(n = n_train, sigma = Omega)  # [n x (p+1)], first column is outcome y

# Save for Stan
saveRDS(
   list(
      Y = yX,            # Data matrix [n x (p+1)], first column is y
      Omega = Omega,     # Covariance used
      Lambda = Lambda    # True factor loadings
   ),
   file = "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen1.rds"
)
