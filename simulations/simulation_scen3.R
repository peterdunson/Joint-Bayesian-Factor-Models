# simulation_scen3.R
# Scenario 3: X from factor model, y weakly related to X (not from factors)

set.seed(123)

# -- Parameters --
n_train <- 1000
p <- 20
k <- 8
rsquared <- 0.1      # Strength of y~X
nonzero_coefs <- 12  # Number of nonzero betas

# -- Simulate X (factor model) --
Lambda <- matrix(rexp(p * k, 1), nrow = p, ncol = k)
Sigma <- diag(rep(0.01, p))
Omega <- Lambda %*% t(Lambda) + Sigma
for (i in 1:p) for (j in 1:p) {
   Omega[i, j] <- Omega[i, j] / sqrt(Omega[i, i] * Omega[j, j])
}

# -- Sparse regression coefficients --
set.seed(42)
nonzero_idx <- sample(1:p, nonzero_coefs)
coefs <- rep(0, p)
coefs[nonzero_idx] <- rnorm(nonzero_coefs)

# -- Normalize for desired marginal R^2 --
unnorm_var_y <- as.numeric(t(coefs) %*% Omega %*% coefs)
beta <- coefs / sqrt(unnorm_var_y) * sqrt(rsquared)

# -- Generate training data --
library(mvtnorm)
X <- rmvnorm(n_train, sigma = Omega)
epsilon <- rnorm(n_train, sd = sqrt(1 - rsquared))
y <- as.vector(X %*% beta) + epsilon
yX <- cbind(y, X)   # [n_train x (p+1)], col 1 = y

# -- Save for Stan/TEB-FAR --
saveRDS(
   list(
      Y = yX,            # [n_train x (p+1)], col 1 = y
      X = X,
      y = y,
      Omega = Omega,
      Lambda = Lambda,
      beta = beta
   ),
   file = "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen3_1000.rds"
)

