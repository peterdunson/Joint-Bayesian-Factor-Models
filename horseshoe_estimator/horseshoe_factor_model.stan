// horseshoe_factor_model.stan
// Sparse Factor Model with Horseshoe Prior (all features symmetric)

data {
  int<lower=1> N;         // number of observations
  int<lower=1> P;         // number of variables
  int<lower=1> K;         // number of factors
  matrix[N, P] Y;         // centered data
}

parameters {
  matrix[P, K] Lambda;              // factor loadings
  matrix[N, K] eta;                 // factor scores

  vector<lower=0>[P] psi;           // residual precisions

  real<lower=0> lambda_global;           // global shrinkage for loadings
  matrix<lower=0>[P, K] lambda_local;    // local shrinkage for loadings
}

model {
  // Hyperparameters for residual variances
  real a_psi = 1.0;
  real b_psi = 0.3;

  // Priors for residual variances
  psi ~ gamma(a_psi, b_psi);

  // Horseshoe priors
  lambda_global ~ cauchy(0, 1);
  to_vector(lambda_local) ~ cauchy(0, 1);

  // Prior for Lambda (horseshoe)
  for (p in 1:P)
    for (k in 1:K)
      Lambda[p, k] ~ normal(0, lambda_global * lambda_local[p, k]);

  // Prior for factors
  to_vector(eta) ~ normal(0, 1);

  // Likelihood
  {
    matrix[N, P] mu = eta * Lambda';
    for (p in 1:P)
      Y[, p] ~ normal(mu[, p], sqrt(1.0 / psi[p]));
  }
}



