// horseshoe_factor_model.stan
// Sparse Factor Model with Horseshoe Prior (Carvalho et al. 2010)

data {
  int<lower=1> N;         // number of observations
  int<lower=1> P;         // number of variables (first column is outcome)
  int<lower=1> K;         // number of factors (truncation)
  matrix[N, P] Y;         // centered data
  real<lower=0> Sigma1;   // fixed variance for outcome (psi[1])
}

parameters {
  matrix[P, K] Lambda;              // factor loadings
  matrix[N, K] eta;                 // factor scores

  // Residual precisions
  vector<lower=0>[P-1] psi_free;

  // Horseshoe global and local scales
  real<lower=0> tau_global;              // global shrinkage for loadings
  matrix<lower=0>[P, K] tau_local;       // local shrinkage for loadings

  // Horseshoe auxiliary variables for the half-Cauchy priors
  real<lower=0> aux_global;
  matrix<lower=0>[P, K] aux_local;
}

transformed parameters {
  vector<lower=0>[P] psi;
  psi[1] = 1.0 / Sigma1;
  for (j in 2:P)
    psi[j] = psi_free[j-1];
}

model {
  // Hyperparameters for residual variances (same as TEB-FAR for comparability)
  real a_psi = 1.0;
  real b_psi = 0.3;

  // Priors for residual variances
  psi_free ~ gamma(a_psi, b_psi);

  // Horseshoe priors for Lambda (per Carvalho et al. 2010, Section 3.1):
  // Lambda_{jk} ~ N(0, tau_global^2 * tau_local[j,k]^2)

  // Prior: tau_global ~ half-Cauchy(0, 1)
  tau_global ~ normal(0, 1);
  aux_global ~ normal(0, 1);
  target += -log(1 + square(tau_global / aux_global)); // Jacobian adjustment for half-Cauchy

  // Prior: tau_local[j,k] ~ half-Cauchy(0, 1)
  to_vector(tau_local) ~ normal(0, 1);
  to_vector(aux_local) ~ normal(0, 1);
  target += sum(-log(1 + square(to_vector(tau_local) ./ to_vector(aux_local)))); // Jacobian

  // Prior for Lambda
  for (p in 1:P)
    for (k in 1:K)
      Lambda[p, k] ~ normal(0, tau_global * tau_local[p, k]);

  // Prior for factors
  to_vector(eta) ~ normal(0, 1);

  // Likelihood
  {
    matrix[N, P] mu = eta * Lambda';
    for (p in 1:P)
      Y[, p] ~ normal(mu[, p], sqrt(1.0 / psi[p]));
  }
}

