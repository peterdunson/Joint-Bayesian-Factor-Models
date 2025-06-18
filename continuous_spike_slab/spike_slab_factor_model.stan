// spike_slab_factor_model.stan
// ----------------------------------------------------------
// Sparse Bayesian factor analysis with a continuous spike-and-slab 
// (Laplace mixture) prior on loadings Lambda, following:
//   - Rocková (2018), Annals of Statistics, for the prior
//   - Bhattacharya & Dunson (2011), for factor model structure
// Each Lambda[p, k] ~ (1-theta) Laplace(0, 1/lambda0) + theta Laplace(0, 1/lambda1)
// The "spike" (lambda0) is strong shrinkage; the "slab" (lambda1) is weak shrinkage
// No discrete indicators are used (Stan limitation); mixture is marginalized.
// ----------------------------------------------------------

functions {
  // Laplace (double-exponential) log density
  real laplace_lpdf(real x, real mu, real b) {
    return -log(2 * b) - fabs(x - mu) / b;
  }
}
data {
  int<lower=1> N;           // number of observations
  int<lower=1> P;           // number of observed variables
  int<lower=1> K;           // number of factors (truncation level)
  matrix[N, P] Y;           // centered data (N x P)
  real<lower=0> lambda0;    // spike penalty (e.g., 20)
  real<lower=0> lambda1;    // slab penalty (e.g., 0.2)
  real<lower=0,upper=1> theta; // prior inclusion probability (e.g., 0.1)
}
parameters {
  matrix[P, K] Lambda;      // factor loadings
  matrix[N, K] eta;         // factor scores
  vector<lower=0>[P] psi;   // residual precisions
}
model {
  // Prior: spike-and-slab lasso for each loading
  for (p in 1:P)
    for (k in 1:K)
      target += log_mix(
        theta,
        laplace_lpdf(Lambda[p, k] | 0, 1/lambda1),   // slab (weak shrinkage)
        laplace_lpdf(Lambda[p, k] | 0, 1/lambda0));  // spike (strong shrinkage)

  // Prior for factor scores
  to_vector(eta) ~ normal(0, 1);

  // Prior for residual precisions (as in Bhattacharya & Dunson 2011)
  real a_psi = 1.0;
  real b_psi = 0.3;
  psi ~ gamma(a_psi, b_psi);

  // Likelihood
  {
    matrix[N, P] mu = eta * Lambda';
    for (p in 1:P)
      Y[, p] ~ normal(mu[, p], sqrt(1.0 / psi[p]));
  }
}

