// mgps_factor_model.stan
// Sparse Bayesian infinite factor model (Bhattacharya & Dunson, 2011)
// with multiplicative gamma process shrinkage (MGPS) on loadings.

data {
  int<lower=1> N;           // number of observations
  int<lower=1> P;           // number of variables
  int<lower=1> K;           // truncation level (choose large, e.g. 30)
  matrix[N, P] Y;           // centered & scaled data matrix
}

parameters {
  matrix[P, K] Lambda;      // factor loadings
  matrix[N, K]  eta;        // latent factor scores
  vector<lower=0>[P] psi;   // residual precisions (1 / var)
  matrix<lower=0>[P, K] phi;  // local shrinkage scales
  vector<lower=0>[K] delta;   // gamma-process weights
}

transformed parameters {
  vector<lower=0>[K] tau;   // global shrinkage for each factor index
  matrix[P, K] lambda_sd;   // per‐loading SD = (phi * tau)^{-1/2}

  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k - 1] * delta[k];

  for (k in 1:K)
    for (p in 1:P)
      lambda_sd[p, k] = sqrt(1.0 / (phi[p, k] * tau[k]));
}

model {
  // ------------- Hyperparameters (paper defaults) -------------
  real a_psi = 1.0;    // residual precision ~ Gamma(a_psi, b_psi)
  real b_psi = 0.3;
  real nu    = 3.0;    // local shrinkage df
  real a1    = 2.1;    // delta[1] ~ Gamma(a1, 1)
  real b1    = 1.0;
  real a2    = 3.1;    // delta[k>1] ~ Gamma(a2, 1)
  real b2    = 1.0;

  // ------- Priors -------
  psi  ~ gamma(a_psi, b_psi);
  to_vector(phi) ~ gamma(nu / 2, nu / 2);

  delta[1]      ~ gamma(a1, b1);
  for (k in 2:K)
    delta[k]    ~ gamma(a2, b2);

  // Factor loadings with MGPS‐induced SDs
  for (p in 1:P)
    for (k in 1:K)
      Lambda[p, k] ~ normal(0, lambda_sd[p, k]);

  // Factor scores
  to_vector(eta) ~ normal(0, 1);

  // ------- Likelihood -------
  {
    // vectorize over n and p for efficiency
    matrix[N, P] mu = eta * Lambda'; 
    for (p in 1:P)
      Y[, p] ~ normal(mu[, p], sqrt(1.0 / psi[p]));
  }
}

