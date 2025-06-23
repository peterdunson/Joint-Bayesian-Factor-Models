// mgps_factor_model.stan
// Sparse Bayesian infinite factor model (Bhattacharya & Dunson, 2011)

data {
  int<lower=1> N;           // number of observations
  int<lower=1> P;           // number of variables
  int<lower=1> K;           // truncation level
  matrix[N, P] Y;           // centered & scaled data matrix
}

parameters {
  matrix[P, K] Lambda;      // factor loadings
  matrix[N, K] eta;         // latent factor scores
  vector<lower=0>[P] psi;   // residual precisions (1 / var)
  matrix<lower=0>[P, K] phi;// local shrinkage scales
  vector<lower=0>[K] delta; // gamma-process weights
}

transformed parameters {
  vector<lower=0>[K] tau;   // global shrinkage for each factor index
  matrix[P, K] lambda_sd;   // per‐loading SD = (phi * tau)^{-1/2}

  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k - 1] * delta[k];

  for (p in 1:P)
    for (k in 1:K)
      lambda_sd[p, k] = sqrt(1.0 / (phi[p, k] * tau[k]));
}

model {
  // hyperparameters (paper defaults)
  real a_psi = 1.0;
  real b_psi = 0.3;
  real nu    = 3.0;
  real a1    = 2.1;
  real b1    = 1.0;
  real a2    = 3.1;
  real b2    = 1.0;

  // priors
  psi  ~ gamma(a_psi, b_psi);
  to_vector(phi) ~ gamma(nu / 2, nu / 2);
  delta[1]      ~ gamma(a1, b1);
  for (k in 2:K)
    delta[k]    ~ gamma(a2, b2);

  // centered loading prior
  for (p in 1:P)
    for (k in 1:K)
      Lambda[p, k] ~ normal(0, lambda_sd[p, k]);

  // factor scores
  to_vector(eta) ~ normal(0, 1);

  // likelihood
  {
    matrix[N, P] mu = eta * Lambda';
    for (p in 1:P)
      Y[, p] ~ normal(mu[, p], sqrt(1.0 / psi[p]));
  }
}
