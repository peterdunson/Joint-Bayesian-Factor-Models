// tebfar_factor_model.stan
// Targeted Empirical Bayes Factor Analysis Regression (TEB-FAR)
// Follows Palmer & Dunson (2025+), extends Bhattacharya & Dunson (2011)

data {
  int<lower=1> N;           // number of observations
  int<lower=1> P;           // number of variables (first column is outcome)
  int<lower=1> K;           // truncation level
  matrix[N, P] Y;           // centered data matrix (first column = outcome)
  real<lower=0> Sigma1;     // user-provided fixed variance for outcome (psi[1])
}

parameters {
  matrix[P, K] Lambda;             // factor loadings
  matrix[N, K] eta;                // factor scores
  vector<lower=0>[P-1] psi_free;   // residual precisions for predictors (free)
  matrix<lower=0>[P, K] phi;       // local shrinkage
  vector<lower=0>[K] delta;        // gamma-process weights
}

transformed parameters {
  vector<lower=0>[K] tau;
  matrix[P, K] lambda_sd;
  vector<lower=0>[P] psi;

  // psi[1] is fixed; rest come from psi_free
  psi[1] = 1.0 / Sigma1;
  for (j in 2:P)
    psi[j] = psi_free[j-1];

  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k-1] * delta[k];

  for (k in 1:K)
    for (p in 1:P)
      lambda_sd[p, k] = sqrt(1.0 / (phi[p, k] * tau[k]));
}

model {
  // Priors (per paper)
  real a_psi = 1.0;
  real b_psi = 0.3;
  real nu    = 3.0;
  real a1    = 2.1;
  real b1    = 1.0;
  real a2    = 3.1;
  real b2    = 1.0;

  psi_free ~ gamma(a_psi, b_psi);
  to_vector(phi) ~ gamma(nu / 2, nu / 2);

  delta[1] ~ gamma(a1, b1);
  for (k in 2:K)
    delta[k] ~ gamma(a2, b2);

  for (p in 1:P)
    for (k in 1:K)
      Lambda[p, k] ~ normal(0, lambda_sd[p, k]);

  to_vector(eta) ~ normal(0, 1);

  // Likelihood (uses fixed psi[1] for outcome)
  {
    matrix[N, P] mu = eta * Lambda';
    for (p in 1:P)
      Y[, p] ~ normal(mu[, p], sqrt(1.0 / psi[p]));
  }
}
