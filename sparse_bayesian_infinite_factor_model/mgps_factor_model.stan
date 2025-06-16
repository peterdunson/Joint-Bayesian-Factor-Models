// mgps_factor_model.stan
// Bhattacharya & Dunson (2011) Sparse Bayesian Infinite Factor Model

data {
  int<lower=1> N;      // number of observations
  int<lower=1> P;      // number of observed variables
  int<lower=1> K;      // upper bound on number of factors (set large, e.g. 20)
  matrix[N, P] Y;      // data matrix (centered and scaled)
}

parameters {
  // Factor loadings
  matrix[P, K] Lambda;
  // Factor scores
  matrix[N, K] eta;
  // Residual precisions (psi)
  vector<lower=0>[P] psi;
  // Local shrinkage
  matrix<lower=0>[P, K] phi;
  // Multiplicative gamma process: delta (not tau!)
  vector<lower=0>[K] delta;
}

transformed parameters {
  // Compute tau as the cumulative product of delta
  vector<lower=0>[K] tau;
  matrix[P, K] lambda_sd;
  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k-1] * delta[k];
  for (k in 1:K) {
    for (p in 1:P) {
      lambda_sd[p, k] = sqrt(1.0 / (phi[p, k] * tau[k]));
    }
  }
}

model {
  // --- Hyperparameters (use the paper's defaults) ---
  real a1 = 2.1;
  real a2 = 3.1;
  real b1 = 1.0;
  real b2 = 1.0;
  real nu = 3.0;   // df for local shrinkage
  real a_psi = 1.0; 
  real b_psi = 0.3;

  // --- Priors ---
  psi ~ gamma(a_psi, b_psi);

  // Local shrinkage (phi)
  to_vector(phi) ~ gamma(nu / 2.0, nu / 2.0);

  // First delta
  delta[1] ~ gamma(a1, b1);
  // Remaining deltas
  for (k in 2:K)
    delta[k] ~ gamma(a2, b2);

  // Factor loadings
  for (p in 1:P)
    for (k in 1:K)
      Lambda[p, k] ~ normal(0, lambda_sd[p, k]);

  // Factor scores
  for (n in 1:N)
    eta[n] ~ normal(0, 1);

  // --- Likelihood ---
  for (n in 1:N) {
    vector[P] mu = Lambda * eta[n]';
    for (p in 1:P)
      Y[n, p] ~ normal(mu[p], sqrt(1.0 / psi[p]));
  }
}

