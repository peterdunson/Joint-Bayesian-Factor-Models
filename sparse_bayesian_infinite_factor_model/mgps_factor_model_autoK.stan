// mgps_factor_model_autoK.stan

data {
  int<lower=1> N;          // number of observations
  int<lower=1> P;          // number of variables
  int<lower=1> K_max;      // truncation for # factors
  matrix[N, P] Y;          // data matrix
}

parameters {
  // non-centered loadings
  matrix[P, K_max] raw_Lambda;
  // factor scores
  matrix[N, K_max] eta;
  // residual precisions
  vector<lower=0>[P] psi;
  // local shrinkage
  matrix<lower=0>[P, K_max] phi;
  // global shrinkage multipliers
  vector<lower=0>[K_max] delta;
}

transformed parameters {
  vector<lower=0>[K_max] tau;
  matrix[P, K_max]       lambda_sd;
  matrix[P, K_max]       Lambda;

  // Cumulative shrinkage (multiplicative gamma process)
  tau[1] = delta[1];
  for (k in 2:K_max)
    tau[k] = tau[k - 1] * delta[k];

  // Local–global scale for each loading
  for (p in 1:P)
    for (k in 1:K_max) {
      lambda_sd[p, k] = sqrt(1.0 / (phi[p, k] * tau[k]));
      Lambda[p, k]    = raw_Lambda[p, k] * lambda_sd[p, k];
    }
}

model {
  // Hyperpriors (tweak these if you like)
  psi             ~ gamma(2.0, 1.0);       // residual precisions
  to_vector(phi)  ~ gamma(2.5, 2.5);       // local shrinkages
  delta[1]        ~ gamma(3.0, 1.0);
  for (k in 2:K_max)
    delta[k]      ~ gamma(4.0, 1.0);

  // Non-centered loadings and scores
  to_vector(raw_Lambda) ~ normal(0, 1);
  to_vector(eta)        ~ normal(0, 1);

  // Likelihood
  {
    matrix[N, P] mu = eta * Lambda';
    for (p in 1:P)
      Y[, p] ~ normal(mu[, p], sqrt(1.0 / psi[p]));
  }
}


