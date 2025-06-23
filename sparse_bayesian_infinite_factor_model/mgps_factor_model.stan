// mgps_factor_model_nc.stan
// Non-centered MGPS infinite factor model

data {
  int<lower=1> N;           // number of observations
  int<lower=1> P;           // number of variables
  int<lower=1> K;           // truncation level
  matrix[N, P] Y;           // centered & scaled data
}

parameters {
  // Non-centered loadings
  matrix[P, K] Lambda_raw;  
  matrix[N, K] eta;            
  vector<lower=0>[P] psi;       
  matrix<lower=0>[P, K] phi;    
  vector<lower=0>[K] delta;     
}

transformed parameters {
  vector<lower=0>[K] tau;
  matrix[P, K] lambda_sd;
  matrix[P, K] Lambda;         // *centered* loadings

  // build MGPS global weights
  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k-1] * delta[k];

  // per-loading SD and re-scale raw loadings
  for (p in 1:P) {
    for (k in 1:K) {
      lambda_sd[p,k] = sqrt(1.0 / (phi[p,k] * tau[k]));
      Lambda[p,k]   = Lambda_raw[p,k] * lambda_sd[p,k];
    }
  }
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
  to_vector(phi) ~ gamma(nu/2, nu/2);
  delta[1] ~ gamma(a1, b1);
  for (k in 2:K)
    delta[k] ~ gamma(a2, b2);

  // non-centered loading prior
  to_vector(Lambda_raw) ~ normal(0, 1);

  // factor scores
  to_vector(eta) ~ normal(0, 1);

  // likelihood
  {
    matrix[N, P] mu = eta * Lambda';
    for (p in 1:P)
      Y[,p] ~ normal(mu[,p], sqrt(1.0 / psi[p]));
  }
}





