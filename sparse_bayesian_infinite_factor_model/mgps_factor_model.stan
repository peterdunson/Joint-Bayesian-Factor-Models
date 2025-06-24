// mgsp_factor_model_tuned2.stan


data {
  int<lower=1> N;         // # observations
  int<lower=1> P;         // # variables
  int<lower=1> K;         // number of factors
  matrix[N, P] Y;         // centered & scaled data
}

parameters {
  matrix[P, K] raw_Lambda;     // non‐centered loadings
  matrix[N, K] raw_eta;        // non‐centered scores
  vector[P]    log_psi;        // log‐residual precisions
  matrix<lower=0>[P, K] phi;   // local shrinkage
  vector<lower=0>[K] delta;    // global shrinkage
}

transformed parameters {
  vector<lower=0>[P]    psi = exp(log_psi);
  vector<lower=0>[K]    tau;
  matrix[P, K]          lambda_sd;
  matrix[P, K]          Lambda;
  matrix[N, K]          eta = raw_eta;

  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k-1] * delta[k];

  for (p in 1:P)
    for (k in 1:K) {
      lambda_sd[p,k] = sqrt(1.0 / (phi[p,k] * tau[k]));
      Lambda[p,k]    = raw_Lambda[p,k] * lambda_sd[p,k];
    }
}

model {
  // priors (unchanged)
  log_psi         ~ normal(0, 1);         
  to_vector(phi)  ~ gamma(2.5, 2.5);
  delta[1]        ~ gamma(3.0, 1.0);
  for (k in 2:K)
    delta[k]      ~ gamma(4.0, 1.0);

  to_vector(raw_Lambda) ~ normal(0, 1);
  to_vector(raw_eta)    ~ normal(0, 1);

  // vectorized likelihood
  {
    matrix[N,P] mu = eta * Lambda';
    for (n in 1:N)
      Y[n] ~ normal(mu[n], 1.0 / sqrt(psi));
  }
}




