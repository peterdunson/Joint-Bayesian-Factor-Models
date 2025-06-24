// mgsp_factor_model_tuned.stan

data {
  int<lower=1> N;        // # observations
  int<lower=1> P;        // # variables
  int<lower=1> K;        // truncation for # factors
  matrix[N, P] Y;        // data matrix (already centered/scaled)
}

parameters {
  matrix[P, K] raw_Lambda;    // non-centered loadings
  matrix[N, K] raw_eta;       // non-centered factor scores
  vector<lower=0>[P] psi;     // residual precisions
  matrix<lower=0>[P, K] phi;  // local–shrinkage scales
  vector<lower=0>[K] delta;   // global–shrinkage multipliers
}

transformed parameters {
  vector<lower=0>[K] tau;
  matrix[P, K]       lambda_sd;
  matrix[P, K]       Lambda;
  matrix[N, K]       eta;

  // build cumulative shrinkage
  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k-1] * delta[k];

  // non-centered transform for loadings
  for (p in 1:P)
    for (k in 1:K) {
      lambda_sd[p,k] = sqrt(1.0 / (phi[p,k] * tau[k]));
      Lambda[p,k]    = raw_Lambda[p,k] * lambda_sd[p,k];
    }

  // non-centered transform for scores
  eta = raw_eta;
}

model {
  // same MGSP priors as in linearMGSP()
  psi              ~ gamma(1.0, 0.3);
  to_vector(phi)   ~ gamma(1.5, 1.5);
  delta[1]         ~ gamma(2.1, 1.0);
  for (k in 2:K)
    delta[k]       ~ gamma(3.1, 1.0);

  // standard normals on raw parameters
  to_vector(raw_Lambda) ~ normal(0, 1);
  to_vector(raw_eta)    ~ normal(0, 1);

  // likelihood
  {
    matrix[N, P] mu = eta * Lambda';
    for (p in 1:P)
      Y[,p] ~ normal(mu[,p], sqrt(1.0 / psi[p]));
  }
}

