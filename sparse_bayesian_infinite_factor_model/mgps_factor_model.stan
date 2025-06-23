data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> K;
  matrix[N, P] Y;
}

parameters {
  matrix[P, K] Lambda;
  matrix[N, K] eta;
  vector<lower=0>[P] psi;
  matrix<lower=0>[P, K] phi;
  vector<lower=0>[K] delta;
}

transformed parameters {
  vector<lower=0>[K] tau;
  matrix[P, K] lambda_sd;

  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k - 1] * delta[k];

  for (p in 1:P)
    for (k in 1:K)
      lambda_sd[p, k] = sqrt(1.0 / (phi[p, k] * tau[k]));
}

model {
  real a_psi = 1.0;
  real b_psi = 0.3;
  real nu    = 3.0;
  real a1    = 2.1;
  real b1    = 1.0;
  real a2    = 3.1;
  real b2    = 1.0;

  psi  ~ gamma(a_psi, b_psi);
  to_vector(phi) ~ gamma(nu / 2, nu / 2);
  delta[1]      ~ gamma(a1, b1);
  for (k in 2:K)
    delta[k]    ~ gamma(a2, b2);

  for (p in 1:P)
    for (k in 1:K)
      Lambda[p, k] ~ normal(0, lambda_sd[p, k]);

  to_vector(eta) ~ normal(0, 1);

  {
    matrix[N, P] mu = eta * Lambda';
    for (p in 1:P)
      Y[, p] ~ normal(mu[, p], sqrt(1.0 / psi[p]));
  }
}

