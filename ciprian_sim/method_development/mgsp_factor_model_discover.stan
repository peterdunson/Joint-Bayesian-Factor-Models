
data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> K;
  matrix[N, P] Y;
  int<lower=1> n_pen;                  // number of penalized pairs
  int<lower=1, upper=P> row_idx[n_pen]; // row indices of penalized pairs
  int<lower=1, upper=P> col_idx[n_pen]; // col indices of penalized pairs
  real<lower=0> pen_weight;            // penalty weight
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
  matrix[N, P] mu;

  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k - 1] * delta[k];
  for (p in 1:P)
    for (k in 1:K)
      lambda_sd[p, k] = sqrt(1.0 / (phi[p, k] * tau[k]));
  mu = eta * Lambda';
}
model {
  // Standard MGSP priors
  psi  ~ gamma(1.0, 0.3);
  to_vector(phi) ~ gamma(1.5, 1.5); // nu = 3
  delta[1]      ~ gamma(2.1, 1.0);
  for (k in 2:K)
    delta[k]    ~ gamma(3.1, 1.0);

  for (p in 1:P)
    for (k in 1:K)
      Lambda[p, k] ~ normal(0, lambda_sd[p, k]);
  to_vector(eta) ~ normal(0, 1);

  // Standard likelihood
  for (p in 1:P)
    Y[, p] ~ normal(mu[, p], sqrt(1.0 / psi[p]));

  // ---- Efficient Penalty on Selected Residual Correlations ----
  {
    matrix[N, P] resid = Y - mu;
    for (n in 1:n_pen) {
      int i = row_idx[n];
      int j = col_idx[n];
      // Compute correlation for this pair
      real num = dot_product(resid[, i], resid[, j]) / (N - 1);
      real si  = sqrt(dot_product(resid[, i], resid[, i]) / (N - 1));
      real sj  = sqrt(dot_product(resid[, j], resid[, j]) / (N - 1));
      real corr_ij = num / (si * sj + 1e-10);
      target += -pen_weight * square(corr_ij);
    }
  }
}


