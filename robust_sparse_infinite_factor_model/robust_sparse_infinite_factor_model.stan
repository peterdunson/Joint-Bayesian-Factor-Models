// robust_sparse_infinite_factor_model.stan
// Robust sparse Bayesian infinite factor model (t errors, MGP prior)

data {
  int<lower=1> N;            // Number of observations
  int<lower=1> P;            // Number of variables
  int<lower=1> K;            // Number of factors (truncation)
  matrix[N, P] Y;            // Centered and scaled data
  real<lower=2> nu;          // Degrees of freedom for t-distribution
}

parameters {
  matrix[P, K] Lambda;                // Factor loadings
  matrix[N, K] eta;                   // Factor scores
  matrix<lower=0>[P, K] phi;          // Local shrinkage for Lambda
  vector<lower=0>[K] delta;           // Global MGP shrinkage
  vector<lower=0>[P] sigma2;          // Unique variances
  vector<lower=0>[N] gamma;           // t-mixture scale parameters
}

transformed parameters {
  vector<lower=0>[K] tau;
  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k-1] * delta[k];
}

model {
  // MGP global-local shrinkage
  delta[1] ~ gamma(2.1, 1.0);
  for (k in 2:K)
    delta[k] ~ gamma(3.1, 1.0);
  
  for (j in 1:P) {
    sigma2[j] ~ inv_gamma(1.0, 0.3);
    for (h in 1:K) {
      phi[j, h] ~ gamma(3.0/2.0, 3.0/2.0); // for nu=3
      Lambda[j, h] ~ normal(0, sqrt(1.0 / (phi[j, h] * tau[h])));
    }
  }
  for (i in 1:N) {
    gamma[i] ~ gamma(nu / 2, nu / 2);
    eta[i] ~ normal(0, 1);
    Y[i]' ~ multi_normal(Lambda * eta[i]', diag_matrix(sigma2) / gamma[i]);
  }
}





