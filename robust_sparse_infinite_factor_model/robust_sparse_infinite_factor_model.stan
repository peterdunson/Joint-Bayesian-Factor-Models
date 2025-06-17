data {
  int<lower=1> N;            // Number of samples
  int<lower=1> P;            // Number of observed variables
  int<lower=1> K;            // Truncation level (max number of factors)
  matrix[N, P] Y;            // Data (rows: obs, cols: variables)
  real<lower=2> nu;          // Degrees of freedom for t-distribution (fixed, e.g. 3 or 5)
  real<lower=0> a1;          // Hyperparam for delta_1 (usually 2)
  real<lower=0> a2;          // Hyperparam for delta_{l>=2} (usually >3)
  real<lower=0> kappa;       // Hyperparam for local shrinkage (usually 2)
  real<lower=0> a_sigma;     // Hyperparam for error variance
  real<lower=0> b_sigma;     // Hyperparam for error variance
}
transformed data {
  vector[K] ones_K = rep_vector(1.0, K);
}
parameters {
  matrix[P, K] Lambda;             // Factor loadings
  matrix[N, K] eta;                // Latent factors
  matrix<lower=0>[P, K] phi;       // Local shrinkage
  vector<lower=0>[K] delta;        // MGP multiplicative gamma process
  vector<lower=0>[K] tau;          // Cumulative product for global shrinkage
  vector<lower=0>[P] sigma2;       // Unique variances (diagonal of Sigma)
  vector<lower=0>[N] gamma;        // Auxiliary for Student-t
}
model {
  // Shrinkage priors
  delta[1] ~ gamma(a1, 1);
  for (k in 2:K)
    delta[k] ~ gamma(a2, 1);
  tau = cumulative_product(delta);
  
  for (j in 1:P) {
    sigma2[j] ~ inv_gamma(a_sigma, b_sigma);
    for (h in 1:K) {
      phi[j, h] ~ gamma(kappa / 2, kappa / 2);
      Lambda[j, h] ~ normal(0, sqrt(1.0 / (phi[j, h] * tau[h])));
    }
  }
  // Student-t auxiliary variables and latent factors
  for (i in 1:N) {
    gamma[i] ~ gamma(nu / 2, nu / 2);
    eta[i] ~ normal(0, 1);
    Y[i] ~ multi_normal(
      Lambda * eta[i], 
      diag_matrix(sigma2) / gamma[i]
    );
  }
}
generated quantities {
  // Posterior covariance estimate for each sample (optional)
  // (Not required for sampling, but useful for diagnostics)
}
