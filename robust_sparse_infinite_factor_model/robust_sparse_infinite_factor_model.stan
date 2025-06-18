// robust_sparse_infinite_factor_model.stan
// Robust sparse Bayesian infinite factor model (Lee, Jo, Lee, 2022, Comp Stat)
// t-distributed errors via scale mixture, MGP prior for factor loadings

data {
  int<lower=1> N;            // Number of samples
  int<lower=1> P;            // Number of observed variables
  int<lower=1> K;            // Truncation level (should be large, e.g. 30)
  matrix[N, P] Y;            // Data (rows: obs, cols: variables)
  real<lower=2> nu;          // Degrees of freedom for t-distribution (suggest 3 or 5)
}

parameters {
  matrix[P, K] Lambda;                // Factor loadings
  matrix[N, K] eta;                   // Latent factors (each row i is eta[i,])
  matrix<lower=0>[P, K] phi;          // Local shrinkage
  vector<lower=0>[K] delta;           // MGP global shrinkage
  vector<lower=0>[P] sigma2;          // Unique variances (diagonal of Sigma)
  vector<lower=0>[N] gamma;           // Scale mixture auxiliary (for Student-t)
}

transformed parameters {
  vector<lower=0>[K] tau;
  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k-1] * delta[k];
}

model {
  // -- Shrinkage priors (as in the paper, Sec 3.2)
  delta[1] ~ gamma(2.1, 1.0);                 // delta_1 ~ Gamma(a1, b1)
  for (k in 2:K)
    delta[k] ~ gamma(3.1, 1.0);               // delta_{k>1} ~ Gamma(a2, b2)
  
  for (j in 1:P) {
    sigma2[j] ~ inv_gamma(1.0, 0.3);          // Unique variances
    for (h in 1:K) {
      phi[j, h] ~ gamma(3.0 / 2.0, 3.0 / 2.0); // phi_{jh} ~ Gamma(nu/2, nu/2), nu=3
      Lambda[j, h] ~ normal(0, sqrt(1.0 / (phi[j, h] * tau[h])));
    }
  }
  // Student-t auxiliary and latent factors
  for (i in 1:N) {
    gamma[i] ~ gamma(nu / 2, nu / 2); // scale mixture for t-errors
    eta[i] ~ normal(0, 1);
    Y[i]' ~ multi_normal(Lambda * eta[i]', diag_matrix(sigma2) / gamma[i]);
    // eta[i]' is a column vector, Lambda * eta[i]' is (P x 1)
    // Y[i]' is a row vector, Stan will match dims.
  }
}

generated quantities {
  // Optionally add posterior predictive quantities, e.g. y_rep, etc.
}
