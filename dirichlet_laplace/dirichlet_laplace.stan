data {
  int<lower=1> N;        // # observations
  int<lower=1> P;        // # variables
  int<lower=1> K;        // # factors
  matrix[N, P] Y;        // Data
}

parameters {
  matrix[P, K] Lambda;           // Factor loadings
  matrix[N, K] eta;              // Factor scores
  vector<lower=0>[P] psi;        // Specific variances

  simplex[P] phi[K];             // Dirichlet scales for each factor (K x P simplex)
  vector<lower=0>[K] tau;        // Global scale for each factor
}

model {
  // Hyperparameters for Dirichlet-Laplace
  real a_psi = 1.0;
  real b_psi = 0.3;
  real a_tau = 1.0;       // You can tune a_tau for global shrinkage

  // Priors
  psi ~ gamma(a_psi, b_psi);

  for (k in 1:K) {
    tau[k] ~ gamma(a_tau, a_tau);    // Global scale, can be half-Cauchy for more shrinkage
    phi[k] ~ dirichlet(rep_vector(1.0, P));  // Uniform Dirichlet, can use alpha < 1 for more sparsity

    for (p in 1:P) {
      // Laplace = scale mixture of normals
      Lambda[p, k] ~ double_exponential(0, phi[k][p] * tau[k]);
      // If Stan version <2.29, use normal_lpdf + exponential prior for scale mixture
    }
  }
  to_vector(eta) ~ normal(0, 1);

  // Likelihood
  for (n in 1:N) {
    Y[n] ~ normal(eta[n] * Lambda', sqrt(inv(psi)));
  }
}
