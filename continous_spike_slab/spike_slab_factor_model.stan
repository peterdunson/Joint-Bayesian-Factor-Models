// spike_slab_factor_model.stan
// -----------------------------------------------------------
// Implements Bayesian Factor Analysis using a continuous 
// spike-and-slab prior on the factor loadings (λ_jk) as described in:
//   Lu, Zhao-Hua, Sy-Miin Chow, & Eric Loken (2016)
//   "Bayesian Factor Analysis as a Variable-Selection Problem: 
//    Alternative Priors and Consequences"
//   Multivariate Behavioral Research, 51:4, 519-539.
//
// - Each factor loading λ_jk (except main loadings, see below) is given a
//   continuous spike-and-slab prior: 
//       λ_jk ~ θ * Normal(0, sigma_slab) + (1-θ) * Normal(0, sigma_spike)
//   where sigma_spike is a small value (narrow, strong shrinkage to zero),
//   sigma_slab is large (weak shrinkage), and θ is the prior inclusion probability.
//
// - The "spike" is not a true point mass at zero (as in discrete spike-and-slab) 
//   due to Stan's limitations. Instead, we use a very small variance Normal as the spike.
// - The indicator variable (r_jk in the paper) is marginalized out using log_mix().
// - Main loadings (primary, theory-supported) can be excluded from spike-and-slab and 
//   given a slab prior only. 
// - This code is for structure learning/comparison—see Lu et al. (2016) for details.
//
// -----------------------------------------------------------

data {
  int<lower=1> N;           // Number of observations
  int<lower=1> P;           // Number of variables (first column is outcome)
  int<lower=1> K;           // Number of factors
  matrix[N, P] Y;           // Centered data matrix
  real<lower=0> Sigma1;     // Fixed variance for outcome
  int main_loadings[P];     // For each variable: which loading is the "main" (0 if none)
}
parameters {
  matrix[P, K] Lambda;             // Factor loadings
  matrix[N, K] eta;                // Factor scores
  vector<lower=0>[P-1] psi_free;   // Residual precisions
  real<lower=0> sigma_slab;        // Slab std (large)
  real<lower=0> sigma_spike;       // Spike std (small)
  real<lower=0, upper=1> theta;    // Prior inclusion probability
}
transformed parameters {
  vector<lower=0>[P] psi;
  psi[1] = 1.0 / Sigma1;
  for (j in 2:P)
    psi[j] = psi_free[j-1];
}
model {
  // Priors for slab/spike
  sigma_slab ~ normal(0, 1.0);    // weakly informative
  sigma_spike ~ normal(0, 0.5);   // weakly informative, center at small value
  theta ~ beta(1, 1);             // uniform prior

  psi_free ~ gamma(1.0, 0.3);
  to_vector(eta) ~ normal(0, 1);

  // Spike-and-slab prior for each loading
  for (j in 1:P) {
    for (k in 1:K) {
      if (main_loadings[j] == k) {
        // Main loading: slab only (unshrunk)
        Lambda[j, k] ~ normal(0, sigma_slab);
      } else {
        // All other loadings: continuous spike-and-slab prior
        target += log_mix(theta,
          normal_lpdf(Lambda[j, k] | 0, sigma_slab),
          normal_lpdf(Lambda[j, k] | 0, sigma_spike));
      }
    }
  }

  // Likelihood
  {
    matrix[N, P] mu = eta * Lambda';
    for (p in 1:P)
      Y[, p] ~ normal(mu[, p], sqrt(1.0 / psi[p]));
  }
}

