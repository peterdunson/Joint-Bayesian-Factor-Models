// =======================================================================
// NOTE: This Stan model implements the "sparse" hierarchical prior
// structure from BASS (Zhao et al., JMLR 2016), using a three-level
// global-factor-local shrinkage prior (i.e., structured horseshoe/TBP).
// -----------------------------------------------------------------------
// The FULL BASS model uses a *mixture* prior:
//    - Each factor is either "sparse", "dense", or "null", controlled
//      by a discrete indicator z_k (sampled in Gibbs/EM, not Stan).
//    - "Dense" factors are regularized only at the column (factor) level.
//    - "Sparse" factors use element-wise shrinkage (as below).
//    - "Null" factors are dropped (all zeros).
// This Stan implementation treats ALL factors as "sparse", i.e. with
// element-wise adaptive shrinkage. Mixture indicator variables are NOT
// included, since Stan cannot sample discrete indicators.
// =======================================================================

data {
  int<lower=1> N;   // number of samples
  int<lower=1> P;   // number of features
  int<lower=1> K;   // number of latent factors
  matrix[N, P] Y;   // centered data matrix
}

parameters {
  // Latent factors
  matrix[N, K] X;

  // Loadings matrix and shrinkage parameters
  matrix[P, K] Lambda;

  // Three-level global-factor-local shrinkage (structured horseshoe/TBP)
  vector<lower=0>[K] tau_global;      // global shrinkage (one per factor)
  vector<lower=0>[K] tau_factor;      // factor-specific shrinkage
  matrix<lower=0>[P, K] tau_local;    // element-specific shrinkage

  // Residual variances
  vector<lower=0>[P] sigma2;
}

transformed parameters {
  // Standard deviations for loadings prior
  matrix<lower=0>[P, K] lambda_sd;
  for (k in 1:K)
    for (j in 1:P)
      lambda_sd[j, k] = sqrt(tau_global[k] * tau_factor[k] * tau_local[j, k]);
}

model {
  // Shrinkage hyperpriors (per BASS paper, Section 4.1/Table 1: a = b = c = d = e = f = 0.5)
  tau_global ~ gamma(0.5, 0.5);                  // global
  tau_factor ~ gamma(0.5, tau_global);           // factor-specific
  for (k in 1:K)
    for (j in 1:P)
      tau_local[j, k] ~ gamma(0.5, tau_factor[k]); // element-wise

  // Residual variances (same as other models for comparability)
  sigma2 ~ inv_gamma(1, 0.3);

  // Loadings
  for (k in 1:K)
    for (j in 1:P)
      Lambda[j, k] ~ normal(0, lambda_sd[j, k]);

  // Latent factors
  to_vector(X) ~ normal(0, 1);

  // Data likelihood
  for (i in 1:N)
    Y[i, ] ~ normal(Lambda * X[i, ]', sqrt(sigma2)');
}


