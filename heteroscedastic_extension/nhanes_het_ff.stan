// File: nhanes_het_ff.stan
// Heteroscedastic Functional Factor Model for NHANES PA Data
// Basis‐function representation with smooth noise‐variance curve

data {
  int<lower=1> N;                // number of subjects
  int<lower=1> n_obs;            // total number of observed (i,t) points
  int<lower=1> H;                // number of latent factors
  int<lower=1> M;                // number of basis functions

  vector[n_obs] y;               // concatenated PA measurements
  int<lower=1,upper=N> subj[n_obs];  // subject index for each observation
  matrix[n_obs, M] B;            // B‐spline basis evaluated at each obs’s time
}

parameters {
  matrix[H, M] Beta;             // factor‐loading coefficients in basis space
  matrix[N, H] eta;              // factor scores per subject
  vector[M] gamma;               // basis coefficients for log‐noise curve

  real<lower=0> sigma_beta;      // prior scale for Beta
  real<lower=0> sigma_gamma;     // prior scale for gamma
}

transformed parameters {
  vector[n_obs] mu;              // mean prediction for each obs
  vector[n_obs] tau;             // noise std‐dev for each obs

  // reconstruct mean and noise at each observation
  for (n in 1:n_obs) {
    // compute factor‐curve values f[h] = sum_m B[n,m] * Beta[h,m]
    vector[H] f;
    for (h in 1:H)
      f[h] = dot_product(row(B, n), Beta[h]);

    // predicted mean = sum_h f[h] * eta[subj[n], h]
    mu[n] = dot_product(f, eta[subj[n]]);

    // log‐noise variance at this time = sum_m B[n,m] * gamma[m]
    tau[n] = exp(0.5 * dot_product(row(B, n), gamma));
  }
}

model {
  // 1) Smoothness priors on loading curves and noise curve
  to_vector(Beta)  ~ normal(0, sigma_beta);
  gamma           ~ normal(0, sigma_gamma);

  // 2) Standard normal priors on factor scores
  to_vector(eta)  ~ normal(0, 1);

  // 3) Hyperpriors on smoothness scales
  sigma_beta  ~ cauchy(0, 2);
  sigma_gamma ~ cauchy(0, 2);

  // 4) Likelihood with time‐varying noise
  y ~ normal(mu, tau);
}

