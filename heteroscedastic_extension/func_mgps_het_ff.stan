// File: func_mgps_het_ff.stan
// Fully Bayesian functional factor model with:
//  1) covariate effects in the mean,
//  2) heteroskedastic error variance τ²(t),
//  3) MGPS shrinkage on the latent‐factor loadings.

data {
  int<lower=1> N;                // number of subjects
  int<lower=1> n_obs;            // total number of observations
  int<lower=1> J;                // number of covariates
  int<lower=1> H;                // number of latent factors
  int<lower=1> M;                // number of spline basis functions

  vector[n_obs]       y;         // observed values (stacked over subjects & times)
  int<lower=1,upper=N> subj[n_obs];   // subject index for each observation
  matrix[N, J]        Z;         // covariate matrix (one row per subject)
  matrix[n_obs, M]    B;         // spline basis evaluated at each obs’s time
}

parameters {
  matrix[J, M]        Beta;      // basis coeffs for covariate functions
  matrix[H, M]        Theta;     // basis coeffs for loading functions
  vector[M]           gamma;     // basis coeffs for log‐noise curve

  real<lower=0>       sigma_beta;   // overall scale for Beta
  real<lower=0>       sigma_gamma;  // overall scale for gamma

  matrix[N, H]        eta;       // subject factor scores

  matrix<lower=0>[H, M] phi;      // local shrinkage scales for loadings
  vector<lower=0>[H]  delta;     // global shrinkage increments
}

transformed parameters {
  vector[H]           tau_mgps;    // cumulative MGPS weights
  matrix[H, M]        lambda_sd;   // per‐loading SDs
  vector[n_obs]       mu;          // predicted means
  vector[n_obs]       tau_obs;     // predicted noise SDs

  // 1) build MGPS weights τₕ = ∏_{r=1}^h δ_r
  tau_mgps[1] = delta[1];
  for (h in 2:H)
    tau_mgps[h] = tau_mgps[h - 1] * delta[h];

  // 2) per‐loading SD = 1 / sqrt(φ[h,m] * τₕ)
  for (h in 1:H)
    for (m in 1:M)
      lambda_sd[h,m] = sqrt(1.0 / (phi[h,m] * tau_mgps[h]));

  // 3) reconstruct means and noise‐SDs for each observation
  for (n in 1:n_obs) {
    vector[M] b = B[n]';                 // basis at this time
    // covariate part: Z[subj]⋅(Beta * b)
    vector[J]  cb = Beta * b;            
    real       cov_part = dot_product(row(Z, subj[n]), cb);
    // factor part:    η[subj]⋅(Theta * b)
    vector[H]  fb = Theta * b;
    real       fac_part = dot_product(eta[subj[n]], fb);
    mu[n]      = cov_part + fac_part;

    // noise‐SD: exp(½⋅b⋅γ)
    tau_obs[n] = exp(0.5 * dot_product(b, gamma));
  }
}

model {
  // hyperparameters for MGPS
  real nu = 3.0;
  real a1 = 2.1;
  real a2 = 3.1;

  // 1) smoothness priors on covariate & noise‐basis coeffs
  to_vector(Beta)  ~ normal(0, sigma_beta);
  gamma            ~ normal(0, sigma_gamma);
  sigma_beta       ~ cauchy(0, 2);
  sigma_gamma      ~ cauchy(0, 2);

  // 2) MGPS shrinkage priors for loadings
  to_vector(phi)   ~ gamma(nu/2, nu/2);
  delta[1]         ~ gamma(a1, 1);
  for (h in 2:H)
    delta[h]       ~ gamma(a2, 1);

  // 3) loading‐basis priors with shrinkage SDs
  for (h in 1:H)
    for (m in 1:M)
      Theta[h,m]  ~ normal(0, lambda_sd[h,m]);

  // 4) standard normal priors on factor scores
  to_vector(eta)   ~ normal(0, 1);

  // 5) heteroskedastic normal likelihood
  y ~ normal(mu, tau_obs);
}


