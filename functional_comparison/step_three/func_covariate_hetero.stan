// Model 3: covariate‐adjusted residuals + heteroskedastic errors
// Functional MGPS factor model with time‐varying noise

data {
  int<lower=1> N;                // number of subjects
  int<lower=1> n_obs;            // total number of observations
  int<lower=1> J;                // number of covariates
  int<lower=1> H;                // max number of latent factors
  int<lower=1> M;                // number of spline basis functions

  vector[n_obs]       y;         // observed values
  int<lower=1,upper=N> subj[n_obs]; // subject index for each obs
  matrix[N, J]        Z;         // covariate matrix (subjects × covariates)
  matrix[n_obs, M]    B;         // spline basis at each observation’s time
}

parameters {
  // 1) covariate effects
  matrix[J, M]        Beta;        // basis coeffs for each covariate‐effect curve
  real<lower=0>       sigma_beta;  // scale for Beta

  // 2) factor loadings + MGPS
  matrix[H, M]        Theta;       // basis coeffs for each loading‐function
  matrix<lower=0>[H, M] phi;       // local MGPS scales
  vector<lower=0>[H]  delta;       // global MGPS increments

  // 3) factor scores
  matrix[N, H]        eta;         // subject scores

  // 4) heteroskedastic noise
  vector[M]           gamma;       // basis coeffs for log‐variance curve
  real<lower=0>       sigma_gamma; // scale for gamma
}

transformed parameters {
  vector[H]    tau;        // cumulative MGPS weights
  matrix[H, M] ls;         // per‐loading SDs
  vector[n_obs] mu;        // predicted mean
  vector[n_obs] tau_obs;   // predicted noise SD

  // build MGPS weights τ[h] = ∏_{r=1}^h δ[r]
  tau[1] = delta[1];
  for (h in 2:H)
    tau[h] = tau[h-1] * delta[h];

  // loading SD = 1 / sqrt(phi * tau)
  for (h in 1:H)
    for (m in 1:M)
      ls[h,m] = sqrt(1.0 / (phi[h,m] * tau[h]));

  // reconstruct mu[n] and noise SD tau_obs[n]
  for (n in 1:n_obs) {
    vector[M] b = B[n]';
    // covariate part
    real cov_part = dot_product(row(Z, subj[n]), Beta * b);
    // factor part
    vector[H] fb; 
    for (h in 1:H)
      fb[h] = dot_product(Theta[h], b);
    real fac_part = dot_product(eta[subj[n]], fb);
    mu[n] = cov_part + fac_part;
    // heteroskedastic SD
    tau_obs[n] = exp(0.5 * dot_product(b, gamma));
  }
}

model {
  // hyperparameters for MGPS
  real nu = 3.0;
  real a1 = 2.1;
  real a2 = 3.1;

  // 1) priors on covariate coefficients
  to_vector(Beta)  ~ normal(0, sigma_beta);
  sigma_beta       ~ cauchy(0, 2);

  // 2) MGPS shrinkage priors for loadings
  to_vector(phi)   ~ gamma(nu/2, nu/2);
  delta[1]         ~ gamma(a1, 1);
  for (h in 2:H)
    delta[h]       ~ gamma(a2, 1);

  // 3) loading‐basis priors
  for (h in 1:H)
    for (m in 1:M)
      Theta[h,m]  ~ normal(0, ls[h,m]);

  // 4) factor‐score priors
  to_vector(eta)   ~ normal(0, 1);

  // 5) noise‐curve priors
  to_vector(gamma) ~ normal(0, sigma_gamma);
  sigma_gamma      ~ cauchy(0, 2);

  // 6) likelihood with heteroskedastic errors
  y ~ normal(mu, tau_obs);
}

