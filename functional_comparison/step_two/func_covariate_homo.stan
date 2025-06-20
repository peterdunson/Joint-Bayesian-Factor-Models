// Functional MGPS factor model with covariate adjustment (homoscedastic)

// Data:
//  N       = number of subjects
//  n_obs   = total observations across all subjects
//  J       = number of covariates
//  H       = max number of latent factors
//  M       = number of spline basis functions
//  y[n]    = observed value at each (subject,time)
//  subj[n] = subject index for observation n (1…N)
//  Z[N,J]  = covariate matrix (one row per subject)
//  B[n,M]  = spline basis at each observation’s time

data {
  int<lower=1> N;
  int<lower=1> n_obs;
  int<lower=1> J;
  int<lower=1> H;
  int<lower=1> M;
  vector[n_obs]       y;
  int<lower=1,upper=N> subj[n_obs];
  matrix[N, J]        Z;
  matrix[n_obs, M]    B;
}

parameters {
  matrix[J, M]        Beta;      // covariate‐effect basis coeffs
  real<lower=0>       sigma_beta;

  matrix[H, M]        Theta;     // loading‐function basis coeffs
  matrix<lower=0>[H, M] phi;     // local MGPS scales
  vector<lower=0>[H]  delta;     // global MGPS increments

  matrix[N, H]        eta;       // factor scores
  real<lower=0>       sigma_y;   // residual SD
}

transformed parameters {
  vector[H]    tau;       // global shrinkage per factor
  matrix[H, M] ls;        // per‐loading SD
  vector[n_obs] mu;       // predicted mean

  // 1) build MGPS weights τ[h] = ∏_{r=1}^h δ[r]
  tau[1] = delta[1];
  for (h in 2:H)
    tau[h] = tau[h-1] * delta[h];

  // 2) per‐loading SD = 1/sqrt(phi * tau)
  for (h in 1:H)
    for (m in 1:M)
      ls[h,m] = sqrt(1.0 / (phi[h,m] * tau[h]));

  // 3) reconstruct mu[n] = covariate part + factor part
  for (n in 1:n_obs) {
    vector[M] b = B[n]';
    // covariate part: Z[subj]·(Beta * b)
    vector[J] cb = Beta * b;
    real cov_part = dot_product(row(Z, subj[n]), cb);
    // factor part: η[subj]·(Theta * b)
    vector[H] fb = Theta * b;
    real fac_part = dot_product(eta[subj[n]], fb);
    mu[n] = cov_part + fac_part;
  }
}

model {
  // hyper‐parameters
  real nu = 3.0;
  real a1 = 2.1;
  real a2 = 3.1;

  // 1) priors on covariate‐effects
  to_vector(Beta) ~ normal(0, sigma_beta);
  sigma_beta      ~ cauchy(0, 2);

  // 2) MGPS shrinkage priors
  to_vector(phi)   ~ gamma(nu/2, nu/2);
  delta[1]         ~ gamma(a1, 1);
  for (h in 2:H)
    delta[h]       ~ gamma(a2, 1);

  // 3) loading‐basis priors with shrinkage SD
  for (h in 1:H)
    for (m in 1:M)
      Theta[h,m]  ~ normal(0, ls[h,m]);

  // 4) factor scores
  to_vector(eta)   ~ normal(0, 1);

  // 5) residual SD prior
  sigma_y         ~ cauchy(0, 2);

  // 6) likelihood
  y ~ normal(mu, sigma_y);
}

