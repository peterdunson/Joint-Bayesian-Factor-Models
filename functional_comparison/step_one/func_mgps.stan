// File: func_mgps.stan
// Direct functional factor model with MGPS shrinkage on loading‐function basis

data {
  int<lower=1> N;                
  int<lower=1> n_obs;            // total number of observations across all subjects
  int<lower=1> H;                // truncation level
  int<lower=1> M;                // number of spline basis functions

  vector[n_obs]       y;         // stacked observations y[i,k]
  int<lower=1,upper=N> subj[n_obs]; // subject index for each observation
  matrix[n_obs, M]    B;         // spline basis evaluated at each obs time
}

parameters {
  matrix[H, M]        Theta;     // basis coefficients for each loading‐function
  matrix[N, H]        eta;       // factor scores per subject
  matrix<lower=0>[H, M] phi;     // local shrinkage scales (MGPS)
  vector<lower=0>[H]  delta;     // global shrinkage increments
  real<lower=0>       sigma_y;   // homoscedastic residual SD
}

transformed parameters {
  vector[H]    tau;       // cumulative MGPS weights
  matrix[H, M] ls;        // per‐loading SD = 1/sqrt(phi * tau)
  vector[n_obs] mu;       // predicted mean for each obs

  // 1) build global shrinkage weights τ[h] = ∏_{r=1}^h δ[r]
  tau[1] = delta[1];
  for (h in 2:H)
    tau[h] = tau[h-1] * delta[h];

  // 2) compute per‐loading SDs
  for (h in 1:H)
    for (m in 1:M)
      ls[h,m] = sqrt(1.0 / (phi[h,m] * tau[h]));

  // 3) reconstruct mu[n] = sum_{h=1}^H η[subj[n],h] * (Θ[h,] · B[n,])
  for (n in 1:n_obs) {
    vector[M] b = B[n]';
    vector[H] fb;
    for (h in 1:H)
      fb[h] = dot_product(Theta[h], b);
    mu[n] = dot_product(eta[subj[n]], fb);
  }
}

model {
  // hyper‐parameters (paper defaults)
  real nu = 3.0;
  real a1 = 2.1;
  real a2 = 3.1;

  // 1) MGPS shrinkage priors
  to_vector(phi)   ~ gamma(nu/2, nu/2);
  delta[1]         ~ gamma(a1, 1);
  for (h in 2:H)
    delta[h]       ~ gamma(a2, 1);

  // 2) loading‐basis priors with MGPS SD
  for (h in 1:H)
    for (m in 1:M)
      Theta[h,m]  ~ normal(0, ls[h,m]);

  // 3) factor‐score priors
  to_vector(eta)   ~ normal(0, 1);

  // 4) residual SD prior
  sigma_y         ~ cauchy(0, 2);

  // 5) likelihood
  y ~ normal(mu, sigma_y);
}

