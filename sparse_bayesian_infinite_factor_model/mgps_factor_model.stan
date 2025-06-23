

// mgps_factor_model.stan
// Sparse Bayesian infinite factor model (Bhattacharya & Dunson, 2011)
// + non-centered Λ, log-ψ & log-φ with Jacobian adjustments

data {
  int<lower=1> N;           // number of observations
  int<lower=1> P;           // number of variables
  int<lower=1> K;           // truncation level
  matrix[N, P] Y;           // centered & scaled data matrix
}

parameters {
  matrix[P, K] Lambda_raw;     // non-centered loadings
  matrix[N, K] eta;            // latent factor scores
  vector[P]      log_psi;      // unconstrained log residual precisions
  matrix[P, K]   log_phi;      // unconstrained log local shrinkage scales
  vector<lower=0>[K] delta;    // global shrinkage weights
}

transformed parameters {
  vector<lower=0>[P] psi;         // residual precisions
  matrix<lower=0>[P, K] phi;      // local shrinkage scales
  vector<lower=0>[K] tau;         
  matrix[P, K]      lambda_sd;    
  matrix[P, K]      Lambda;       

  // 1) build psi = exp(log_psi)
  for (p in 1:P)
    psi[p] = exp(log_psi[p]);

  // 2) build phi = exp(log_phi)
  for (p in 1:P)
    for (k in 1:K)
      phi[p,k] = exp(log_phi[p,k]);

  // 3) MGPS global weights
  tau[1] = delta[1];
  for (k in 2:K)
    tau[k] = tau[k-1] * delta[k];

  // 4) per-loading SD and non-centered scaling
  for (k in 1:K) {
    lambda_sd[,k] = sqrt(1.0 ./ (phi[,k] * tau[k]));
    Lambda[,k]    = Lambda_raw[,k] .* lambda_sd[,k];
  }
}

model {
  // hyperparameters (paper defaults)
  real a_psi = 1.0, b_psi = 0.3;
  real nu    = 3.0;
  real a1    = 2.1, b1 = 1.0;
  real a2    = 3.1, b2 = 1.0;

  // 1) log-ψ prior + Jacobian for psi ~ Gamma(a_psi, b_psi)
  for (p in 1:P) {
    target += gamma_lpdf(psi[p] | a_psi, b_psi)
            + log_psi[p];  // Jacobian: dψ/dlog_ψ = ψ
  }

  // 2) log-φ prior + Jacobian for phi ~ Gamma(nu/2, nu/2)
  for (p in 1:P)
    for (k in 1:K) {
      target += gamma_lpdf(phi[p,k] | nu/2, nu/2)
              + log_phi[p,k];  // Jacobian: dφ/dlog_φ = φ
    }

  // 3) global shrinkage δ ~ Gamma(a1,b1) & Gamma(a2,b2)
  delta[1] ~ gamma(a1, b1);
  for (k in 2:K)
    delta[k] ~ gamma(a2, b2);

  // 4) non-centered loading & factor priors
  to_vector(Lambda_raw) ~ normal(0, 1);
  to_vector(eta)        ~ normal(0, 1);

  // 5) likelihood: loop over variables
  {
    matrix[N, P] mu = eta * Lambda';
    for (p in 1:P) {
      Y[,p] ~ normal(mu[,p], sqrt(1.0 / psi[p]));
    }
  }
}

