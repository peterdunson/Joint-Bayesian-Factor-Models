// LowFR.stan
// Finalized model from "Low-rank longitudinal factor regression" paper
// This implementation does not include additional covariates and cannot handle missing data.

functions{
  matrix kronecker(matrix A, matrix B) {
    matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
    int m;
    int n;
    int p;
    int q;
    m = rows(A);
    n = cols(A);
    p = rows(B);
    q = cols(B);
    for (i in 1:m) {
      for (j in 1:n) {
        int row_start;
        int row_end;
        int col_start;
        int col_end;
        row_start = (i - 1) * p + 1;
        row_end = (i - 1) * p + p;
        col_start = (j - 1) * q + 1;
        col_end = (j - 1) * q + q;
        C[row_start:row_end, col_start:col_end] = A[i, j] * B;
      }
    }
    return C;
  }
}

data {
  int<lower=0> N;         // number of subjects
  int<lower=0> p;         // number of exposures
  int<lower=0> k;         // number of latent factors
  int<lower=0> TT;        // number of time points
  int<lower=0> H;         // rank hyperparameter (default: min(p, TT))
  vector[N] y;            // outcome
  matrix[N, p*TT] X;      // longitudinal exposures matrix (vectorized per subject)
}

parameters {
  real mu;                    // regression intercept
  matrix[k,H] beta;           // main effects: factors × rank
  matrix[TT,H] omega;         // main effects: time × rank
  matrix[k,k] B;              // quadratic (interaction) factors
  matrix[TT,TT] W;            // quadratic (interaction) time
  real<lower=0> sigma2;       // outcome variance

  // Factor model
  matrix[p,k] Lambda;
  vector<lower=0>[p] Sigma;
  real<lower=0, upper=1> phi;
  matrix[N, k*TT] Eta;        // subject-specific factor trajectories

  // Multiplicative gamma process prior (main effects)
  real<lower=0> delta[H];
  real<lower=0> xi[k+TT,H];
  real<lower=0> a1;
  real<lower=0> a2;

  // Multiplicative gamma process prior (interactions)
  real<lower=0> tau_int;
  real<lower=0> xi_int[k*k+TT*TT];
  real<lower=0> a1_int;
}

transformed parameters {
  vector[k*TT] theta;
  matrix[k*TT, k*TT] Omega;
  matrix[TT, TT] Phi_mat;
  real<lower=0> tau[H];

  // Construct theta (main effect: time × factors)
  theta = to_vector(omega * beta');

  // Construct Omega (interaction effect: Kronecker of factor and time)
  Omega = kronecker(B, W);
  Omega = (Omega + Omega') / 2;

  // Compound symmetric temporal correlation
  Phi_mat = rep_matrix(0, TT, TT);
  for (i in 1:TT) {
    for (j in 1:TT) {
      Phi_mat[i, j] = (i == j) ? 1 : phi;
    }
  }

  // Multiplicative gamma process
  tau[1] = delta[1];
  for (l in 2:H)
    tau[l] = tau[l-1] * delta[l];
}

model {
  matrix[k*TT, k*TT] I_kron_phi = kronecker(diag_matrix(rep_vector(1, k)), Phi_mat);
  matrix[p*TT, p*TT] Sigma_kron_phi = kronecker(diag_matrix(Sigma), Phi_mat);
  matrix[TT*p, TT*k] Lambda_kron_I = kronecker(Lambda, diag_matrix(rep_vector(1, TT)));

  // Priors for variance terms
  Sigma ~ inv_gamma(1, 1);
  sigma2 ~ inv_gamma(1, 1);

  // Intercept
  mu ~ normal(0, sqrt(10));

  // Multiplicative gamma process priors: main effects
  for (j in 1:k)
    for (l in 1:H)
      beta[j, l] ~ normal(0, 1 / sqrt(xi[j, l] * tau[l]));
  for (t in 1:TT)
    for (l in 1:H)
      omega[t, l] ~ normal(0, 1 / sqrt(xi[k+t, l] * tau[l]));
  delta[1] ~ gamma(a1, 1);
  for (l in 2:H)
    delta[l] ~ gamma(a2, 1);
  for (j in 1:(k+TT))
    for (l in 1:H)
      xi[j, l] ~ gamma(1.5, 1.5);
  a1 ~ gamma(2, 1);
  a2 ~ gamma(2, 1);

  // Interactions: MGP
  tau_int ~ gamma(a1_int, 1);
  xi_int ~ gamma(1.5, 1.5);
  a1_int ~ gamma(2, 1);

  // Quadratic regression terms (interactions)
  for (i in 1:k)
    for (j in 1:k)
      B[i, j] ~ normal(0, 1 / sqrt(xi_int[(i-1)*k + j] * tau_int));
  for (i in 1:TT)
    for (j in 1:TT)
      W[i, j] ~ normal(0, 1 / sqrt(xi_int[k*k + (i-1)*TT + j] * tau_int));

  // Factor model priors
  for (i in 1:N)
    Eta[i, ] ~ multi_normal(rep_vector(0, k*TT), I_kron_phi);
  for (i in 1:p)
    for (j in 1:k)
      Lambda[i, j] ~ normal(0, sqrt(10));
  phi ~ uniform(0, 1);

  // Likelihood: exposures and outcome
  for (i in 1:N) {
    X[i, ] ~ multi_normal(Lambda_kron_I * to_vector(Eta[i, ]), Sigma_kron_phi);
    y[i] ~ normal(mu + Eta[i] * theta + quad_form(Omega, Eta[i, ]'), sqrt(sigma2));
  }
}

generated quantities {
  // Induced regression effects of y on observed X
  matrix[k*TT, p*TT] A;
  matrix[k*TT, k*TT] V;
  matrix[p, p] Sigma_inv;
  real alpha_0;
  vector[p*TT] alpha;
  matrix[p*TT, p*TT] Gamma;

  Sigma_inv = rep_matrix(0, p, p);
  for (j in 1:p)
    Sigma_inv[j, j] = 1 / Sigma[j];
  V = kronecker(inverse_spd(Lambda' * Sigma_inv * Lambda + diag_matrix(rep_vector(1, k))), Phi_mat);
  A = V * kronecker(Lambda' * Sigma_inv, inverse_spd(Phi_mat));
  alpha_0 = mu + trace(Omega * V);
  alpha = A' * theta;
  Gamma = A' * Omega * A;
}
