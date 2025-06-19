// file: hetfofm.stan
data {
  int<lower=1> N;         // number of curves (subjects)
  int<lower=1> p;         // number of time‐points
  int<lower=1> H;         // number of latent factors
  matrix[N,p] X;          // observed curves
}
parameters {
  matrix[H,p]     Lambda;        // loading functions evaluated at each t_k
  matrix[N,H]     eta;           // factor scores for each subject
  vector[p]       log_tau2;      // log noise‐variance at each t_k
  real<lower=0>   sigma_rw;      // RW scale for log_tau2
}
model {
  // 1) Priors on loadings & scores (standard normal here)
  to_vector(Lambda) ~ normal(0,1);
  to_vector(eta)    ~ normal(0,1);

  // 2) Random‐walk prior on log‐variance curve
  sigma_rw ~ cauchy(0,1);            // weak prior for RW‐scale
  for (k in 2:p)
    log_tau2[k] ~ normal(log_tau2[k-1], sigma_rw);
  // you could also pin down log_tau2[1] ~ normal(0,1)

  // 3) Likelihood with time‐varying noise
  for (i in 1:N) {
    for (k in 1:p) {
      real mu = dot_product(Lambda[,k], eta[i]);
      X[i,k] ~ normal(mu, exp(0.5 * log_tau2[k]));
    }
  }
}

