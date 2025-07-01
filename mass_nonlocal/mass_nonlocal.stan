// spike_pmom_factor_model.stan
// Sparse Bayesian factor model with MGSP prior on loadings and
// spike-and-pMOM prior on factor scores (continuous spike, Stan compatible)

data {
    int<lower=1> N;            // Number of samples
    int<lower=1> P;            // Number of observed variables
    int<lower=1> K;            // Number of latent factors (overestimate, allow shrinkage)
    matrix[N, P] Y;            // Data matrix (centered, optionally standardized)
    
    // Hyperparameters
    real<lower=0> a_theta;     // Beta prior shape for theta
    real<lower=0> b_theta;
    real<lower=0> a_sigma;     // Gamma prior shape for sigma^2
    real<lower=0> b_sigma;
    real<lower=0> nu;          // Degrees of freedom for local shrinkage (MGPS)
    real<lower=0> a1;          // MGPS prior
    real<lower=0> a2;          // MGPS prior
    real<lower=0> psi;         // Dispersion parameter for pMOM
    real<lower=0> spike_sd;    // Standard deviation for the "continuous spike"
}

parameters {
    // Factor loadings: P x K
    matrix[P, K] Lambda;

    // Factor scores: N x K
    matrix[N, K] eta;

    // Local and global shrinkage for MGPS
    matrix<lower=0>[P, K] phi;     // local
    vector<lower=0>[K] delta;      // global, cumulative product

    // Variance
    vector<lower=0>[P] sigma2;

    // Mixture weight for factor scores
    vector<lower=0, upper=1>[K] theta;  // Bernoulli prob for each factor
}

transformed parameters {
    vector<lower=0>[K] tau;
    tau[1] = delta[1];
    for (h in 2:K) {
        tau[h] = tau[h-1] * delta[h];
    }
}

model {
    // --- Priors ---
    // MGPS prior for Lambda
    for (j in 1:P) {
        for (h in 1:K) {
            Lambda[j, h] ~ normal(0, sqrt(1 / (phi[j, h] * tau[h])));
            phi[j, h] ~ gamma(nu/2, nu/2);
        }
    }
    delta[1] ~ gamma(a1, 1);
    for (h in 2:K) delta[h] ~ gamma(a2, 1);

    // Residual variances
    sigma2 ~ gamma(a_sigma, b_sigma);

    // Theta (sparsity)
    theta ~ beta(a_theta, b_theta);

    // --- Factor scores: mixture prior (continuous spike + pMOM slab) ---
    for (h in 1:K) {
        for (i in 1:N) {
            target += log_mix(theta[h],
                // pMOM slab: eta_ih ~ pMOM(psi)
                // Density: (1/sqrt(2 * pi * psi^3)) * exp(-eta^2/(2*psi)) * eta^2
                log(1 / sqrt(2 * pi() * psi^3)) - 0.5 * square(eta[i, h]) / psi + 2 * log(fabs(eta[i, h])),
                // Continuous spike at zero:
                normal_lpdf(eta[i, h] | 0, spike_sd)
            );
        }
    }

    // --- Likelihood ---
    for (i in 1:N) {
        Y[i] ~ multi_normal(Lambda * eta[i]', diag_matrix(sigma2));
    }
}



