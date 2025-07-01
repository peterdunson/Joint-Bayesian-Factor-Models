# fit_and_save_spike_pmom_NHANES.R
library(rstan)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/mass_nonlocal")

# ---- Load and prepare NHANES data ----
# dat <- ...  # Load your NHANES data as a matrix (n x p)
Y   <- scale(dat, center = TRUE, scale = TRUE)
n   <- nrow(Y)
p   <- ncol(Y)
K   <- 1  # Number of factors

# ---- Set hyperparameters ----
a_theta <- 1    # Beta prior for theta
b_theta <- 5
a_sigma <- 1.0  # Gamma prior for sigma2
b_sigma <- 0.3
nu      <- 3    # Degrees of freedom for local shrinkage (MGPS)
a1      <- 2.1  # MGPS prior (delta[1])
a2      <- 3.1  # MGPS prior (delta[2:K])
psi     <- 1    # Dispersion parameter for pMOM
spike_sd <- 1e-4  # Continuous spike std dev

stan_data_j <- list(
   N = n,
   P = p,
   K = K,
   Y = Y,
   a_theta = a_theta,
   b_theta = b_theta,
   a_sigma = a_sigma,
   b_sigma = b_sigma,
   nu = nu,
   a1 = a1,
   a2 = a2,
   psi = psi,
   spike_sd = spike_sd
)

mod <- stan_model("spike_pmom_factor_model.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/applications/mass_nonlocal")

fit_j <- sampling(
   object       = mod,
   data         = stan_data_j,
   chains       = 4,
   iter         = 16000,
   warmup       = 8000,
   seed         = 19,
   init         = "random",
   init_r       = 2,
   control      = list(
      adapt_delta   = 0.99,
      max_treedepth = 15
   )
)

post_j       <- extract(fit_j)
Lambda_j_hat <- apply(post_j$Lambda, c(2,3), mean)

saveRDS(
   list(
      fit        = fit_j,
      posterior  = post_j,
      Lambda_hat = Lambda_j_hat
   ),
   file = "fit_spikepmom_NHANES1718_k1.rds"
)

sum_j      <- summary(fit_j)$summary
max_rhat_j <- max(sum_j[,"Rhat"],   na.rm = TRUE)
min_ess_j  <- min(sum_j[,"n_eff"],  na.rm = TRUE)
bfmi_j     <- rstan::get_bfmi(fit_j)

cat("=== Spike-pMOM model diagnostics (NHANES) ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_j))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_j))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_j, na.rm = TRUE)))
