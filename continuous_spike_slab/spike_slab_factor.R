# --------------------------------------------------------------
# Baseline spike-and-slab factor model (Laplace mixture) for
# sparse Bayesian factor analysis via Stan.
# Follows pipeline for MGPS factor model.
# --------------------------------------------------------------

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# If using Mac, keep this:
Sys.setenv(SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX15.5.sdk")

# ---- CHOOSE SCENARIO ----
scenario <- 1  # Change as needed

# ---- Path to simulation data ----
sim_path <- sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds", scenario)
sim <- readRDS(sim_path)

Y <- sim$Y  # n x p matrix

# ---- Center columns to mean 0 (as required for factor models)
Y <- scale(Y, center = TRUE, scale = FALSE)
n <- nrow(Y)
p <- ncol(Y)

# ---- Stan model setup ----
K <- 5  # Or as appropriate for your data

# ---- Set spike-and-slab prior hyperparameters (adjust as needed)
lambda0 <- 20     # spike (strong shrinkage)
lambda1 <- 0.2    # slab (weak shrinkage)
theta   <- 0.1    # prior inclusion probability

stan_data <- list(
   N = n,
   P = p,
   K = K,
   Y = Y,
   lambda0 = lambda0,
   lambda1 = lambda1,
   theta = theta
)

#set working directory to where the Stan model file is located
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/continuous_spike_slab")

# ---- Compile and sample
mod <- stan_model("spike_slab_factor_model.stan")

# Set working directory to where the simulation data is stored
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations")

fit <- sampling(
   object = mod,
   data = stan_data,
   chains = 4,
   iter = 4000,
   warmup = 2000,
   seed = 42,
   control = list(adapt_delta = 0.99, max_treedepth = 15)
)

print(fit, pars = c("psi", "Lambda"), probs = c(0.1, 0.5, 0.9))

saveRDS(fit, file = sprintf("ssl_factor_fit_scen%d_%d.rds", scenario, K))
