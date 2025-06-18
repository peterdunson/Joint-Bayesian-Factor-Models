# --------------------------------------------------------------
# Robust Sparse Bayesian Infinite Factor Model (t errors, MGP)
# Lee, Jo, Lee (2022, Comp Stat)
# Follows pipeline for MGPS, SSL, etc.
# --------------------------------------------------------------

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# For Mac only (if needed)
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
K <- 5        # Truncation level for number of factors (set as appropriate)
nu <- 3       # Degrees of freedom for t-distribution (default 3 or 5)
a1 <- 2.1     # MGP prior hyperparameter (delta[1])
a2 <- 3.1     # MGP prior hyperparameter (delta[2:K])
kappa <- 3    # Local shrinkage prior (usually 2 or 3)
a_sigma <- 1.0
b_sigma <- 0.3

stan_data <- list(
   N = n,
   P = p,
   K = K,
   Y = Y,
   nu = nu,
   a1 = a1,
   a2 = a2,
   kappa = kappa,
   a_sigma = a_sigma,
   b_sigma = b_sigma
)

#set working directory to where the Stan model file is located
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/robust_sparse_infinite_factor_model")

# ---- Compile and sample
mod <- stan_model("robust_sparse_infinite_factor_model.stan")

# Set working directory to where the simulation data is stored
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations")

fit <- sampling(
   object = mod,
   data = stan_data,
   chains = 4,
   iter = 3000,
   warmup = 1500,
   seed = 42,
   control = list(adapt_delta = 0.99, max_treedepth = 15)
)

print(fit, pars = c("sigma2", "Lambda"), probs = c(0.1, 0.5, 0.9))

saveRDS(fit, file = sprintf("robust_factor_fit_scen%d_%d.rds", scenario, K))
