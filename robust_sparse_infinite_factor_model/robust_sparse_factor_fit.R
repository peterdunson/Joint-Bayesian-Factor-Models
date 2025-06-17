# fit_robust_sparse_factor_model.R
# --------------------------------------------------------
# Fits the robust sparse Bayesian infinite factor model (Student-t errors + MGPS prior) via Stan
# --------------------------------------------------------

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Only needed for Mac:
Sys.setenv(SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX15.5.sdk")

# ---- CHOOSE SCENARIO ----
scenario <- 1   # Change as needed

# ---- Path to simulation data ----
sim_path <- sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds", scenario)
sim <- readRDS(sim_path)

Y <- sim$Y  # (n x p+1) matrix or (n x p), first column y if present
n <- nrow(Y)
p <- ncol(Y)

# Center columns (as in TEB-FAR)
Y <- scale(Y, center = TRUE, scale = FALSE)

# ---- Set Stan model parameters ----******************************************
K        <- 5     # Truncation level (upper bound for factors)
nu       <- 3     # Degrees of freedom for Student-t
a1       <- 2
a2       <- 3
kappa    <- 2
a_sigma  <- 2
b_sigma  <- 2

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

# set working directory to where the Stan model file is located
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/robust_sparse_infinite_factor_model")

mod <- stan_model("robust_sparse_infinite_factor_model.stan")

# set working directory to where the simulation data is stored (for output)
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

# ---- Save results ----
saveRDS(fit, sprintf("robust_sparse_factor_fit_scen%d_5.rds", scenario))

# ---- (Optional) Quick summary ----
print(fit, pars = c("sigma2", "tau"), probs = c(0.1, 0.5, 0.9))

