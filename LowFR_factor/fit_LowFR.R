# --------------------------------------------------------------
# Low-Rank Longitudinal Factor Regression (LowFR) model
# Palmer, Herring & Dunson (2023, arXiv)
# Follows pipeline for MGPS, SSL, Horseshoe, etc.
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

X <- sim$X  # n x (p*TT) matrix (vectorized exposures)
y <- sim$y  # length n outcome

# ---- Center columns to mean 0 (as required for factor models)
X <- scale(X, center = TRUE, scale = FALSE)
y <- as.vector(scale(y, center = TRUE, scale = FALSE))

n <- nrow(X)
p <- sim$p      # Number of exposures (or set manually if not in sim)
TT <- sim$TT    # Number of time points (or set manually)
K <- 5          # Number of latent factors (set as needed)
H <- min(p, TT) # Rank hyperparameter (paper default)

# ---- Stan model setup ----
stan_data <- list(
   N = n,
   p = p,
   k = K,
   TT = TT,
   H = H,
   y = y,
   X = X
)

# Set working directory to where the Stan model file is located
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/LowFR_factor")

# ---- Compile and sample ----
mod <- stan_model("LowFR.stan")

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

print(fit, pars = c("alpha_0", "sigma2"), probs = c(0.1, 0.5, 0.9))

saveRDS(fit, file = sprintf("lowfr_fit_scen%d_%d.rds", scenario, K))


