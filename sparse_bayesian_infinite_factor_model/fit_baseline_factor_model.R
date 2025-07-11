# -----------------------------------------
# Baseline MGPS factor model from Bhattacharya & Dunson (2011) using Stan.
# Updated to use simulation data from scenario directory.
# -----------------------------------------

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# If using Mac, keep this:
Sys.setenv(SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX15.5.sdk")

# ---- CHOOSE SCENARIO ----
scenario <- 2  # Change this to 1, 2, or 3 as needed

# ---- Path to simulation data ----
sim_path <- sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds", scenario) #CHANGE BASED ON NTRAIN
sim <- readRDS(sim_path)

Y <- sim$Y  # n x p matrix

# ---- Center columns to mean 0 (as required for factor models)
Y <- scale(Y, center = TRUE, scale = FALSE)
n <- nrow(Y)
p <- ncol(Y)

# ---- Stan model setup ----
K <- 5  # Or as appropriate for your data

stan_data <- list(
   N = n,
   P = p,
   K = K,
   Y = Y
)

#set working directory to where the Stan model file is located
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")

# ---- Compile and sample
mod <- stan_model("mgps_factor_model.stan")

# Set working directory to where the simulation data is stored
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations")

fit <- sampling(
   object = mod,
   data = stan_data,
   chains = 4, 
   iter = 6000,
   warmup = 3000,
   seed = 12,
   control = list(adapt_delta = 0.999, max_treedepth = 15)
)


print(fit, pars = c("tau", "psi", "delta"), probs = c(0.1, 0.5, 0.9))

saveRDS(fit, file = sprintf("mgps_fit_scen%d_5.rds", scenario))

