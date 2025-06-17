# fit_horseshoe_factor_model.R
# --------------------------------------------------------
# Fits the horseshoe factor model (Carvalho et al., 2010) via Stan.
# The outcome variance (first variable) is fixed to user-provided Sigma1.
# --------------------------------------------------------

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Only needed for Mac (change if needed or comment out if not on Mac):
Sys.setenv(SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX15.5.sdk")

# ---- CHOOSE SCENARIO ----
scenario <- 1   # Change this as needed

# ---- Path to simulation data ----
sim_path <- sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds", scenario)
sim <- readRDS(sim_path)

Y <- sim$Y
n <- nrow(Y)
p <- ncol(Y)

Y <- scale(Y, center = TRUE, scale = FALSE) # Center columns

K <- 5           # truncation for factors
Sigma1 <- 1.0    # fixed variance for outcome

stan_data <- list(
   N = n,
   P = p,
   K = K,
   Y = Y,
   Sigma1 = Sigma1
)

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/horseshoe_estimator")

mod <- stan_model("horseshoe_factor_model.stan")

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

saveRDS(fit, sprintf("stan_horseshoe_fit_scen%d_5.rds", scenario))
print(fit, pars = c("psi", "lambda_global", "lambda_local"), probs = c(0.1, 0.5, 0.9))

