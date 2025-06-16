# fit_tebfar_factor_model.R
# --------------------------------------------------------
# Fits the TEB-FAR model (Palmer & Dunson, 2025+) via Stan.
# The outcome variance (first variable) is fixed to user-provided Sigma1.
# --------------------------------------------------------

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Only needed for Mac:
Sys.setenv(SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX15.5.sdk")

# ---- CHOOSE SCENARIO ----
scenario <- 1   # Change this to 1, 2, or 3 as needed

# ---- Path to simulation data ----
sim_path <- sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d.rds", scenario)
sim <- readRDS(sim_path)

Y <- sim$Y  # (n x p+1) matrix, first column is y
n <- nrow(Y)
p <- ncol(Y)

# Center columns (outcome and predictors)
Y <- scale(Y, center = TRUE, scale = FALSE)

# ---- Set Stan model parameters ----
K      <- 30         # upper bound for number of factors
Sigma1 <- 1.0        # fixed variance for outcome (can be set as desired)

stan_data <- list(
   N      = n,
   P      = p,
   K      = K,
   Y      = Y,
   Sigma1 = Sigma1
)

# ---- Compile and fit Stan model ----
mod <- stan_model("tebfar_factor_model.stan")

fit <- sampling(
   object  = mod,
   data    = stan_data,
   chains  = 4,
   iter    = 2000,
   warmup  = 1000,
   seed    = 42,
   control = list(adapt_delta = 0.95)
)

# ---- Save results ----
saveRDS(fit, sprintf("stan_tebfar_fit_scen%d.rds", scenario))

# ---- (Optional) Quick summary ----
print(fit, pars = c("psi", "tau"), probs = c(0.1, 0.5, 0.9))

