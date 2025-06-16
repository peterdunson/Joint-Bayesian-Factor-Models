# mgps_factor_analysis.R
# --------------------------------------------------------
# Implements the multiplicative gamma process shrinkage (MGPS)
# factor model from Bhattacharya & Dunson (2011) via Stan.
# --------------------------------------------------------

# 1. Setup -------------------------------------------------------------------
library(rstan)

# Allow Stan to reuse compiled models and run chains in parallel
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# On macOS with CommandLineTools only:
Sys.setenv(SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX15.5.sdk")


# 2. Load & preprocess data -------------------------------------------------
# simdata.rds should contain a list with element Y: an (n × p) matrix of observations
sim  <- readRDS("simdata.rds")
Y    <- sim$Y               # our data matrix
n    <- nrow(Y)             # number of samples
p    <- ncol(Y)             # number of observed variables

# MGPS priors work best when Y is centered (zero mean) across rows:
Y <- scale(Y, center = TRUE, scale = FALSE)


# 3. Prepare data for Stan --------------------------------------------------
# N: number of observations
# P: dimension of each observation
# K: upper bound on latent factors (choose conservatively)
K <- 30

stan_data <- list(
   N = n,
   P = p,
   K = K,
   Y = Y
)


# 4. Compile & fit the model ------------------------------------------------
# mgps_factor_model.stan implements the Bayesian infinite factor model with
# multiplicative gamma process shrinkage on the loadings.
mod <- stan_model("mgps_factor_model.stan")

fit <- sampling(
   object      = mod,
   data        = stan_data,
   chains      = 4,
   iter        = 2000,     # total iterations per chain
   warmup      = 1000,     # iterations used for adaptation (discarded)
   seed        = 42,
   control     = list(adapt_delta = 0.95)
)

# Quick check: print a summary of key parameters
print(fit, pars = c("tau", "phi"), probs = c(0.1, 0.5, 0.9))


# 5. Save results -----------------------------------------------------------
# Persist the fitted object for downstream analysis (e.g., posterior predictive checks)
saveRDS(fit, "stan_fit.rds")

