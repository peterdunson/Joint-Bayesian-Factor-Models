# fit_tebfar_factor_model.R
# --------------------------------------------------------
# Fits the TEB-FAR model (Palmer & Dunson, 2025+) via Stan.
# The outcome variance (first variable) is fixed to user-provided Sigma1.
# --------------------------------------------------------

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

Sys.setenv(SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX15.5.sdk")

# 1. Load & preprocess data ---------------------------------------------
sim  <- readRDS("simdata.rds")
Y    <- sim$Y               # (n x p) data matrix; first column = outcome y
n    <- nrow(Y)
p    <- ncol(Y)

# Center the data
Y <- scale(Y, center = TRUE, scale = FALSE)

# 2. Set model parameters -----------------------------------------------
K      <- 30         # upper bound for number of factors
Sigma1 <- 1.0        # <- set this to your chosen outcome variance

# 3. Prepare Stan data --------------------------------------------------
stan_data <- list(
   N      = n,
   P      = p,
   K      = K,
   Y      = Y,
   Sigma1 = Sigma1
)

# 4. Fit the TEB-FAR model ----------------------------------------------
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

# 5. Save results -------------------------------------------------------
saveRDS(fit, "stan_tebfar_fit.rds")

# 6. (Optional) Quick summary
print(fit, pars = c("psi", "tau"), probs = c(0.1, 0.5, 0.9))
