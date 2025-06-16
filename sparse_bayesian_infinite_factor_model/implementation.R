# -----------------------------------------
# Simple implementation of the MGPS factor model from
# Bhattacharya & Dunson (2011) using Stan.
# -----------------------------------------

# 1. Load packages
library(rstan)

# 2. (macOS only) This fixes compiler issues for Stan on Mac.
Sys.setenv(SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX15.5.sdk")

# 3. Load your data
# Expect a file called "simdata.rds" containing a list with an n x p matrix Y
sim <- readRDS("simdata.rds")
Y <- sim$Y  # This should be a numeric matrix: n rows = samples, p columns = variables

# 4. Preprocess data: Center columns to mean 0
Y <- scale(Y, center = TRUE, scale = FALSE)
n <- nrow(Y)
p <- ncol(Y)

# 5. Set Stan model parameters
K <- 30  # Upper bound on number of factors, per paper's guidance

stan_data <- list(
   N = n,
   P = p,
   K = K,
   Y = Y
)

# 6. Compile the Stan model
mod <- stan_model("mgps_factor_model.stan")

# 7. Run MCMC to fit the model
fit <- sampling(
   object = mod,
   data = stan_data,
   chains = 4,
   iter = 2000,
   warmup = 1000,
   seed = 42,
   control = list(adapt_delta = 0.95)
)

# 8. Print a quick summary of results
print(fit, pars = c("tau", "psi", "delta"), probs = c(0.1, 0.5, 0.9))

# 9. Save the fitted model for future use
saveRDS(fit, file = "mgps_fit.rds")
