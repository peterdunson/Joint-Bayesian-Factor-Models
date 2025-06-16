# implementation.R
library(rstan)

# Load data
sim <- readRDS("simdata.rds")
Y <- sim$Y
n <- nrow(Y)
p <- ncol(Y)
K <- 20  # Upper bound for number of factors

stan_data <- list(
   N = n,
   P = p,
   K = K,
   Y = Y
)

# Fit the Stan model
fit <- stan(
   file = "mgps_factor_model.stan",
   data = stan_data,
   iter = 2000,
   warmup = 1000,
   chains = 4,
   seed = 42,
   control = list(adapt_delta = 0.95)
)

# Save fit object for future use
saveRDS(fit, file = "stan_fit.rds")
