# fit_bass_factor_model.R
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Data setup: centered
scenario <- 1
sim_path <- sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds", scenario)
sim <- readRDS(sim_path)
Y <- sim$Y
Y <- scale(Y, center=TRUE, scale=FALSE)
n <- nrow(Y)
p <- ncol(Y)
K <- 5

stan_data <- list(
   N = n,
   P = p,
   K = K,
   Y = Y
)

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/bass")

mod <- stan_model("bass_factor.stan")

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

saveRDS(fit, sprintf("stan_bass_fit_scen%d_%d.rds", scenario, K))
print(fit, pars = c("sigma2", "tau_global", "tau_factor"), probs = c(0.1, 0.5, 0.9))
