# fit_tebfar_model_updated.R

library(rstan)
library(bayesplot)

set.seed(19)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- CHOOSE SCENARIO ----
scenario <- 2
sim_path <- sprintf(
  "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds",
  scenario
)
sim <- readRDS(sim_path)

# ---- EXTRACT & PREPARE DATA ----
# sim$Y is [n × (p+1)], first column is the outcome y
Y_raw <- sim$Y
n     <- nrow(Y_raw)
P     <- ncol(Y_raw)

# center (but do not scale) so that outcome variance is fixed correctly
Y_cent <- scale(Y_raw, center = TRUE, scale = FALSE)

# ---- MODEL SETTINGS ----
K      <- 5       # truncation / upper bound on # factors
Sigma1 <- 1.0     # fixed outcome variance

stan_data <- list(
  N      = n,
  P      = P,
  K      = K,
  Y      = Y_cent,
  Sigma1 = Sigma1
)

# ---- COMPILE & FIT ----
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/teb_far")
mod <- stan_model("tebfar_factor_model.stan")

fit <- sampling(
  mod,
  data    = stan_data,
  chains  = 4,
  iter    = 6000,
  warmup  = 3000,
  seed    = 42,
  init    = "random",
  init_r  = 5,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# ---- SAVE RESULTS ----
saveRDS(
  list(
    fit       = fit,
    stan_data = stan_data
  ),
  file = sprintf("stan_tebfar_fit_scen%d_K%d.rds", scenario, K)
)

# ---- QUICK DIAGNOSTICS ----
sum_stats  <- summary(fit)$summary
max_rhat   <- max(sum_stats[,"Rhat"],   na.rm = TRUE)
min_ess    <- min(sum_stats[,"n_eff"],  na.rm = TRUE)
min_bfmi   <- min(rstan::get_bfmi(fit), na.rm = TRUE)

cat("=== TEB-FAR model diagnostics ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat))
cat(sprintf("  min n_eff = %.0f\n",  min_ess))
cat(sprintf("  min BFMI  = %.3f\n", min_bfmi))
