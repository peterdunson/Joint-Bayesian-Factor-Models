library(rstan)
library(bayesplot)

set.seed(19)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- CHOOSE SIMULATION ----
sim_path <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/sim_fixed_lambda_k1_2.rds"
sim      <- readRDS(sim_path)

Y <- sim$X
n <- nrow(Y)
p <- ncol(Y)
K <- 1

# ---- COMPILE MODEL ----
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/horseshoe_estimator")
mod <- stan_model("horseshoe_factor_model.stan")

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/fits")

# ---- Joint fit with random initial values ----
stan_data_j <- list(N = n, P = p, K = K, Y = Y)

fit_j <- sampling(
   object  = mod,
   data    = stan_data_j,
   chains  = 4,
   iter    = 16000,
   warmup  = 8000,
   seed    = 19,
   init    = "random",
   init_r  = 2,
   control = list(
      adapt_delta   = 0.99,
      max_treedepth = 15
   )
)

# ---- Extract & summarize ----
post_j       <- extract(fit_j)
Lambda_j_hat <- apply(post_j$Lambda, c(2,3), mean)

saveRDS(
   list(
      fit        = fit_j,
      posterior  = post_j,
      Lambda_hat = Lambda_j_hat
   ),
   file = "fit_horseshoe_k1_cipsim_2.rds"
)

# ---- Diagnostics ----
sum_j      <- summary(fit_j)$summary
max_rhat_j <- max(sum_j[,"Rhat"],   na.rm = TRUE)
min_ess_j  <- min(sum_j[,"n_eff"],  na.rm = TRUE)
bfmi_j     <- rstan::get_bfmi(fit_j)

cat("=== Horseshoe model diagnostics ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_j))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_j))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_j, na.rm = TRUE)))
