# fit_and_save_randominit.R

library(rstan)
library(bayesplot)

set.seed(19)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- CHOOSE SCENARIO ----
scenario <- 2
sim_path <- sprintf(
   "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000_1.rds",
   scenario
)
sim <- readRDS(sim_path)

# ---- CENTER & SCALE DATA ----
Y   <- scale(sim$Y,   center = TRUE, scale = TRUE)
n   <- nrow(Y)
p   <- ncol(Y)
K   <- 1

# ---- COMPILE MODEL ----
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
mod <- stan_model("mgps_factor_model.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment")

# ---- 2) Joint fit with random initial values ----
stan_data_j <- list(N = n, P = p, K = K, Y = Y)

fit_j <- sampling(
   object       = mod,
   data         = stan_data_j,
   chains       = 4,
   iter         = 12000,
   warmup       = 6000,
   seed         = 19,
   init         = "random",    # random inits
   init_r       = 2,           # uniform(−2,2)
   control      = list(
      adapt_delta   = 0.99,
      max_treedepth = 15
   )
)

# ---- 3) Extract & summarize ----
post_j       <- extract(fit_j)
Lambda_j_hat <- apply(post_j$Lambda, c(2,3), mean)

saveRDS(
   list(
      fit        = fit_j,
      posterior  = post_j,
      Lambda_hat = Lambda_j_hat
   ),
   file = sprintf("fit_Joint_scen%d_scale_all_randominit_3.rds", scenario)
)

# ---- 4) Diagnostics ----
sum_j      <- summary(fit_j)$summary
max_rhat_j <- max(sum_j[,"Rhat"],   na.rm = TRUE)
min_ess_j  <- min(sum_j[,"n_eff"], na.rm = TRUE)
bfmi_j     <- rstan::get_bfmi(fit_j)

cat("=== Joint model diagnostics ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_j))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_j))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_j, na.rm = TRUE)))

# ---- 5) Launch ShinyStan ----
library(shinystan)
launch_shinystan(fit_j)

