# fit_and_save.R
library(rstan)
library(bayesplot)

set.seed(19)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- CHOOSE SCENARIO ----
scenario <- 2
sim_path <- sprintf( "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds",
   scenario
)
sim <- readRDS(sim_path)

# ---- CENTER & SCALE DATA ----
Y   <- scale(sim$Y, center=TRUE, scale=TRUE)
X   <- Y[, -1]
n   <- nrow(Y)
p   <- ncol(Y)
p_x <- ncol(X)
K   <- 5
K_max <- 5

# ---- COMPILE MODEL ----
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
mod <- stan_model("mgps_factor_model.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment")



# ---- 2) Joint fit ----
stan_data_j <- list(N=n, P=p, K=K, Y=Y)
fit_j <- sampling(
   mod,
   data       = stan_data_j,
   chains     = 4,
   iter       = 6000,
   warmup     = 3000,
   seed       = 12,
   init       = "random",
   init_r     = 5,    # ← Uniform(−5,5) inits
   control    = list(
      adapt_delta   = 0.99,
      max_treedepth = 15
   )
)

# extract & summarize
post_j       <- extract(fit_j)
Lambda_j_hat <- apply(post_j$Lambda, c(2,3), mean)

# one file, everything in it
saveRDS(
   list(
      fit        = fit_j,
      posterior  = post_j,
      Lambda_hat = Lambda_j_hat
   ),
   file = sprintf("fit_Joint_scen%d_scale_all_15k.rds", scenario)
)

# ---- Diagnostics: Joint ----
sum_j      <- summary(fit_j)$summary
max_rhat_j <- max(sum_j[,"Rhat"],   na.rm=TRUE)
min_ess_j  <- min(sum_j[,"n_eff"],  na.rm=TRUE)
bfmi_j     <- rstan::get_bfmi(fit_j)

cat("=== Joint model diagnostics ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_j))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_j))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_j, na.rm=TRUE)))



library(shinystan)

# Launch shinyStan (this will open a web app in your browser)
launch_shinystan(fit_j)
