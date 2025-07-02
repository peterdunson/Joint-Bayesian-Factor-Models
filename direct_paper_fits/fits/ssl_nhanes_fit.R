# --------------------------------------------------------------
# Baseline spike-and-slab factor model for sparse Bayesian FA via Stan
# Follows the same pipeline as your other fits (center/scale, diagnostics, save list)
# --------------------------------------------------------------

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- DATA: Prepare & scale ----
dat <- # <-- load or assign your data matrix here
   Y <- scale(dat, center = TRUE, scale = TRUE)
n <- nrow(Y)
p <- ncol(Y)
K <- 2  

# ---- MODEL ----
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/continuous_spike_slab")
mod <- stan_model("spike_slab_factor_model.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/direct_paper_fits/storing_fit/spike_slab")

# ---- HYPERPARAMETERS ----
lambda0 <- 20     # spike (strong shrinkage)
lambda1 <- 0.2    # slab (weak shrinkage)
theta   <- 0.1    # prior inclusion probability

stan_data_j <- list(
   N = n,
   P = p,
   K = K,
   Y = Y,
   lambda0 = lambda0,
   lambda1 = lambda1,
   theta   = theta
)

# ---- FIT ----
fit_j <- sampling(
   object       = mod,
   data         = stan_data_j,
   chains       = 4,
   iter         = 16000,
   warmup       = 8000,
   seed         = 19,
   init         = "random",
   init_r       = 2,
   control      = list(
      adapt_delta   = 0.99,
      max_treedepth = 15
   )
)

post_j       <- extract(fit_j)
Lambda_j_hat <- apply(post_j$Lambda, c(2,3), mean)

saveRDS(
   list(
      fit        = fit_j,
      posterior  = post_j,
      Lambda_hat = Lambda_j_hat
   ),
   file = "fit_SSL_NHANES1718_k2.rds"
)

sum_j      <- summary(fit_j)$summary
max_rhat_j <- max(sum_j[,"Rhat"],   na.rm = TRUE)
min_ess_j  <- min(sum_j[,"n_eff"],  na.rm = TRUE)
bfmi_j     <- rstan::get_bfmi(fit_j)

cat("=== Spike-and-slab model diagnostics ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_j))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_j))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_j, na.rm = TRUE)))

# # If you want ShinyStan:
# library(shinystan)
# launch_shinystan(fit_j)
