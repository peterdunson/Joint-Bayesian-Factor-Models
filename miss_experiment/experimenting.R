# fit_and_save.R
library(rstan)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- CHOOSE SCENARIO ----
scenario <- 2
sim_path <- sprintf(
   "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds",
   scenario
)
sim <- readRDS(sim_path)

# center & scale to unit variance
Y   <- scale(sim$Y, center = TRUE, scale = TRUE)
X   <- Y[, -1]
n   <- nrow(Y)
p   <- ncol(Y)
p_x <- ncol(X)
K   <- 5

# compile the non-centered model
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
mod <- stan_model("mgps_factor_model.stan")

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment")

# ---- 1) X-only fit ----
stan_data_x <- list(N = n, P = p_x, K = K, Y = X)
fit_x       <- sampling(
   mod, data    = stan_data_x,
   chains  = 4, iter = 6000, warmup = 3000,
   seed    = 12,
   control = list(adapt_delta=0.99, max_treedepth=15)
)

# extract & save
post_x       <- extract(fit_x)
Lambda_x_hat <- apply(post_x$Lambda, c(2,3), mean)
saveRDS(
   list(post       = post_x,
        Lambda_hat = Lambda_x_hat),
   file = sprintf("fit_Xonly_scen%d_scale.rds", scenario)
)

# ---- Diagnostics: X-only ----
sum_x      <- summary(fit_x)$summary
max_rhat_x <- max(sum_x[,"Rhat"], na.rm=TRUE)
min_ess_x  <- min(sum_x[,"n_eff"], na.rm=TRUE)
bfmi_x     <- bfmi(fit_x)

cat("\nX-only diagnostics:\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_x))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_x))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_x, na.rm=TRUE)))



# ---- 2) Joint fit ----
stan_data_j <- list(N = n, P = p,   K = K, Y = Y)
fit_j       <- sampling(
   mod, data    = stan_data_j,
   chains  = 4, iter = 6000, warmup = 3000,
   seed    = 12,
   control = list(adapt_delta=0.99, max_treedepth=15)
)

# extract & save
post_j       <- extract(fit_j)
Lambda_j_hat <- apply(post_j$Lambda, c(2,3), mean)
saveRDS(
   list(post       = post_j,
        Lambda_hat = Lambda_j_hat),
   file = sprintf("fit_Joint_scen%d_scale.rds", scenario)
)

# ---- Diagnostics: Joint ----
sum_j      <- summary(fit_j)$summary
max_rhat_j <- max(sum_j[,"Rhat"], na.rm=TRUE)
min_ess_j  <- min(sum_j[,"n_eff"], na.rm=TRUE)
bfmi_j     <- bfmi(fit_j)

cat("\nJoint model diagnostics:\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_j))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_j))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_j, na.rm=TRUE)))

