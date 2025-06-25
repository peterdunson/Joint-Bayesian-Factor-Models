# load_wine.R

#Wine178 wine samples × 13 chemical measurements

# ─── 1) Load & Prepare Wine Data ──────────────────────────────────────────
library(rattle)

data(wine, package = "rattle")
df_wine <- wine[, setdiff(names(wine), "Type")]  # drop the Class/Type column

# Standardize (center + scale)
Y <- scale(df_wine, center = TRUE, scale = TRUE)
n <- nrow(Y)
p <- ncol(Y)
K <- 3   # choose number of latent factors

# ─── 2) Compile MGPS Stan Model ───────────────────────────────────────────
library(rstan)
library(bayesplot)

set.seed(22)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
mod <- stan_model("mgps_factor_model.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment")

# ─── 3) Fit with Random Initial Values ────────────────────────────────────
stan_data <- list(N = n, P = p, K = K, Y = Y)

fit_wine <- sampling(
   object       = mod,
   data         = stan_data,
   chains       = 4,
   iter         = 6000,
   warmup       = 3000,
   seed         = 22,
   init         = "random",
   init_r       = 2,
   control      = list(adapt_delta = 0.99, max_treedepth = 15)
)

# ─── 4) Extract Posterior & Compute Posterior-Mean Loadings ───────────────
post_wine    <- extract(fit_wine)
Lambda_hat   <- apply(post_wine$Lambda, c(2,3), mean)

# ─── 5) Save Results ──────────────────────────────────────────────────────
saveRDS(
   list(
      fit         = fit_wine,
      posterior   = post_wine,
      Lambda_hat  = Lambda_hat
   ),
   file = "fit_Wine_mgps.rds"
)

# ─── 6) Diagnostics ────────────────────────────────────────────────────────
sum_stats   <- summary(fit_wine)$summary
cat("max R̂    =", max(sum_stats[,"Rhat"], na.rm=TRUE), "\n")
cat("min n_eff =", min(sum_stats[,"n_eff"], na.rm=TRUE), "\n")
cat("min BFMI  =", min(rstan::get_bfmi(fit_wine), na.rm=TRUE), "\n")

# ─── 7) (Optional) ShinyStan ───────────────────────────────────────────────
# library(shinystan)
# launch_shinystan(fit_wine)
