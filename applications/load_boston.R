# load_boston.R

#Boston Housing506 neighborhoods × 13 socioeconomic and housing variables

# ─── 1) Load & Prepare Boston Data ────────────────────────────────────────
# install.packages("MASS")
library(MASS)

data("Boston", package = "MASS")
df_boston       <- Boston[, setdiff(names(Boston), "medv")]
df_boston_scaled <- scale(df_boston)

# ─── 2) Set up for Stan ───────────────────────────────────────────────────
library(rstan)
library(bayesplot)

set.seed(19)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

Y <- df_boston_scaled
n <- nrow(Y)
p <- ncol(Y)
K <- 3   # choose number of factors

# ─── 3) Compile MGPS Stan Model ──────────────────────────────────────────
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
mod <- stan_model("mgps_factor_model.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment")

# ─── 4) Fit with Random Inits ─────────────────────────────────────────────
stan_data <- list(N = n, P = p, K = K, Y = Y)

fit_boston <- sampling(
   object   = mod,
   data     = stan_data,
   chains   = 4,
   iter     = 12000,
   warmup   = 6000,
   seed     = 19,
   init     = "random",
   init_r   = 2,
   control  = list(adapt_delta = 0.99, max_treedepth = 15)
)

# ─── 5) Extract, Save & Diagnose ──────────────────────────────────────────
post      <- extract(fit_boston)
Lambda_hat <- apply(post$Lambda, c(2,3), mean)

saveRDS(
   list(fit         = fit_boston,
        posterior   = post,
        Lambda_hat  = Lambda_hat),
   file = "fit_Boston_mgps.rds"
)

sum_stats  <- summary(fit_boston)$summary
cat("max R̂    =", max(sum_stats[,"Rhat"], na.rm=TRUE), "\n")
cat("min n_eff =", min(sum_stats[,"n_eff"], na.rm=TRUE), "\n")
cat("min BFMI  =", min(rstan::get_bfmi(fit_boston), na.rm=TRUE), "\n")

# ─── 6) (Optional) Launch ShinyStan ───────────────────────────────────────
# library(shinystan)
# launch_shinystan(fit_boston)
