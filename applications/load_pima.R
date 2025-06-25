# load_pima.R

#Pima Indians Diabetes 
#392 women × 9 variables: 8 numeric clinical measurements+ diabetes.

# ─── 1) Load & Prepare Pima Indians Diabetes Data ─────────────────────────
library(mlbench)

data("PimaIndiansDiabetes", package = "mlbench")
df_pima <- na.omit(PimaIndiansDiabetes)

# Separate outcome and predictors, then standardize predictors
outcome <- df_pima$diabetes
Y       <- scale(df_pima[ , setdiff(names(df_pima), "diabetes")])
n       <- nrow(Y)
p       <- ncol(Y)
K       <- 3   # choose number of latent factors

# ─── 2) Compile MGPS Stan Model ───────────────────────────────────────────
library(rstan)
library(bayesplot)

set.seed(23)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
mod <- stan_model("mgps_factor_model.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment")

# ─── 3) Fit with Random Initial Values ────────────────────────────────────
stan_data <- list(N = n, P = p, K = K, Y = Y)

fit_pima <- sampling(
   object       = mod,
   data         = stan_data,
   chains       = 4,
   iter         = 6000,
   warmup       = 3000,
   seed         = 23,
   init         = "random",
   init_r       = 2,
   control      = list(adapt_delta = 0.99, max_treedepth = 15)
)

# ─── 4) Extract Posterior & Compute Posterior-Mean Loadings ───────────────
post_pima   <- extract(fit_pima)
Lambda_hat  <- apply(post_pima$Lambda, c(2,3), mean)

# ─── 5) Save Results ──────────────────────────────────────────────────────
saveRDS(
   list(
      fit         = fit_pima,
      posterior   = post_pima,
      Lambda_hat  = Lambda_hat
   ),
   file = "fit_Pima_mgps.rds"
)

# ─── 6) Diagnostics ────────────────────────────────────────────────────────
sum_stats  <- summary(fit_pima)$summary
cat("max R̂    =", max(sum_stats[,"Rhat"], na.rm=TRUE), "\n")
cat("min n_eff =", min(sum_stats[,"n_eff"], na.rm=TRUE), "\n")
cat("min BFMI  =", min(rstan::get_bfmi(fit_pima), na.rm=TRUE), "\n")

# ─── 7) (Optional) ShinyStan ───────────────────────────────────────────────
# library(shinystan)
# launch_shinystan(fit_pima)


