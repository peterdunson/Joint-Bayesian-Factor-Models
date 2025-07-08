library(rstan)
library(bayesplot)

set.seed(19)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- USE YOUR NHANES DATAFRAME `dat` ----
# `dat` is your data.frame of size n × p
Y <- as.matrix(dat)               # convert to matrix
Y <- scale(Y, center = TRUE, scale = TRUE)  # optionally center & scale
n <- nrow(Y)
p <- ncol(Y)
K <- 5

# ---- COMPILE MGSP FACTOR MODEL ----
mod <- stan_model("/path/to/sparse_bayesian_infinite_factor_model/mgps_factor_model.stan")

# ---- FIT WITH RANDOM INITS ----
stan_data <- list(N = n, P = p, K = K, Y = Y)

fit_nhanes <- sampling(
   object  = mod,
   data    = stan_data,
   chains  = 4,
   iter    = 16000,
   warmup  = 8000,
   seed    = 19,
   init    = "random",
   init_r  = 2,
   control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# ---- EXTRACT & SAVE ----
post        <- extract(fit_nhanes)
Lambda_hat  <- apply(post$Lambda, c(2,3), mean)

saveRDS(
   list(
      fit        = fit_nhanes,
      posterior  = post,
      Lambda_hat = Lambda_hat
   ),
   file = "fit_mgsp_k1_nhanes.rds"
)

# ---- DIAGNOSTICS ----
sum_stats <- summary(fit_nhanes)$summary
cat("=== NHANES MGSP-1 Diagnostics ===\n")
cat(sprintf("  max R̂    = %.3f\n", max(sum_stats[,"Rhat"],   na.rm = TRUE)))
cat(sprintf("  min n_eff = %.0f\n",  min(sum_stats[,"n_eff"], na.rm = TRUE)))
cat(sprintf("  min BFMI  = %.3f\n", min(rstan::get_bfmi(fit_nhanes), na.rm = TRUE)))

