# Center and scale (columnwise z-score)

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/nhanes")
dat <- readRDS("nhanes_phthalates_adults.rds")
dat <- log1p(dat)

Y <- scale(dat, center = TRUE, scale = TRUE)
n <- nrow(Y)
p <- ncol(Y)
K <- 1  # Number of factors

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/horseshoe_estimator")
library(rstan)
mod <- stan_model("horseshoe_factor_model.stan")  # <-- horseshoe model
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/nhanes")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan_data_j <- list(N = n, P = p, K = K, Y = Y)

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
   file = "fit_HS_NHANES1718_k1_log.rds"
)

sum_j      <- summary(fit_j)$summary
max_rhat_j <- max(sum_j[,"Rhat"],   na.rm = TRUE)
min_ess_j  <- min(sum_j[,"n_eff"],  na.rm = TRUE)
bfmi_j     <- rstan::get_bfmi(fit_j)

cat("=== Horseshoe model diagnostics ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_j))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_j))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_j, na.rm = TRUE)))

# If you want to use ShinyStan, uncomment:
# library(shinystan)
# launch_shinystan(fit_j)
