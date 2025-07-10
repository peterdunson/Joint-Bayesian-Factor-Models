# Center and scale (columnwise z-score)
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/nhanes")
dat <- readRDS("nhanes_phthalates_adults.rds")

dat <- log1p(dat)

Y <- scale(dat, center = TRUE, scale = TRUE)
n <- nrow(Y)
p <- ncol(Y)
K <- 1  # Set as desired

# --- Compile Stan model (Dirichlet-Laplace prior) ---
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/dirichlet_laplace")
library(rstan)
mod <- stan_model("dirichlet_laplace.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/method_development")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan_data_dl <- list(N = n, P = p, K = K, Y = Y)

fit_dl <- sampling(
   object       = mod,
   data         = stan_data_dl,
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

post_dl       <- extract(fit_dl)
Lambda_dl_hat <- apply(post_dl$Lambda, c(2,3), mean)

saveRDS(
   list(
      fit        = fit_dl,
      posterior  = post_dl,
      Lambda_hat = Lambda_dl_hat
   ),
   file = "fit_DL_NHANES1718_k1_log.rds"
)

sum_dl      <- summary(fit_dl)$summary
max_rhat_dl <- max(sum_dl[,"Rhat"],   na.rm = TRUE)
min_ess_dl  <- min(sum_dl[,"n_eff"], na.rm = TRUE)
bfmi_dl     <- rstan::get_bfmi(fit_dl)

cat("=== Dirichlet-Laplace model diagnostics ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_dl))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_dl))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_dl, na.rm = TRUE)))
