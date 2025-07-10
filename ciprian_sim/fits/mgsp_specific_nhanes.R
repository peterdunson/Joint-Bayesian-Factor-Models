# ---- Load & Prepare Data ----
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/nhanes")
dat <- readRDS("nhanes_phthalates_adults.rds")
dat <- log1p(dat)  # Optional: log-transform

# Only these 3 variables
mgsp_vars <- c("URXECP", "URXMHH", "URXMOH")
stopifnot(all(mgsp_vars %in% colnames(dat)))
dat3 <- dat[, mgsp_vars, drop = FALSE]

# Z-score scale
Y3 <- scale(dat3, center = TRUE, scale = TRUE)
n  <- nrow(Y3)
p3 <- ncol(Y3)
K  <- 1

# ---- Fit MGSP Model (Stan) ----
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
library(rstan)
mod <- stan_model("mgps_factor_model.stan")

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/nhanes")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan_data <- list(N = n, P = p3, K = K, Y = Y3)
fit <- sampling(
   object       = mod,
   data         = stan_data,
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
post       <- extract(fit)
Lambda_hat <- apply(post$Lambda, c(2,3), mean)   # P3 x K

saveRDS(
   list(
      fit        = fit,
      posterior  = post,
      Lambda_hat = Lambda_hat
   ),
   file = "fit_MGSP_NHANES1718_k1_3var.rds"
)

sum_fit      <- summary(fit)$summary
max_rhat     <- max(sum_fit[,"Rhat"],   na.rm = TRUE)
min_ess      <- min(sum_fit[,"n_eff"],  na.rm = TRUE)
bfmi         <- rstan::get_bfmi(fit)

cat("=== MGSP fit diagnostics ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat))
cat(sprintf("  min n_eff = %.0f\n",  min_ess))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi, na.rm = TRUE)))

