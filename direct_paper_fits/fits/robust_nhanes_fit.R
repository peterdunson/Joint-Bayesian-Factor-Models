library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(19)

# ---- Load and prepare NHANES data ----
# dat <- ...  # Your NHANES data matrix (n x p)
Y   <- scale(dat, center = TRUE, scale = TRUE)
n   <- nrow(Y)
p   <- ncol(Y)
K   <- 1

# ---- Stan model hyperparameters ----
nu     <- 3    # t-distribution df for robust errors
a1     <- 2.1  # MGP prior hyperparam (delta[1])
a2     <- 3.1  # MGP prior hyperparam (delta[2:K])
kappa  <- 3    # Local shrinkage prior
a_sigma <- 1.0
b_sigma <- 0.3

stan_data_j <- list(
   N = n,
   P = p,
   K = K,
   Y = Y,
   nu = nu,
   a1 = a1,
   a2 = a2,
   kappa = kappa,
   a_sigma = a_sigma,
   b_sigma = b_sigma
)

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/robust_sparse_infinite_factor_model")
mod <- stan_model("robust_sparse_infinite_factor_model.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/applications/robust")

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
   file = "fit_robust_NHANES1718_1.rds"
)

sum_j      <- summary(fit_j)$summary
max_rhat_j <- max(sum_j[,"Rhat"],   na.rm = TRUE)
min_ess_j  <- min(sum_j[,"n_eff"],  na.rm = TRUE)
bfmi_j     <- rstan::get_bfmi(fit_j)

cat("=== Robust model diagnostics (NHANES) ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_j))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_j))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_j, na.rm = TRUE)))

# # ---- Launch ShinyStan (optional) ----
# if (interactive()) {
#    library(shinystan)
#    launch_shinystan(fit_j)
# }
