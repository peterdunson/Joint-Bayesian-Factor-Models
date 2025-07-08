# ---- 0. Setup ----
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/nhanes")

dat <- readRDS("nhanes_phthalates_adults.rds")
Y <- scale(dat, center = TRUE, scale = TRUE)
n <- nrow(Y)
p <- ncol(Y)
K <- 1  # Number of factors

# ---- 1. Create structured prior SD matrix ----
corY <- cor(Y)
offdiag_idx <- which(upper.tri(corY), arr.ind = TRUE)
cor_vals <- abs(corY[upper.tri(corY)])

# Identify the top 20% highest |correlations|
thresh <- quantile(cor_vals, 0.90)
top_idx <- offdiag_idx[cor_vals >= thresh, , drop = FALSE]
important_vars <- unique(c(top_idx[,1], top_idx[,2]))

# R code to create lambda_prior_sd
lambda_prior_sd <- matrix(1, nrow = p, ncol = K)
lambda_prior_sd[important_vars, ] <- 10

# ---- 2. Stan Data List ----
stan_data_j <- list(
   N = n,
   P = p,
   K = K,
   Y = Y,
   lambda_prior_sd = lambda_prior_sd
)

# ---- 3. Stan Model Fitting ----
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/method_development")
library(rstan)
mod <- stan_model("targeted_shrinkage.stan") # <-- updated stan model with new prior
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/method_development")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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

# ---- 4. Extract Posterior Mean of Lambda ----
post_j       <- extract(fit_j)
Lambda_j_hat <- apply(post_j$Lambda, c(2,3), mean)

# ---- 5. Save Results ----
saveRDS(
   list(
      fit        = fit_j,
      posterior  = post_j,
      Lambda_hat = Lambda_j_hat
   ),
   file = "fit_Joint_NHANES1718_k1_structuredprior.rds"
)

# ---- 6. Diagnostics ----
sum_j      <- summary(fit_j)$summary
max_rhat_j <- max(sum_j[,"Rhat"],   na.rm = TRUE)
min_ess_j  <- min(sum_j[,"n_eff"], na.rm = TRUE)
bfmi_j     <- rstan::get_bfmi(fit_j)

cat("=== Joint model diagnostics ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_j))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_j))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_j, na.rm = TRUE)))
