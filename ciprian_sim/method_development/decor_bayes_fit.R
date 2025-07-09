library(rstan)
library(bayesplot)

set.seed(19)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/nhanes")

dat <- readRDS("nhanes_phthalates_adults.rds")

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/method_development")

# ---- Data ----
Y <- as.matrix(dat)
Y <- scale(Y, center = TRUE, scale = TRUE)
n <- nrow(Y)
p <- ncol(Y)
K <- 1 # or any K you want

# ---- Construct Data-Aware Priors ----
corY <- cor(Y)
eig  <- eigen(corY)
lambda_mu  <- matrix(0, nrow = p, ncol = K)
lambda_tau <- matrix(1, nrow = p, ncol = K)
for (k in 1:K) {
   lambda_mu[, k] <- eig$vectors[, k]
   # Optionally: make lambda_tau larger for higher prior flexibility
}

stan_data <- list(
   N = n,
   P = p,
   K = K,
   Y = Y,
   lambda_mu = lambda_mu,
   lambda_tau = lambda_tau
)

# ---- Compile and Fit Model ----
mod <- stan_model("decor_bayes.stan")

fit <- sampling(
   object  = mod,
   data    = stan_data,
   chains  = 4,
   iter    = 16000,
   warmup  = 8000,
   seed    = 19,
   control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# ---- Extract Posterior ----
post        <- extract(fit)
Lambda_hat  <- apply(post$Lambda, c(2,3), mean)

saveRDS(
   list(
      fit        = fit,
      posterior  = post,
      Lambda_hat = Lambda_hat
   ),
   file = sprintf("fit_decor_k%d_nhanes.rds", K)
)

# ---- Diagnostics ----
sum_stats <- summary(fit)$summary
cat(sprintf("=== Decorrelating Factor Model (K=%d) Diagnostics ===\n", K))
cat(sprintf("  max R̂    = %.3f\n", max(sum_stats[,"Rhat"],   na.rm = TRUE)))
cat(sprintf("  min n_eff = %.0f\n",  min(sum_stats[,"n_eff"], na.rm = TRUE)))
cat(sprintf("  min BFMI  = %.3f\n", min(rstan::get_bfmi(fit), na.rm = TRUE)))
