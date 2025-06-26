# -------------------------------------------------------
# Toy Example: 3-Factor Model with 9 Variables
# - Simulates data at low and high noise
# - Fits Stan 3-factor model
# - Plots correlation/residual histograms and Fisher-z
# - Computes DSC (mean, sd, skew, kurtosis) on Fisher-z
# - Launches ShinyStan for each fit
# -------------------------------------------------------

library(rstan)
library(shinystan)

# ----- Simulation -----
set.seed(42)
n <- 800     # Number of samples
p <- 10      # Number of variables (3 per factor, could be more)
K_true <- 3  # True number of factors

# Strong loadings, structure: each factor loads strongly on 3 variables
Lambda <- matrix(0, nrow=p, ncol=K_true)
for (k in 1:K_true) {
   Lambda[((k-1)*3 + 1):(k*3), k] <- rnorm(3, mean=2, sd=0.4)
}
# Add some weak cross-loadings if you want more realism
# Lambda[Lambda == 0] <- rnorm(sum(Lambda==0), mean=0, sd=0.2)

Eta <- matrix(rnorm(n*K_true), nrow=n, ncol=K_true)
Y_signal <- Eta %*% t(Lambda)  # n x p

noise_levels <- c(0.3, 1.2)  # Low and high noise
all_data <- list()
for (noise in noise_levels) {
   Y <- Y_signal + matrix(rnorm(n*p, sd=noise), nrow=n, ncol=p)
   Y <- scale(Y, center=TRUE, scale=TRUE)
   all_data[[as.character(noise)]] <- Y
}

# ----- Stan Model -----
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
mod <- stan_model("mgps_factor_model.stan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fits <- list()
residuals_list <- list()
par(mfrow=c(2, 4))

for (i in seq_along(noise_levels)) {
   noise <- noise_levels[i]
   cat("\nNoise SD:", noise, "\n")
   Y <- all_data[[as.character(noise)]]
   
   # --- Stan fit (3-factor) ---
   stan_data <- list(N=n, P=p, K=3, Y=Y)
   fit <- sampling(
      object = mod,
      data = stan_data,
      chains = 4,
      iter = 10000,
      warmup = 5000,
      init = "random",
      init_r = 2,  # uniform(-2, 2)
      control = list(adapt_delta = 0.99, max_treedepth = 15)
   )
   fits[[i]] <- fit
   
   # --- Posterior mean estimates ---
   post <- extract(fit)
   Lambda_hat <- apply(post$Lambda, c(2,3), mean)
   eta_hat <- Y %*% Lambda_hat %*% solve(t(Lambda_hat) %*% Lambda_hat)
   Y_hat <- eta_hat %*% t(Lambda_hat)
   resid_mat <- Y - Y_hat
   residuals_list[[i]] <- resid_mat
   
   # --- Correlation histograms ---
   R_obs <- cor(Y)
   R_resid <- cor(resid_mat)
   off_diag_obs <- R_obs[lower.tri(R_obs)]
   off_diag_resid <- R_resid[lower.tri(R_resid)]
   
   hist(off_diag_obs, breaks=15, main=paste0("Obs Corr (noise=", noise, ")"), xlab="Correlation")
   hist(off_diag_resid, breaks=15, main=paste0("Resid Corr (noise=", noise, ")"), xlab="Correlation")
   
   # --- Fisher z-transform ---
   fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
   z_obs <- fisher_z(off_diag_obs)
   z_resid <- fisher_z(off_diag_resid)
   hist(z_obs, breaks=15, main=paste0("Fisher-z Obs (", noise, ")"), xlab="z(corr)")
   hist(z_resid, breaks=15, main=paste0("Fisher-z Resid (", noise, ")"), xlab="z(corr)")
   
   # --- DSC (mean, sd, skew, kurtosis) ---
   dsc <- function(corrs, mu2 = 0, sd2 = NULL, sk2 = 0, ku2 = 3, n = NULL) {
      mu1 <- mean(corrs)
      sd1 <- sd(corrs)
      if (sd1 < 1e-12) {
         sk1 <- 0; ku1 <- 0
      } else {
         sk1 <- mean((corrs - mu1)^3) / sd1^3
         ku1 <- mean((corrs - mu1)^4) / sd1^4
      }
      if (is.null(sd2)) {
         if (is.null(n)) stop("Need n (sample size) if sd2 is not given")
         sd2 <- 1 / sqrt(n - 3)
      }
      sk_diff <- (sign(sk1) * abs(sk1)^(1/3)) - (sign(sk2) * abs(sk2)^(1/3))
      ku_diff <- (sign(ku1) * abs(ku1)^(1/4)) - (sign(ku2) * abs(ku2)^(1/4))
      D <- sqrt((mu1 - mu2)^2 + (sd1 - sd2)^2 + sk_diff^2 + ku_diff^2)
      return(list(DSC = D, mu1 = mu1, sd1 = sd1, sk1 = sk1, ku1 = ku1))
   }
   dsc_obs <- dsc(z_obs, n=n)
   dsc_resid <- dsc(z_resid, n=n)
   cat(sprintf(
      "DSC obs (noise %.1f): %.3f (mean=%.3f, sd=%.3f, skew=%.3f, kurtosis=%.3f)\n",
      noise, dsc_obs$DSC, dsc_obs$mu1, dsc_obs$sd1, dsc_obs$sk1, dsc_obs$ku1))
   cat(sprintf(
      "DSC resid (noise %.1f): %.3f (mean=%.3f, sd=%.3f, skew=%.3f, kurtosis=%.3f)\n",
      noise, dsc_resid$DSC, dsc_resid$mu1, dsc_resid$sd1, dsc_resid$sk1, dsc_resid$ku1))

}

par(mfrow=c(1,1))

# To launch ShinyStan (for the first fit):
launch_shinystan(fits[[2]])



