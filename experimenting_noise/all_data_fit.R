# -------------------------------------------------------
# This script:
# - Simulates data from a factor model with varying noise
# - Fits 1-factor and 3-factor Bayesian factor models using Stan
# - Calculates and plots diagnostics:
#     - Fraction of compatible covariance triplets (bootstrapped)
#     - Mean |correlation| and estimated loadings
#     - Residual correlations after model fit
#     - DSC distance (mean, SD, skewness, kurtosis) using Fisher z-transformed correlations
# - Compares how model fits and correlation structure change with noise
# -------------------------------------------------------


library(rstan)

# Set noise levels, data dims, and diagnostics storage
noise_levels <- seq(0.1, 2, by=0.4)
p <- ncol(all_data[[1]])
n <- nrow(all_data[[1]])

diagnostics <- data.frame(
   noise_sd = numeric(),
   mean_offdiag_cor = numeric(),
   mean_resid_cor_1f = numeric(),
   mean_resid_cor_3f = numeric(),
   compat_triplets = numeric()
)

# Prep to save residuals for DSC
all_resid_1f <- list()
all_resid_3f <- list()

# Compile/load Stan model
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
mod <- stan_model("mgps_factor_model.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/experimenting_noise")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

par(mfrow=c(3,2))  # 3 rows, 2 columns of plots per noise level

for (noise in noise_levels) {
   cat(sprintf("\n=== Noise SD: %.2f ===\n", noise))
   Y <- all_data[[as.character(noise)]]
   
   # --- BOOTSTRAP: Distribution of fraction of compatible triplets ---
   B <- 100
   compat_fraction_boot <- numeric(B)
   for (b in 1:B) {
      idx <- sample(1:n, replace=TRUE)
      Yb <- Y[idx, ]
      Cb <- cov(Yb)
      trip_compat <- sapply(combn(p, 3, simplify=FALSE), function(trip) {
         cpr <- Cb[trip[1], trip[2]]
         cps <- Cb[trip[1], trip[3]]
         crs <- Cb[trip[2], trip[3]]
         cpr*cps*crs >= 0
      })
      compat_fraction_boot[b] <- mean(trip_compat)
   }
   hist(compat_fraction_boot, breaks=30, main=sprintf("Noise=%.1f Boot Compat Triplets", noise),
        xlab="Fraction compatible (bootstrap)")
   
   # --- Observed off-diagonal correlations ---
   R_obs <- cor(Y)
   off_diag_corrs_obs <- R_obs[lower.tri(R_obs)]
   mean_offdiag_cor <- mean(abs(off_diag_corrs_obs))
   hist(off_diag_corrs_obs, breaks=30, main=sprintf("Noise=%.1f Obs Corr", noise), xlab="Correlation")
   
   # --- Method-of-moments lambda estimates (classic) ---
   C_obs <- cov(Y)
   lambda_1 <- sqrt(abs(C_obs[1,2])) # Choose a pair as reference
   lambda_mom <- c(lambda_1, sapply(2:p, function(j) C_obs[1, j] / lambda_1))
   hist(lambda_mom, breaks=30, main=sprintf("Noise=%.1f MoM Lambdas", noise), xlab="Lambda (MoM)")
   
   # --- Fraction of compatible triplets (for summary/diagnostics) ---
   trip_compat <- sapply(combn(p, 3, simplify=FALSE), function(trip) {
      cpr <- C_obs[trip[1], trip[2]]
      cps <- C_obs[trip[1], trip[3]]
      crs <- C_obs[trip[2], trip[3]]
      cpr*cps*crs >= 0
   })
   compat_triplets <- mean(trip_compat)
   
   # --- Fit 1-factor model ---
   cat("Fitting 1-factor model...\n")
   stan_data_1f <- list(N=n, P=p, K=1, Y=Y)
   fit_1f <- sampling(
      object       = mod,
      data         = stan_data_1f,
      chains       = 4,
      iter         = 4000,
      warmup       = 2000,
      seed         = 123 + round(100*noise),
      init         = "random",
      init_r       = 2,  # uniform(-2, 2)
      control      = list(
         adapt_delta = 0.99,
         max_treedepth = 15
      )
   )
   
   post_1f <- extract(fit_1f)
   Lambda_hat_1f <- matrix(apply(post_1f$Lambda, c(2,3), mean), ncol=1)
   eta_hat_1f <- Y %*% Lambda_hat_1f / sum(Lambda_hat_1f^2)
   Y_hat_1f <- eta_hat_1f %*% t(Lambda_hat_1f)
   resid_1f <- Y - Y_hat_1f
   all_resid_1f[[as.character(noise)]] <- resid_1f
   mean_resid_cor_1f <- mean(abs(cor(resid_1f)[lower.tri(cor(resid_1f))]))
   hist(cor(resid_1f)[lower.tri(cor(resid_1f))], breaks=30, main=sprintf("Noise=%.1f 1f Resid Corr", noise), xlab="Correlation")
   
   # --- Fit 3-factor model ---
   cat("Fitting 3-factor model...\n")
   stan_data_3f <- list(N=n, P=p, K=3, Y=Y)
   fit_3f <- sampling(
      object       = mod,
      data         = stan_data_3f,
      chains       = 4,
      iter         = 4000,
      warmup       = 2000,
      seed         = 456 + round(100*noise),
      init         = "random",
      init_r       = 2,  # uniform(-2, 2)
      control      = list(
         adapt_delta = 0.99,
         max_treedepth = 15
      )
   )
   post_3f <- extract(fit_3f)
   Lambda_hat_3f <- apply(post_3f$Lambda, c(2,3), mean)
   eta_hat_3f <- Y %*% Lambda_hat_3f %*% solve(t(Lambda_hat_3f) %*% Lambda_hat_3f)
   Y_hat_3f <- eta_hat_3f %*% t(Lambda_hat_3f)
   resid_3f <- Y - Y_hat_3f
   all_resid_3f[[as.character(noise)]] <- resid_3f
   mean_resid_cor_3f <- mean(abs(cor(resid_3f)[lower.tri(cor(resid_3f))]))
   hist(cor(resid_3f)[lower.tri(cor(resid_3f))], breaks=30, main=sprintf("Noise=%.1f 3f Resid Corr", noise), xlab="Correlation")
   
   # --- Fraction compatible triplets plot (as a bar) ---
   barplot(compat_triplets, names.arg=sprintf("Compat %.2f", compat_triplets), ylim=c(0,1),
           main=sprintf("Noise=%.1f Compat Triplets", noise))
   
   # --- Store diagnostics
   diagnostics <- rbind(
      diagnostics,
      data.frame(
         noise_sd=noise,
         mean_offdiag_cor=mean_offdiag_cor,
         mean_resid_cor_1f=mean_resid_cor_1f,
         mean_resid_cor_3f=mean_resid_cor_3f,
         compat_triplets=compat_triplets
      )
   )
}

par(mfrow=c(1,1))

# --- DSC function: mean, SD, skewness, kurtosis ---
dsc <- function(corrs, mu2 = 0, sd2 = NULL, sk2 = 0, ku2 = 3, n = NULL) {
   mu1 <- mean(corrs)
   sd1 <- sd(corrs)
   # Skewness and kurtosis
   if (sd1 < 1e-12) {
      sk1 <- 0
      ku1 <- 0
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
   return(list(DSC = D, mu1 = mu1, sd1 = sd1, sk1 = sk1, ku1 = ku1,
               mu2 = mu2, sd2 = sd2, sk2 = sk2, ku2 = ku2))
}

# === DSC Distance Calculation (always with Fisher-z correlations) ===

fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))

dsc_obs <- numeric(length(noise_levels))
dsc_resid_1f <- numeric(length(noise_levels))
dsc_resid_3f <- numeric(length(noise_levels))

for (i in seq_along(noise_levels)) {
   noise <- noise_levels[i]
   n <- nrow(all_data[[as.character(noise)]])
   
   # 1. Observed Fisher-z correlations
   R_obs <- cor(all_data[[as.character(noise)]])
   z_obs <- fisher_z(R_obs[lower.tri(R_obs)])
   dsc_obs[i] <- dsc(z_obs, n=n)$DSC
   
   # 2. Residual Fisher-z after 1-factor fit
   resid_1f <- all_resid_1f[[as.character(noise)]]
   R_resid_1f <- cor(resid_1f)
   z_resid_1f <- fisher_z(R_resid_1f[lower.tri(R_resid_1f)])
   dsc_resid_1f[i] <- dsc(z_resid_1f, n=n)$DSC
   
   # 3. Residual Fisher-z after 3-factor fit
   resid_3f <- all_resid_3f[[as.character(noise)]]
   R_resid_3f <- cor(resid_3f)
   z_resid_3f <- fisher_z(R_resid_3f[lower.tri(R_resid_3f)])
   dsc_resid_3f[i] <- dsc(z_resid_3f, n=n)$DSC
}

# Plot DSC distances (Fisher-z, all four moments)
plot(noise_levels, dsc_obs, type='b', col='black', pch=16,
     ylim=range(0, dsc_obs, dsc_resid_1f, dsc_resid_3f, na.rm=TRUE),
     xlab="Noise SD", ylab="DSC distance (Fisher-z, all moments)",
     main="DSC (mean, sd, skew, kurtosis) vs. Noise")
lines(noise_levels, dsc_resid_1f, type='b', col='red', pch=16)
lines(noise_levels, dsc_resid_3f, type='b', col='blue', pch=16)
legend("bottomleft", legend=c("Observed", "1-factor residual", "3-factor residual"),
       col=c("black", "red", "blue"), lty=1, pch=16)



