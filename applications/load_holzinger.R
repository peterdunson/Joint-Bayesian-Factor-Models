# load_holzinger.R

#Holzinger–Swineford 1939145 students × 9 cognitive-test scores

# ─── 1) Load & Prepare Holzinger–Swineford Data ───────────────────────────
library(lavaan)

data("HolzingerSwineford1939", package = "lavaan")
df_hs <- subset(
   HolzingerSwineford1939,
   school == "Pasteur",
   select = paste0("x", 1:9)
)

# Standardize (center + scale)
Y <- scale(df_hs, center = TRUE, scale = TRUE)
n <- nrow(Y)
p <- ncol(Y)
K <- 3   # choose number of latent factors

# ─── 2) Compile MGPS Stan Model ────────────────────────────────────────────
library(rstan)
library(bayesplot)

set.seed(21)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
mod <- stan_model("mgps_factor_model.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/applications")

# ─── 3) Fit with Random Initial Values ────────────────────────────────────
stan_data <- list(N = n, P = p, K = K, Y = Y)

fit_hs <- sampling(
   object       = mod,
   data         = stan_data,
   chains       = 4,
   iter         = 6000,
   warmup       = 3000,
   seed         = 21,
   init         = "random",
   init_r       = 2,
   control      = list(adapt_delta = 0.99, max_treedepth = 15)
)

# ─── 4) Extract Posterior & Compute Posterior-Mean Loadings ───────────────
post_hs     <- extract(fit_hs)
Lambda_hat  <- apply(post_hs$Lambda, c(2,3), mean)

# ─── 5) Save Results ──────────────────────────────────────────────────────
saveRDS(
   list(
      fit         = fit_hs,
      posterior   = post_hs,
      Lambda_hat  = Lambda_hat
   ),
   file = "fit_Holzinger_mgps.rds"
)

# ─── 6) Diagnostics ────────────────────────────────────────────────────────
sum_stats   <- summary(fit_hs)$summary
cat("max R̂    =", max(sum_stats[,"Rhat"], na.rm=TRUE), "\n")
cat("min n_eff =", min(sum_stats[,"n_eff"], na.rm=TRUE), "\n")
cat("min BFMI  =", min(rstan::get_bfmi(fit_hs), na.rm=TRUE), "\n")

# ─── 7) (Optional) ShinyStan ───────────────────────────────────────────────
# library(shinystan)
# launch_shinystan(fit_hs)














# ---------------------------------------------------------
# DSC for Holzinger–Swineford (Pasteur, 9 variables, 3 factors)
# Empirical Null for Comparison
# ---------------------------------------------------------

# -- 1) Load fit and data --
fit_path <- "fit_Holzinger_mgps.rds"

fit_obj <- readRDS(fit_path)
Lambda_hat <- fit_obj$Lambda_hat

library(lavaan)
data("HolzingerSwineford1939", package = "lavaan")
df_hs <- subset(HolzingerSwineford1939, school == "Pasteur", select = paste0("x", 1:9))
Y <- scale(df_hs, center = TRUE, scale = TRUE)
n <- nrow(Y)
p <- ncol(Y)

# -- 2) Compute factor scores, recon, residuals --
eta_hat <- Y %*% Lambda_hat %*% solve(t(Lambda_hat) %*% Lambda_hat)
Y_hat <- eta_hat %*% t(Lambda_hat)
resid_mat <- Y - Y_hat

# -- 3) Compute Fisher-z off-diagonal correlations --
fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
R_obs <- cor(Y)
R_resid <- cor(resid_mat)
z_obs <- fisher_z(R_obs[lower.tri(R_obs)])
z_resid <- fisher_z(R_resid[lower.tri(R_resid)])

# -- 4) DSC function --
dsc <- function(corrs, mu2 = 0, sd2 = NULL, sk2 = 0, ku2 = 3, n = NULL) {
   mu1 <- mean(corrs)
   sd1 <- sd(corrs)
   if (sd1 < 1e-12) { sk1 <- 0; ku1 <- 0 }
   else {
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

# -- 5) Compute DSC for residuals --
dsc_resid <- dsc(z_resid, n = n)

cat(sprintf(
   "Residual DSC: %.3f (mean=%.3f, sd=%.3f, skew=%.3f, kurtosis=%.3f)\n",
   dsc_resid$DSC, dsc_resid$mu1, dsc_resid$sd1, dsc_resid$sk1, dsc_resid$ku1))

# -- 6) Empirical Null Distribution of DSC (noise only) --
set.seed(999)
B <- 1000
dsc_null <- numeric(B)
for (b in 1:B) {
   Y_null <- matrix(rnorm(n * p), nrow = n, ncol = p)
   R_null <- cor(Y_null)
   z_null <- fisher_z(R_null[lower.tri(R_null)])
   dsc_null[b] <- dsc(z_null, n = n)$DSC
}

# -- 7) Plot null + observed --
hist(dsc_null, breaks = 40, main = "Empirical Null Distribution of DSC (Holzinger)", 
     xlab = "DSC (noise-only)", col = "gray", border = "white")
abline(v = dsc_resid$DSC, col = "blue", lwd = 3)
abline(v = quantile(dsc_null, c(0.025, 0.975)), col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Observed residual DSC", "95% null envelope"),
       col = c("blue", "red"), lwd = c(3,2), lty = c(1,2))

cat(sprintf("Empirical null DSC: Mean = %.3f, SD = %.3f\n", mean(dsc_null), sd(dsc_null)))
cat(sprintf("95%% null interval: %.3f to %.3f\n", 
            quantile(dsc_null, 0.025), quantile(dsc_null, 0.975)))
