# evaluate_mgps_from_df_all.R

# ─── 1) Libraries ─────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(mvtnorm)   # for cor()

# ─── 2) Load simulation truth & df_all ────────────────────────────────────
sim         <- readRDS("~/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Lambda_true <- sim$Lambda    # (p+1) × K
Omega_true  <- sim$Omega
# assume df_all is already in your workspace, produced by run_mgps_chains_and_trace.R

# ─── 3) Build posterior-median Λ_est from df_all ───────────────────────────
# df_all has columns: Variable, Factor, Iteration, Loading, Chain
Lambda_est_df <- df_all %>%
   group_by(Variable, Factor) %>%
   summarize(median_loading = median(Loading), .groups = "drop")

# convert back to matrix in same order as Lambda_true
p <- nrow(Lambda_true); K <- ncol(Lambda_true)
Lambda_est <- matrix(
   Lambda_est_df$median_loading,
   nrow = p, ncol = K,
   byrow = FALSE,
   dimnames = list(
      rownames(Lambda_true),
      colnames(Lambda_true)
   )
)

# ─── 4) MSE between Λ_est and Λ_true ───────────────────────────────────────
mse_lambda <- mean((Lambda_est - Lambda_true)^2)
cat("MSE(Lambda_est, Lambda_true) =", round(mse_lambda, 5), "\n")

# ─── 5) Implied covariance Ω_est & R² ─────────────────────────────────────
Sigma_true <- diag(rep(0.2, p))           # same noise as in sim
Omega_est  <- Lambda_est %*% t(Lambda_est) + Sigma_true
R2_cov     <- cor(as.vector(Omega_est),
                  as.vector(Omega_true))^2
cat("R²(Omega_est, Omega_true)    =", round(R2_cov, 5), "\n")

# ─── 6) Return results (invisible) ───────────────────────────────────────
invisible(list(
   Lambda_MSE = mse_lambda,
   Omega_R2   = R2_cov
))

