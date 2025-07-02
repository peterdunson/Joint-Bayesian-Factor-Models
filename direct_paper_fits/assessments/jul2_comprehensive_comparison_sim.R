# Full evaluation script: load simulation, fits, compute estimators on all variables, compare to truth, and plot

# ---- CHOOSE SCENARIO ----
scenario <- 2
sim_path <- sprintf(
   "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000_1.rds",
   scenario
)
sim <- readRDS(sim_path)

# ---- CENTER & SCALE FULL DATA (y + predictors) ----
Y_full <- scale(sim$Y, center = TRUE, scale = TRUE)   # (n × (p+1))
n      <- nrow(Y_full)
p_full <- ncol(Y_full)  # includes y variable

# ---- Load Bayesian fits ----
fit_dir  <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/direct_paper_fits/storing_fit"
fit_MGSP <- readRDS(file.path(fit_dir, "MGSP",      "fit_joint_scen2_k1.rds"))
fit_HS   <- readRDS(file.path(fit_dir, "Horseshoe", "fit_horseshoe_scen2_k1.rds"))
fit_SSL  <- readRDS(file.path(fit_dir, "spike_slab","fit_spikeslab_scen2_k1.rds"))

# ---- MoM Trios estimator for lambda_1 on all variables ----
C_full <- cor(Y_full)
lambda2_trios <- numeric(p_full)
for (pp in 1:p_full) {
   vals <- c()
   others <- setdiff(1:p_full, pp)
   for (qq in others) {
      for (rr in setdiff(others, qq)) {
         cpq <- C_full[pp, qq]; cpr <- C_full[pp, rr]; cqr <- C_full[qq, rr]
         if (!is.na(cqr) && abs(cqr) > 1e-8) {
            vals <- c(vals, (cpq * cpr) / cqr)
         }
      }
   }
   lambda2_trios[pp] <- mean(vals, na.rm = TRUE)
}
lambda_est_trios_all <- sign(lambda2_trios) * sqrt(abs(lambda2_trios))
b_hat_trios_all     <- as.numeric(Y_full %*% lambda_est_trios_all) / sum(lambda_est_trios_all^2)

# ---- MoM Eigen (PC) estimator for lambda_1 on all variables ----
S_full            <- cov(Y_full)
e_full            <- eigen(S_full)
d_full            <- e_full$values
V_full            <- e_full$vectors
sigma2_hat_full   <- mean(d_full[-1])
lambda_mag_full   <- sqrt(max(d_full[1] - sigma2_hat_full, 0))
lambda1_eigen_all <- lambda_mag_full * V_full[, 1]
b_hat_eigen_all   <- as.numeric(Y_full %*% lambda1_eigen_all) / sum(lambda1_eigen_all^2)

# ---- Extract raw Bayesian loadings & scores ----
extract_lambda_b <- function(fit_obj, Y) {
   λ <- as.numeric(fit_obj$Lambda_hat[, 1])
   b <- as.numeric(Y %*% λ) / sum(λ^2)
   list(lambda1 = λ, b_hat = b)
}
out_MGSP <- extract_lambda_b(fit_MGSP, Y_full)
out_HS   <- extract_lambda_b(fit_HS,   Y_full)
out_SSL  <- extract_lambda_b(fit_SSL,  Y_full)

# ---- Normalize Bayesian loadings & recompute scores ----
normalize <- function(v) v / sqrt(sum(v^2))
lambda_MGSP_norm <- normalize(out_MGSP$lambda1)
lambda_HS_norm   <- normalize(out_HS$lambda1)
lambda_SSL_norm  <- normalize(out_SSL$lambda1)

b_MGSP_norm <- as.numeric(Y_full %*% lambda_MGSP_norm)
b_HS_norm   <- as.numeric(Y_full %*% lambda_HS_norm)
b_SSL_norm  <- as.numeric(Y_full %*% lambda_SSL_norm)

# ---- Truth from simulation ----
Lambda_true <- sim$Lambda[, 1]  # length p+1
eta_true    <- sim$eta[, 1]     # length n

# ---- Gather into lists ----
lambda_est_list <- list(
   True      = Lambda_true,
   MoM_trios = lambda_est_trios_all,
   MoM_eigen = lambda1_eigen_all,
   MGSP      = lambda_MGSP_norm,
   HS        = lambda_HS_norm,
   SSL       = lambda_SSL_norm
)

b_hat_list <- list(
   True      = eta_true,
   MoM_trios = b_hat_trios_all,
   MoM_eigen = b_hat_eigen_all,
   MGSP      = b_MGSP_norm,
   HS        = b_HS_norm,
   SSL       = b_SSL_norm
)

# ---- Print loadings & scores ----
cat("=== First-Factor Loadings (λ1) ===\n")
print(as.data.frame(lambda_est_list))
cat("\n=== First-Factor Scores (b_i) [first 20] ===\n")
print(as.data.frame(b_hat_list)[1:20, ])

# ---- Heatmap of loadings ----
library(pheatmap)
loadings_df <- as.data.frame(lambda_est_list)
rownames(loadings_df) <- NULL
pheatmap(
   loadings_df,
   main          = "True vs. Estimated First-Factor Loadings",
   cluster_rows  = FALSE,
   cluster_cols  = FALSE,
   show_rownames = FALSE,
   border_color  = NA,
   fontsize_col  = 10
)

# ---- Density of residual correlations for each method ----
get_offdiag <- function(M) M[lower.tri(M)]
resid_cor <- function(b, λ) {
   R <- Y_full - outer(b, λ)
   get_offdiag(cor(R))
}
r_trio  <- resid_cor(b_hat_trios_all, lambda_est_trios_all)
r_eigen <- resid_cor(b_hat_eigen_all, lambda1_eigen_all)
r_mgsp  <- resid_cor(b_MGSP_norm,      lambda_MGSP_norm)




# ---- Density overlay of residual off-diagonal correlations (including observed) ----
get_offdiag <- function(M) M[lower.tri(M)]

# observed residuals = just the raw correlations
r_obs <- get_offdiag(cor(Y_full))

# method residuals
r_trio  <- resid_cor(b_hat_trios_all,  lambda_est_trios_all)
r_eigen <- resid_cor(b_hat_eigen_all,  lambda1_eigen_all)
r_mgsp  <- resid_cor(b_MGSP_norm,       lambda_MGSP_norm)

# plot
plot(density(r_obs),  lwd = 2, col = "grey50",
     main = "Residual Off-Diag Correlations (k=1)",
     xlab = "Correlation", xlim = c(-1,1))
lines(density(r_trio),  lwd = 2, col = "black")
lines(density(r_eigen), lwd = 2, col = "blue")
lines(density(r_mgsp),  lwd = 2, col = "darkgreen")
legend("topright",
       legend = c("Observed","MoM-Trio","MoM-Eigen","MGSP"),
       col    = c("grey50","black","blue","darkgreen"),
       lwd    = 2)

