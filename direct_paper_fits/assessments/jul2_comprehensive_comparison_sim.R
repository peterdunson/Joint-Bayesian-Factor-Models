# Full evaluation script: load simulation, fits, compute estimators, compare to truth,
# plot heatmap, and compare λ₁·bᵢ & residual variances

# ---- CHOOSE SCENARIO ----
scenario <- 2
sim_path <- sprintf(
   "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000_1.rds",
   scenario
)
sim <- readRDS(sim_path)

# ---- CENTER & SCALE DATA ----
# Full data (including y)
Y_full <- scale(sim$Y, center = TRUE, scale = TRUE)   # (n × (p+1))
# Predictor-only data
Y_pred <- Y_full[, -1]                                  # (n × p)
n <- nrow(Y_full)
p <- ncol(Y_pred)

# ---- Load fits ----
fit_dir  <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/direct_paper_fits/storing_fit"
fit_MGSP <- readRDS(file.path(fit_dir, "MGSP",      "fit_joint_scen2_k1.rds"))
fit_HS   <- readRDS(file.path(fit_dir, "Horseshoe", "fit_horseshoe_scen2_k1.rds"))
fit_SSL  <- readRDS(file.path(fit_dir, "spike_slab","fit_spikeslab_scen2_k1.rds"))

# ---- MoM Trios estimator for lambda_1 (predictors only) ----
C <- cor(Y_pred)
lambda2_est_trios <- numeric(p)
for (pp in 1:p) {
   vals <- c()
   for (qq in setdiff(1:p, pp)) {
      for (rr in setdiff(1:p, c(pp, qq))) {
         c_pq <- C[pp, qq]; c_pr <- C[pp, rr]; c_qr <- C[qq, rr]
         if (!is.na(c_qr) && abs(c_qr) > 1e-8) {
            vals <- c(vals, (c_pq * c_pr) / c_qr)
         }
      }
   }
   lambda2_est_trios[pp] <- mean(vals, na.rm = TRUE)
}
lambda_est_mom_trios <- sign(lambda2_est_trios) * sqrt(abs(lambda2_est_trios))
b_hat_mom_trios     <- as.numeric(Y_pred %*% lambda_est_mom_trios) / sum(lambda_est_mom_trios^2)

# ---- MoM Eigen (PC) estimator for lambda_1 (predictors only) ----
S                 <- cov(Y_pred)
e                 <- eigen(S)
d                 <- e$values
V                 <- e$vectors
sigma2_hat        <- mean(d[-1])
lambda_mag        <- sqrt(max(d[1] - sigma2_hat, 0))
lambda1_hat_eigen <- lambda_mag * V[, 1]
b_hat_mom_eigen   <- as.numeric(Y_pred %*% lambda1_hat_eigen) / sum(lambda1_hat_eigen^2)

# ---- Extract raw Bayesian loadings & scores (full data) ----
extract_lambda_b <- function(fit_obj, Y) {
   λ    <- as.numeric(fit_obj$Lambda_hat[,1])      # length p+1
   b    <- as.numeric(Y %*% λ) / sum(λ^2)
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
   MoM_trios = c(NA, lambda_est_mom_trios),       # pad with NA for y-row
   MoM_eigen = c(NA, lambda1_hat_eigen),
   MGSP      = lambda_MGSP_norm,
   HS        = lambda_HS_norm,
   SSL       = lambda_SSL_norm
)

b_hat_list <- list(
   True      = eta_true,
   MoM_trios = b_hat_mom_trios,
   MoM_eigen = b_hat_mom_eigen,
   MGSP      = b_MGSP_norm,
   HS        = b_HS_norm,
   SSL       = b_SSL_norm
)

# ---- Print first-factor loadings & first 20 scores ----
cat("=== First-Factor Loadings (λ1) ===\n")
print(as.data.frame(lambda_est_list))

cat("\n=== First-Factor Scores (b_i) [first 20] ===\n")
print(as.data.frame(b_hat_list)[1:20, ])

# ---- Heatmap of True vs. Estimated Loadings (no row names) ----
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

# ---- Compare λ₁·bᵢ & residual variances across methods ----
lam_list <- list(
   MoM_trios = lambda_est_mom_trios,
   MoM_eigen = lambda1_hat_eigen,
   MGSP      = lambda_MGSP_norm[-1],  # drop y-row
   HS        = lambda_HS_norm[-1],
   SSL       = lambda_SSL_norm[-1]
)
b_list <- list(
   MoM_trios = b_hat_mom_trios,
   MoM_eigen = b_hat_mom_eigen,
   MGSP      = b_MGSP_norm,
   HS        = b_HS_norm,
   SSL       = b_SSL_norm
)

proj_mean_mat <- matrix(NA, nrow = p, ncol = length(lam_list),
                        dimnames = list(paste0("Var",1:p), names(lam_list)))
resid_var_mat <- proj_mean_mat

for (m in names(lam_list)) {
   λ   <- lam_list[[m]]
   b   <- b_list[[m]]
   Ym  <- outer(b, λ)
   R   <- Y_full[, -1] - Ym    # residuals on predictors only
   proj_mean_mat[, m] <- colMeans(Ym)
   resid_var_mat[, m] <- apply(R, 2, var)
}

cat("=== Mean(λ₁·bᵢ) by variable & method ===\n")
print(round(proj_mean_mat, 4))

cat("\n=== Residual variances by variable & method ===\n")
print(round(resid_var_mat, 4))

pheatmap(proj_mean_mat,
         main = "Mean Contribution: λ₁·bᵢ",
         cluster_rows  = FALSE,
         cluster_cols  = FALSE,
         show_rownames = FALSE,
         border_color  = NA)

pheatmap(resid_var_mat,
         main = "Residual Variance",
         cluster_rows  = FALSE,
         cluster_cols  = FALSE,
         show_rownames = FALSE,
         border_color  = NA)


# 1. Compute proportion explained
explained_mat <- 1 - resid_var_mat
# (each entry is fraction of Var(Y_j) captured by the first factor)

# 2. Plot heat map
library(pheatmap)
pheatmap(
   explained_mat,
   main          = "Proportion of Variance Explained\nby First Factor (per Variable & Method)",
   cluster_rows  = FALSE,
   cluster_cols  = FALSE,
   show_rownames = FALSE,
   border_color  = NA,
   fontsize_col  = 10
)




# 1. Pull out the raw simulated predictors and true factor
Y_raw    <- sim$Y[, -1]         # n × p, unscaled
lambda   <- sim$Lambda[-1, 1]   # true loadings on raw scale
eta_true <- sim$eta[, 1]        # true scores

# 2. Form the raw residuals
E_raw <- Y_raw - outer(eta_true, lambda)

# 3. Compute each variable’s variance
var_raw <- apply(E_raw, 2, var)

# 4. Print and compare to 0.2
print(round(var_raw, 3))



















# 1. Raw predictors
Y_raw   <- sim$Y[, -1]          # n × p
orig_sds <- apply(Y_raw, 2, sd) # length p

# 2. Align MGSP loadings by dropping the y‐row
lambda_mgsp_scaled <- out_MGSP$lambda1       # length p+1
lambda_mgsp_pred   <- lambda_mgsp_scaled[-1] # drop the first

# 3. Un‐scale
lambda_raw_mgsp <- lambda_mgsp_pred / orig_sds

# 4. Now you can recompute raw scores and residual variances
b_raw_mgsp <- as.numeric(Y_raw %*% lambda_raw_mgsp) / sum(lambda_raw_mgsp^2)

get_var <- function(Y, b, lam) {
   R <- Y - outer(b, lam)
   apply(R, 2, var)
}

resid_var_raw_mgsp <- get_var(Y_raw, b_raw_mgsp, lambda_raw_mgsp)
print(round(resid_var_raw_mgsp, 3))


