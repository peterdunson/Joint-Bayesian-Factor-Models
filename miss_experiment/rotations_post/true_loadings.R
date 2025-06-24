# align_and_plot_loadings.R

# ─── 1) Libraries ─────────────────────────────────────────────────────────
library(pheatmap)    # install.packages("pheatmap") if needed

# ─── 2) Load true simulated loadings ───────────────────────────────────────
sim <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Lambda_true <- sim$Lambda    # (p+1) × K

# ─── 3) Load estimated loadings ────────────────────────────────────────────
# 3a) Unrotated posterior‐mean
fit_unrot <- dat
Lambda_hat <- fit_unrot$Lambda_hat    # (p+1) × K

# 3b) Varimax‐rotated posterior‐mean
varimax_res    <- readRDS("fit_Joint_scen2_scale_all_rotated_fullposterior.rds")
Lambda_vm_hat  <- varimax_res$Lambda_vm_hat

# 3c) Procrustes‐rotated posterior‐mean
pr_res         <- readRDS("fit_Joint_scen2_scale_all_procrustes_rotated.rds")
Lambda_pr_hat  <- pr_res$Lambda_pr_hat

# ─── 4) Compute sign‐correction vectors & align ─────────────────────────────
# Function to align columns of est to true by sign of column‐wise dot product
align_sign <- function(L_true, L_est) {
   sgn <- sign(colSums(L_true * L_est))
   # if any column sums are zero, leave sign = +1
   sgn[sgn == 0] <- 1
   sweep(L_est, 2, sgn, `*`)
}

Lambda_hat_aligned  <- align_sign(Lambda_true, Lambda_hat)
Lambda_vm_aligned   <- align_sign(Lambda_true, Lambda_vm_hat)
Lambda_pr_aligned   <- align_sign(Lambda_true, Lambda_pr_hat)

# ─── 5) Plot heatmaps ──────────────────────────────────────────────────────
cols <- colorRampPalette(c("white","grey43","grey5"))(50)

# True
pheatmap(
   Lambda_true,
   color        = cols,
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   main         = "True Loadings"
)

# Unrotated (sign‐aligned)
pheatmap(
   Lambda_hat_aligned,
   color        = cols,
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   main         = "Unrotated Posterior-Mean Loadings\n(sign‐aligned)"
)

# Varimax (sign‐aligned)
pheatmap(
   Lambda_vm_aligned,
   color        = cols,
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   main         = "Varimax-Rotated Loadings\n(sign‐aligned)"
)

# Procrustes (sign‐aligned)
pheatmap(
   Lambda_pr_aligned,
   color        = cols,
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   main         = "Procrustes-Rotated Loadings\n(sign‐aligned)"
)
