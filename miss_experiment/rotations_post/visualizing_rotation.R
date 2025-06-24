# visualize_rotated_loadings.R

# ─── 0) Libraries ─────────────────────────────────────────────────────────
library(pheatmap)    # install.packages("pheatmap")
# If you prefer ggplot2 instead:
# library(ggplot2)
# library(reshape2)

# ─── 1) Load your varimax‐rotated summaries ────────────────────────────────
varimax_res <- readRDS("fit_Joint_scen2_scale_all_rotated_fullposterior.rds")
Lambda_vm_hat <- varimax_res$Lambda_vm_hat   # p × K matrix

# ─── 2) Load your Procrustes‐rotated summaries ────────────────────────────
procrustes_res <- readRDS("fit_Joint_scen2_scale_all_procrustes_rotated.rds")
Lambda_pr_hat <- procrustes_res$Lambda_pr_hat  # p × K matrix

# ─── 3A) Simple heatmaps with pheatmap ────────────────────────────────────
cols <- colorRampPalette(c("grey5","grey43","white"))(50)

pheatmap(
   Lambda_vm_hat,
   color        = cols,
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   main         = "Varimax‐Rotated Loadings"
)

pheatmap(
   Lambda_pr_hat,
   color        = cols,
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   main         = "Procrustes‐Rotated Loadings"
)


