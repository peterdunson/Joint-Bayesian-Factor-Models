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
cols <- colorRampPalette(c("darkgreen","white","darkred"))(50)

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

# ─── 3B) (Optional) If you’d rather use ggplot2 ────────────────────────────
# # helper to ggplot‐heat
# plot_heatmap <- function(L, title) {
#   mat <- as.matrix(L)
#   df  <- melt(mat)
#   names(df) <- c("Variable","Factor","Loading")
#   ggplot(df, aes(x = Factor, y = Variable, fill = Loading)) +
#     geom_tile() +
#     scale_fill_gradient2(low="darkgreen", mid="white", high="darkred", midpoint=0) +
#     theme_minimal() +
#     labs(title = title, x = "Factor", y = "Variable")
# }
#
# plot_heatmap(Lambda_vm_hat, "Varimax‐Rotated Loadings")
# plot_heatmap(Lambda_pr_hat, "Procrustes‐Rotated Loadings")
