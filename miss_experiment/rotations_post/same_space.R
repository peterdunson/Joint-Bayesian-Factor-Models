# 1) Clean up the true‐loadings matrix
Lambda_true_mat <- matrix(
   as.numeric(unlist(sim$Lambda)),
   nrow = nrow(sim$Lambda),
   byrow = FALSE
)
storage.mode(Lambda_true_mat) <- "double"

# 2) Run GPArotation’s Varimax
library(GPArotation)
vm_ref <- GPArotation::Varimax(Lambda_true_mat)

# 3) Extract the K×K orthogonal rotation
R_ref <- vm_ref$Th
dim(R_ref)  # should be K × K


# align_and_matchall_loadings.R
# align_raw_then_rotate_all.R

# ─── 1) Libraries ─────────────────────────────────────────────────────────
library(GPArotation)    # install.packages("GPArotation")
library(infinitefactor)  # install.packages("infinitefactor")
library(pheatmap)       # install.packages("pheatmap")

# ─── 2) Load true & estimated loadings ───────────────────────────────────
sim             <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Lambda_true     <- as.matrix(sim$Lambda)       # (p+1) x K
storage.mode(Lambda_true) <- "double"

fit_unrot       <- readRDS("fit_Joint_scen2_scale_all_6k.rds")
Lambda_hat      <- as.matrix(fit_unrot$Lambda_hat)

varimax_res     <- readRDS("fit_Joint_scen2_scale_all_rotated_fullposterior_6k.rds")
Lambda_vm_hat   <- as.matrix(varimax_res$Lambda_vm_hat)

pr_res          <- readRDS("fit_Joint_scen2_scale_all_procrustes_rotated_6k.rds")
Lambda_pr_hat   <- as.matrix(pr_res$Lambda_pr_hat)

# ─── 3) Match‐align only the raw (unrotated) estimate to the true loadings ─
Lambda_unrot_ma <- msf(Lambda_hat, Lambda_true)

# ─── 4) Compute reference varimax rotation on the true loadings ──────────
vm_ref <- GPArotation::Varimax(Lambda_true)
R_ref  <- vm_ref$Th    # K × K rotation matrix

# ─── 5) Rotate everything into the true’s varimax basis ──────────────────
Lambda_true_rot     <- Lambda_true      %*% R_ref
Lambda_unrot_ma_rot <- Lambda_unrot_ma  %*% R_ref
Lambda_vm_rot       <- Lambda_vm_hat    %*% R_ref
Lambda_pr_rot       <- Lambda_pr_hat    %*% R_ref

# ─── 6) Plot all four in the common reference frame ───────────────────────
cols <- colorRampPalette(c("white","grey43","grey5"))(50)

pheatmap(
   Lambda_true_rot, 
   color        = cols, 
   cluster_rows = FALSE, 
   cluster_cols = FALSE,
   main         = "True Loadings (Ref‐Varimax)"
)

pheatmap(
   Lambda_unrot_ma_rot, 
   color        = cols, 
   cluster_rows = FALSE, 
   cluster_cols = FALSE,
   main         = "Raw Posterior‐Mean → MatchAln → Ref‐Varimax"
)

pheatmap(
   Lambda_vm_rot, 
   color        = cols, 
   cluster_rows = FALSE, 
   cluster_cols = FALSE,
   main         = "Varimax Estimate → Ref‐Varimax"
)

pheatmap(
   Lambda_pr_rot, 
   color        = cols, 
   cluster_rows = FALSE, 
   cluster_cols = FALSE,
   main         = "Procrustes Estimate → Ref‐Varimax"
)
