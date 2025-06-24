# compare_all_five_heatmaps.R

library(GPArotation)
library(infinitefactor)
library(pheatmap)

# 1) Load and prepare matrices
sim               <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Lambda_true       <- as.matrix(sim$Lambda)

fit_unrot         <- readRDS("fit_Joint_scen2_scale_all_6k.rds")
Lambda_hat        <- as.matrix(fit_unrot$Lambda_hat)

varimax_res       <- readRDS("fit_Joint_scen2_scale_all_rotated_fullposterior_6k.rds")
Lambda_vm_hat     <- as.matrix(varimax_res$Lambda_vm_hat)

pr_res            <- readRDS("fit_Joint_scen2_scale_all_procrustes_rotated_6k.rds")
Lambda_pr_hat     <- as.matrix(pr_res$Lambda_pr_hat)

# 2) Match‐align raw output to truth, then rotation
Lambda_unrot_ma   <- msf(Lambda_hat, Lambda_true)

# 3) Reference varimax of true loadings
vm_ref     <- GPArotation::Varimax(Lambda_true)
R_ref      <- vm_ref$Th

# 4) Rotate each into true’s varimax basis
Lambda_true_rot    <- Lambda_true       %*% R_ref
Lambda_unrot_rot   <- Lambda_hat        %*% R_ref
Lambda_unrot_ma_rot<- Lambda_unrot_ma   %*% R_ref
Lambda_vm_rot      <- Lambda_vm_hat     %*% R_ref
Lambda_pr_rot      <- Lambda_pr_hat     %*% R_ref

# 5) Common color scale
all_mats <- list(
   True                  = Lambda_true_rot,
   Raw_matchaligned_rot  = Lambda_unrot_ma_rot,
   Raw_rotated           = Lambda_unrot_rot,
   Varimax_rotated       = Lambda_vm_rot,
   Procrustes_rotated    = Lambda_pr_rot
)
maxval <- max(abs(unlist(all_mats)))
cols   <- colorRampPalette(c("pink","white","grey12"))(100)
brks   <- seq(-maxval, maxval, length.out = length(cols) + 1)

# 6) Plot
for(name in names(all_mats)) {
   pheatmap(
      all_mats[[name]],
      color        = cols,
      breaks       = brks,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      main         = name
   )
}

