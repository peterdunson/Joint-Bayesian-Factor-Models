# align_and_plot_ref_basis.R

# ─── 1) Libraries ─────────────────────────────────────────────────────────
library(GPArotation)   # install.packages("GPArotation")
library(infinitefactor) # install.packages("infinitefactor")
library(pheatmap)       # install.packages("pheatmap")

# ─── 2) Load true & TEB-FAR fit ────────────────────────────────────────────
sim           <- readRDS("~/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Lambda_true   <- as.matrix(sim$Lambda)          # (p+1) × K

#fit_TEBFAR    <- readRDS("~/Desktop/fit_TEBFAR_scen2.rds")
lambda_list   <- fit_TEBFAR$lambdaSamps         # list of posterior draws

# ─── 3) Compute posterior–median unrotated MGPS loadings ──────────────────
p <- nrow(Lambda_true); K <- ncol(Lambda_true)
Lambda_hat <- matrix(NA_real_, p, K)
for (j in 1:p) for (k in 1:K) {
   vals <- vapply(lambda_list,
                  function(mat) if (ncol(mat) >= k) mat[j,k] else NA_real_,
                  FUN.VALUE = 0.0)
   Lambda_hat[j,k] <- median(vals, na.rm = TRUE)
}

# ─── 4) Match-align raw to truth in factor-space ──────────────────────────
Lambda_ma <- msf(Lambda_hat, Lambda_true)

# ─── 5) Get the varimax rotation R from the true loadings ─────────────────
vm_ref <- GPArotation::Varimax(Lambda_true)
R_ref  <- vm_ref$Th   # K × K orthogonal matrix

# ─── 6) Rotate all three matrices into that same basis ────────────────────
Lambda_true_rot <- Lambda_true %*% R_ref
Lambda_hat_rot  <- Lambda_hat    %*% R_ref
Lambda_ma_rot   <- Lambda_ma     %*% R_ref

# ─── 7) Common diverging color scale ─────────────────────────────────────
all_mats <- list(
   True     = Lambda_true_rot,
   Match    = Lambda_ma_rot,
   Raw      = Lambda_hat_rot
)
maxval <- max(abs(unlist(all_mats)))
cols   <- colorRampPalette(c("pink","white","grey4"))(100)
brks   <- seq(-maxval, maxval, length.out = length(cols) + 1)

# ─── 8) Plot them side-by-side ────────────────────────────────────────────
for (nm in names(all_mats)) {
   pheatmap(
      all_mats[[nm]],
      color        = cols,
      breaks       = brks,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      main         = nm
   )
}


