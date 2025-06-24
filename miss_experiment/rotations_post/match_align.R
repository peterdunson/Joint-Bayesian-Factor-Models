# matchalign_thresholded_unrotated_plot.R

# 1) Libraries
library(infinitefactor)  # install.packages("infinitefactor") if needed
library(pheatmap)        # install.packages("pheatmap")

# 2) Load true & unrotated estimate
sim         <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Lambda_true <- sim$Lambda

fit_unrot   <- dat
Lambda_hat  <- fit_unrot$Lambda_hat

# 3) Threshold at 0.05
Lambda_hat_thr <- Lambda_hat
Lambda_hat_thr[abs(Lambda_hat_thr) < 0.05] <- 0

# 4) Match-align the thresholded loadings to truth
Lambda_hat_ma_thr <- msf(Lambda_hat_thr, Lambda_true)

# 5) Plot the match-aligned, thresholded loadings
cols <- colorRampPalette(c("white","grey43","grey5"))(50)

pheatmap(
   Lambda_hat_ma_thr,
   color        = cols,
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   main         = "Unrotated Loadings (Thresholded & Match-Aligned)"
)
