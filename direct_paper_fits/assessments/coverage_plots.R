# ----------------------------------------------
# Coverage Heatmaps for Each Model (Sim2)
# "95% CI Coverage & Interval Endpoints\nlower (true) upper"
# ----------------------------------------------

library(pheatmap)

# 1. Get truth
Lambda_true <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")$Lambda

fit_list <- list(
   MGSP   = fit_MGSP,
   HS     = fit_HS,
   SSL    = fit_SSL,
   MASS   = fit_MASS,
   Robust = fit_RobustMGSP
)

for (mod in names(fit_list)) {
   cat("\n======", mod, "======\n")
   obj  <- fit_list[[mod]]
   post <- obj$posterior
   # Posterior draws: [iter, j, k]
   Lam_draws <- post$Lambda
   p1 <- nrow(Lambda_true)
   K  <- ncol(Lambda_true)
   
   # Preallocate
   cover_mat <- matrix(FALSE,  nrow = p1, ncol = K)
   lower_mat <- matrix(NA_real_, nrow = p1, ncol = K)
   upper_mat <- matrix(NA_real_, nrow = p1, ncol = K)
   
   for (j in seq_len(p1)) {
      for (k in seq_len(K)) {
         draws <- Lam_draws[, j, k]
         ci    <- quantile(draws, c(0.025, 0.975))
         lower_mat[j, k] <- ci[1]
         upper_mat[j, k] <- ci[2]
         cover_mat[j, k] <- (Lambda_true[j, k] >= ci[1] && Lambda_true[j, k] <= ci[2])
      }
   }
   
   label_mat <- matrix(
      paste0(
         sprintf("%.2f", lower_mat), 
         " (", sprintf("%.2f", Lambda_true), ") ", 
         sprintf("%.2f", upper_mat)
      ),
      nrow = p1, ncol = K,
      dimnames = list(
         paste0("V", 1:p1),
         paste0("F", 1:K)
      )
   )
   cover_num <- cover_mat * 1
   
   # Plot and save
   plot_title <- sprintf("%s: 95%% CI Coverage & Endpoints\nlower (true) upper", mod)
   pheatmap(
      cover_num,
      color           = c("grey80", "grey13"),
      cluster_rows    = FALSE,
      cluster_cols    = FALSE,
      display_numbers = label_mat,
      number_color    = "white",
      main            = plot_title,
      fontsize_number = 8
   )
}
