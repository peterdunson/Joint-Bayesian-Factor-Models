library(rstan)
library(pheatmap)

# ---- User settings ----
#fit_file  <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/mgps_fit_scen1_5.rds"  # or tebfar file
sim_file  <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen1_1000.rds"
fit_file <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/stan_tebfar_fit_scen1_5.rds"
model_name <- "TEB-FAR1"   # or "TEB-FAR"

# ---- Load data ----
fit <- readRDS(fit_file)
sim <- readRDS(sim_file)
post <- rstan::extract(fit)
Lambda_true <- sim$Lambda      # [num_vars, K]
Lambda_samples <- post$Lambda  # [iterations, num_vars, K]
p <- dim(Lambda_true)[1]
K <- dim(Lambda_true)[2]

# ---- Calculate 95% CI coverage matrix ----
coverage_matrix <- matrix(FALSE, nrow = p, ncol = K)
for (j in 1:p) {
   for (k in 1:K) {
      ci <- quantile(Lambda_samples[, j, k], c(0.025, 0.975))
      coverage_matrix[j, k] <- (Lambda_true[j, k] >= ci[1]) & (Lambda_true[j, k] <= ci[2])
   }
}

# ---- Plot coverage heatmap ----
pheatmap(
   coverage_matrix * 1,  # Convert logical to numeric (0/1)
   color = c("grey", "grey12"),
   cluster_rows = FALSE, cluster_cols = FALSE,
   legend_breaks = c(0,1), legend_labels = c("Not covered", "Covered"),
   main = sprintf("%s Per-Variable/Factor 95%% Coverage", model_name),
   labels_row = paste0("Var", 1:p),
   labels_col = paste0("F", 1:K)
)


