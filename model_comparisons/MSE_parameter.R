
library(rstan)
library(pheatmap)

# ---- File paths and settings ----
sim_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations"
scenario <- 1
n_train <- 1000
K <- 5

sim_file <- sprintf("%s/sim_scen%d_1000.rds", sim_dir, scenario)
mgps_file <- sprintf("%s/mgps_fit_scen%d_5.rds", sim_dir, scenario)
tebfar_file <- sprintf("%s/stan_tebfar_fit_scen%d_5.rds", sim_dir, scenario)

# ---- Load simulation and true loadings ----
sim <- readRDS(sim_file)
Lambda_true <- sim$Lambda  # [p+1, K] (usually first row is y)

# Drop outcome if needed:
Lambda_true_pred <- Lambda_true[-1, , drop=FALSE] # predictors only, [p x K]

# ---- Helper: Extract posterior mean loadings ----
get_mean_loadings <- function(fit_file) {
   fit <- readRDS(fit_file)
   post <- rstan::extract(fit)
   Lambda_mean <- apply(post$Lambda, c(2,3), mean)
   # Drop outcome row if present
   if (nrow(Lambda_mean) == nrow(Lambda_true)) {
      Lambda_mean <- Lambda_mean[-1, , drop=FALSE]
   }
   return(Lambda_mean)
}

# ---- Extract posterior mean loadings for each model ----
Lambda_mgps   <- get_mean_loadings(mgps_file)
Lambda_tebfar <- get_mean_loadings(tebfar_file)

# ---- Compute per-parameter MSE matrices ----
mse_matrix_mgps   <- (Lambda_mgps   - Lambda_true_pred)^2
mse_matrix_tebfar <- (Lambda_tebfar - Lambda_true_pred)^2

# ---- Plot heatmaps of per-parameter MSE ----
pheatmap(mse_matrix_mgps,
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "MGPS Per-Parameter MSE",
         labels_row = paste0("Var", 1:nrow(mse_matrix_mgps)),
         labels_col = paste0("F", 1:K),
         color = colorRampPalette(c("white", "orange", "red"))(50))

pheatmap(mse_matrix_tebfar,
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "TEB-FAR Per-Parameter MSE",
         labels_row = paste0("Var", 1:nrow(mse_matrix_tebfar)),
         labels_col = paste0("F", 1:K),
         color = colorRampPalette(c("white", "orange", "red"))(50))

# ---- Optionally, summarize mean MSE per variable/factor ----
cat("MGPS mean MSE per variable:", round(rowMeans(mse_matrix_mgps), 3), "\n")
cat("TEB-FAR mean MSE per variable:", round(rowMeans(mse_matrix_tebfar), 3), "\n")
cat("MGPS overall MSE:", mean(mse_matrix_mgps), "\n")
cat("TEB-FAR overall MSE:", mean(mse_matrix_tebfar), "\n")




