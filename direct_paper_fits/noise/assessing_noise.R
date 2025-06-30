getwd()

# ----- Load results -----
fit_res <- readRDS("fit_Joint_scen2_1000_noise_1.rds")
Lambda_hat <- fit_res$Lambda_hat
post_j <- fit_res$posterior

# ----- Plot 1: Heatmap of Posterior Mean Loadings -----
library(pheatmap)
pheatmap(
   as.matrix(Lambda_hat),
   color = colorRampPalette(c("blue", "white", "red"))(100),
   cluster_rows = FALSE, cluster_cols = FALSE,
   main = "Posterior Mean Loadings (K=1, Pure Noise)"
)

# ----- Plot 2: Boxplot of Posterior Mean Loadings -----
boxplot(as.vector(Lambda_hat), main = "Distribution of Posterior Mean Loadings", ylab = "Value")

# ----- Plot 3: Histogram of All Posterior Draws (Flattened) -----
Lambda_draws <- as.vector(post_j$Lambda)
hist(Lambda_draws, breaks = 30, main = "All Posterior Draws of Loadings", xlab = "Loading Value", col = "gray")

# ----- Plot 4: Residual Correlations Before & After -----
Y <- scale(readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000_noise.rds")$Y)
cor_before <- cor(Y)
get_offdiag <- function(mat) mat[lower.tri(mat)]
r_before <- get_offdiag(cor_before)

# Project out K=1 factor
eta_hat1 <- Y %*% Lambda_hat %*% solve(t(Lambda_hat) %*% Lambda_hat)
Y_hat1   <- eta_hat1 %*% t(Lambda_hat)
resid1   <- Y - Y_hat1
cor_after1 <- cor(resid1)
r_after1 <- get_offdiag(cor_after1)

hist(r_before, breaks = 30, col = rgb(0.2, 0.4, 1, 0.5), main = "Off-diag Correlations\nBefore (blue) / After K=1 (orange)", xlab = "Correlation", border = "white")
hist(r_after1, breaks = 30, col = rgb(1, 0.7, 0.2, 0.5), add = TRUE, border = "white")
legend("topright", legend = c("Before", "After K=1"), fill = c(rgb(0.2, 0.4, 1, 0.5), rgb(1, 0.7, 0.2, 0.5)), border = NA)

# ----- Optional: Print summary of loadings -----
cat("Summary of posterior mean loadings:\n")
print(summary(as.vector(Lambda_hat)))

