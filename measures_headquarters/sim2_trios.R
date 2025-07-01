# -------------------------------------------------------------------
# Method-of-Moments 1-Factor Analysis: Simulation 2
# Implements Ciprian's notes (MoM estimation, identifiability checks)
# -------------------------------------------------------------------

# 1. Load simulation datax
sim_path <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000_1.rds"
sim2 <- readRDS(sim_path)
X <- scale(sim2$Y, center = TRUE, scale = TRUE)  # n x p data matrix

n <- nrow(X)
p <- ncol(X)

# 2. Compute covariance matrix
S <- cov(X)   # p x p

# 3. Check trio constraint: for all trios (p, q, r), compute c_{pq} c_{pr} c_{qr}
num_trios <- 0
num_neg   <- 0
for (p1 in 1:(p-2)) for (p2 in (p1+1):(p-1)) for (p3 in (p2+1):p) {
   c12 <- S[p1, p2]
   c13 <- S[p1, p3]
   c23 <- S[p2, p3]
   if (!is.na(c12) && !is.na(c13) && !is.na(c23)) {
      prod <- c12 * c13 * c23
      num_trios <- num_trios + 1
      if (prod < 0) num_neg <- num_neg + 1
   }
}
cat(sprintf("Fraction of negative trios (should be near 0 for true 1-factor): %.4f\n", num_neg / num_trios))

# 4. MoM Estimate: lambda_p^2 for each variable
lambda2_est <- rep(NA, p)
for (pp in 1:p) {
   inds_qr <- setdiff(1:p, pp)
   vals <- c()
   for (q in inds_qr) for (r in inds_qr) if (q < r && S[q, r] != 0) {
      numer <- S[pp, q] * S[pp, r]
      denom <- S[q, r]
      vals <- c(vals, numer / denom)
   }
   # Remove any infinite or NA values
   vals <- vals[is.finite(vals)]
   lambda2_est[pp] <- mean(vals)
}
# (Optionally force non-negativity)
lambda2_est[lambda2_est < 0] <- 0

# 5. Recover lambda_p (set sign using principal eigenvector if you want, else just sqrt)
lambda_est <- sqrt(lambda2_est)
names(lambda_est) <- paste0("V", 1:p)

# 6. Estimate b_i for each subject: b_i = sum_p lambda_p * X_{i,p} / sum_p lambda_p^2
denom <- sum(lambda_est^2)
b_hat <- as.vector(X %*% lambda_est) / denom

# 7. Plots & Outputs

# 7a. Plot estimated loadings
barplot(lambda_est, main = "MoM Estimated Factor Loadings (sqrt)", xlab = "Variable", ylab = "Loading", col = "skyblue")

# 7b. Histogram of b_i (latent score estimates)
hist(b_hat, breaks = 40, main = "Estimated Latent Scores (b_i) via MoM", xlab = "b_i", col = "orange")

# 7c. Overlay normal density for visual non-Gaussianity
curve(dnorm(x, mean(b_hat), sd(b_hat)) * length(b_hat) * diff(hist(b_hat, plot = FALSE)$breaks[1:2]), 
      add = TRUE, col = "darkblue", lwd = 2)

# 8. Print summary stats and comments
cat("First 10 estimated loadings (MoM):\n")
print(round(lambda_est[1:10], 4))

cat("Mean, SD, skewness, kurtosis of b_i:\n")
library(moments)
cat(sprintf("Mean: %.3f  SD: %.3f  Skew: %.3f  Kurtosis: %.3f\n",
            mean(b_hat), sd(b_hat), skewness(b_hat), kurtosis(b_hat)))



# --- Correlation matrix: print and heatmap (using your color scheme) ---
cor_mat <- cor(X)
rounded_cor <- round(cor_mat, 3)

cat("Empirical correlation matrix:\n")
print(rounded_cor)

library(pheatmap)
pheatmap(
   rounded_cor,
   color = colorRampPalette(c("white", "pink", "red"))(100),
   main = "Correlation Matrix Heatmap Sim2",
   display_numbers = TRUE,
   number_color = "black",
   cluster_rows = FALSE,
   cluster_cols = FALSE
)


