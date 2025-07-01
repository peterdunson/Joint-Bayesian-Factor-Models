# --- NHANES One-Factor Trio Analysis ---

# (1) Load your data (replace with your path and object if needed)
# dat <- readRDS("/path/to/nhanes_data.rds")

# (2) Standardize: rows = subjects, columns = variables
X <- scale(dat, center = TRUE, scale = TRUE)

# (3) Covariance and correlation
S <- cov(X)
C <- cor(X)
p <- ncol(X)

# (4) Trio constraint: fraction of negative trios
neg_count <- 0
total_trios <- 0

for (i in 1:(p-2)) {
   for (j in (i+1):(p-1)) {
      for (k in (j+1):p) {
         cpq <- C[i, j]
         cpr <- C[i, k]
         cqr <- C[j, k]
         prod <- cpq * cpr * cqr
         if (!is.na(prod)) {
            total_trios <- total_trios + 1
            if (prod < 0) neg_count <- neg_count + 1
         }
      }
   }
}
frac_neg <- neg_count / total_trios
cat(sprintf("Fraction of negative trios (should be near 0 for true 1-factor): %.4f\n", frac_neg))

# (5) Estimate squared loadings λ_p^2 for each variable
lambda2_est <- numeric(p)
for (pp in 1:p) {
   vals <- c()
   for (qq in setdiff(1:p, pp)) {
      for (rr in setdiff(1:p, c(pp, qq))) {
         c_pq <- C[pp, qq]
         c_pr <- C[pp, rr]
         c_qr <- C[qq, rr]
         if (!is.na(c_qr) && abs(c_qr) > 1e-8) {
            vals <- c(vals, (c_pq * c_pr) / c_qr)
         }
      }
   }
   lambda2_est[pp] <- mean(vals, na.rm = TRUE)
}
names(lambda2_est) <- colnames(X)
cat("\nEstimated squared loadings (λ_p^2):\n")
print(round(lambda2_est, 3))

# (6) Estimate latent scores b_i for each subject
# We'll use λ_p = sign(sum(X[,p])) * sqrt(λ_p^2) to pick a direction, or just sqrt(λ_p^2)
lambda_est <- sqrt(abs(lambda2_est)) * sign(lambda2_est)
denom <- sum(lambda_est^2)
b_hat <- as.numeric(X %*% lambda_est) / denom

cat("\nFirst 10 estimated latent scores (b_i):\n")
print(round(head(b_hat, 10), 3))



# Show empirical covariance/correlation matrix for NHANES
cat("Empirical correlation matrix (NHANES):\n")
print(round(cor(X), 3))

# Count negatives
off_diag <- function(M) M[lower.tri(M)]
num_neg <- sum(off_diag(cor(X)) < 0)
frac_neg <- num_neg / length(off_diag(cor(X)))
cat(sprintf("Fraction of negative off-diagonal correlations: %.4f\n", frac_neg))


print(round(cor(X), 3))



# Assume X is your data matrix, already scaled
cor_mat <- cor(X)
rounded_cor <- round(cor_mat, 3)

# Basic heatmap using pheatmap
library(pheatmap)

# Choose a color palette that clearly shows sign (blue = negative, white = zero, red = positive)
pheatmap(
   rounded_cor,
   color = colorRampPalette(c("white", "pink", "red"))(100),
   main = "Correlation Matrix Heatmap",
   display_numbers = TRUE,
   number_color = "black",    # Numbers are visible on all backgrounds
   cluster_rows = FALSE,
   cluster_cols = FALSE
)


