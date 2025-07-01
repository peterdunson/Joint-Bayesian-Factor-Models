# -------------------------------------------------------------------
# One-Factor MoM Diagnostics & Identifiability Checks (Ciprian Notes)
# -------------------------------------------------------------------

# 1. Load and standardize your data
# (Replace with actual data for simulation or NHANES)
# dat <- readRDS("/path/to/data.rds")
X <- scale(dat, center = TRUE, scale = TRUE)   # n x p

n <- nrow(X)
p <- ncol(X)

# 2. Compute covariance and correlation matrices
S <- cov(X)
C <- cor(X)

# 3. Estimate bi for each individual as a linear combination of X_{pi}
#    Using Method-of-Moments estimator for lambda^2

# Estimate lambda_p^2 for each variable p using all valid trios (p, q, r)
lambda2_est <- numeric(p)
lambda2_all <- vector("list", p) # Store all trio values for each p

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
   lambda2_all[[pp]] <- vals
}
names(lambda2_est) <- colnames(X)

cat("\nEstimated squared loadings (lambda_p^2) for each variable:\n")
print(round(lambda2_est, 4))

# 4. Check if any lambda_p^2 violate upper bound: lambda_p^2 > Var(X_{ip}) = 1 (since standardized)
violations <- which(lambda2_est > 1)
if (length(violations) > 0) {
   cat("\nWARNING: The following variables have lambda_p^2 > Var(X_{ip}) = 1:\n")
   print(names(lambda2_est)[violations])
   print(round(lambda2_est[violations], 4))
} else {
   cat("\nNo lambda_p^2 exceed Var(X_{ip}) = 1 (all <= 1, as expected for true one-factor model).\n")
}

# 5. Estimate bi for each subject (using estimated loadings)
# We'll use lambda_p = sign(lambda2_est) * sqrt(abs(lambda2_est)) to retain possible negative direction
lambda_est <- sign(lambda2_est) * sqrt(abs(lambda2_est))
denom <- sum(lambda_est^2)
b_hat <- as.numeric(X %*% lambda_est) / denom

cat("\nFirst 10 estimated latent scores (b_i):\n")
print(round(head(b_hat, 10), 3))

# 6. Check normality of b_i
par(mfrow = c(1,2))
hist(
   b_hat,
   breaks = 40,
   main = "Estimated Latent Scores (b_i)",
   xlab = "b_i",
   col = "orange",
   xlim = range(b_hat)  # <-- this ensures you see all values, negative or positive
)

qqnorm(b_hat, main = "Q-Q Plot of b_i"); qqline(b_hat, col=2)
par(mfrow = c(1,1))
cat("\nShapiro-Wilk test for normality of b_i:\n")
print(shapiro.test(b_hat))

# 7. For all trios (p, s, r), check sign constraint
neg_trios <- 0
total_trios <- 0
for (p1 in 1:(p-2)) for (p2 in (p1+1):(p-1)) for (p3 in (p2+1):p) {
   cpq <- C[p1, p2]
   cps <- C[p1, p3]
   cqs <- C[p2, p3]
   prod <- cpq * cps * cqs
   if (!is.na(prod)) {
      total_trios <- total_trios + 1
      if (prod < 0) neg_trios <- neg_trios + 1
   }
}
frac_neg_trios <- neg_trios / total_trios
cat(sprintf("\nFraction of negative trios (should be near 0 for true 1-factor): %.4f\n", frac_neg_trios))

# 8. Report summary of lambda^2 distributions for each variable (optional)
cat("\nSummary statistics for all trio-based lambda^2 estimates (by variable):\n")
for (pp in 1:p) {
   vals <- lambda2_all[[pp]]
   cat(sprintf("%s: mean = %.4f, sd = %.4f, min = %.4f, max = %.4f, n_trios = %d\n",
               colnames(X)[pp], mean(vals, na.rm=TRUE), sd(vals, na.rm=TRUE),
               min(vals, na.rm=TRUE), max(vals, na.rm=TRUE), length(vals)))
}

# Done!




# Assuming X is your n x p data matrix (centered and scaled)

par(mfrow = c(1, 1))  # reset to a single plot


get_offdiag <- function(M) M[lower.tri(M)]

# (1) Compute observed off-diagonal correlations
R_obs <- cor(X)
r_obs <- get_offdiag(R_obs)

# (2) Permutation null: permute each column independently
set.seed(123)
X_perm <- apply(X, 2, sample)
R_perm <- cor(X_perm)
r_perm <- get_offdiag(R_perm)

# (3) Plot overlayed histograms
breaks <- seq(-1, 1, length.out = 41)
hist(r_perm, breaks = breaks, col = rgb(0.2, 0.4, 1, 0.5), border = "white",
     main = "Off-diagonal correlations: Observed vs Permutation Null",
     xlab = "Correlation", ylab = "Frequency", xlim = c(-1, 1))
hist(r_obs, breaks = breaks, col = rgb(1, 0.4, 0.4, 0.5), add = TRUE, border = "white")
legend("topright", legend = c("Permutation Null", "Observed"),
       fill = c(rgb(0.2, 0.4, 1, 0.5), rgb(1, 0.4, 0.4, 0.5)), border = NA)










# Assume you already have: X <- scale(dat, center = TRUE, scale = TRUE)
n <- nrow(X)
p <- ncol(X)
C <- cor(X)

lambda2_est2 <- numeric(p)
for (pp in 1:p) {
   vals <- c()
   for (qq in setdiff(1:p, pp)) {
      c_pq <- C[pp, qq]
      c_qq <- C[qq, qq]
      if (!is.na(c_pq) && abs(c_pq) > 1e-8 && !is.na(c_qq) && abs(c_qq) > 1e-8) {
         vals <- c(vals, (c_pq^2) / c_qq)
      }
   }
   lambda2_est2[pp] <- mean(vals, na.rm = TRUE)
}
names(lambda2_est2) <- colnames(X)
cat("\nSecond estimator of squared loadings (lambda_p^2) for each variable:\n")
print(round(lambda2_est2, 4))


lambda2_est_nz <- numeric(p)
for (pp in 1:p) {
   vals <- c()
   for (qq in setdiff(1:p, pp)) {
      for (rr in setdiff(1:p, c(pp, qq))) {
         c_pq <- C[pp, qq]
         c_pr <- C[pp, rr]
         c_qr <- C[qq, rr]
         if (!is.na(c_pq) && !is.na(c_pr) && !is.na(c_qr) &&
             abs(c_pq) > 1e-8 && abs(c_pr) > 1e-8 && abs(c_qr) > 1e-8) {
            vals <- c(vals, (c_pq * c_pr) / c_qr)
         }
      }
   }
   lambda2_est_nz[pp] <- mean(vals, na.rm = TRUE)
}
names(lambda2_est_nz) <- colnames(X)
cat("\nTrio-based estimator (nonzero correlations only) for lambda_p^2:\n")
print(round(lambda2_est_nz, 4))


# (a) Using first estimator (original trios)
lambda_est1 <- sign(lambda2_est) * sqrt(abs(lambda2_est))
denom1 <- sum(lambda_est1^2)
b_hat1 <- as.numeric(X %*% lambda_est1) / denom1

# (b) Using second estimator
lambda_est2 <- sign(lambda2_est2) * sqrt(abs(lambda2_est2))
denom2 <- sum(lambda_est2^2)
b_hat2 <- as.numeric(X %*% lambda_est2) / denom2

# (c) Using nonzero-correlation trio estimator
lambda_est_nz <- sign(lambda2_est_nz) * sqrt(abs(lambda2_est_nz))
denom_nz <- sum(lambda_est_nz^2)
b_hat_nz <- as.numeric(X %*% lambda_est_nz) / denom_nz

# Summaries
cat("\nFirst 10 estimated latent scores (b_i) for each estimator:\n")
cat("Original trio estimator:\n")
print(round(head(b_hat1, 10), 3))
cat("Second estimator (pairwise ratio):\n")
print(round(head(b_hat2, 10), 3))
cat("Nonzero-correlation trios estimator:\n")
print(round(head(b_hat_nz, 10), 3))


cat("\nWeights for each variable (lambda_p) under different estimators:\n")
weights <- data.frame(
   variable = colnames(X),
   lambda_trio = round(lambda_est1, 4),
   lambda_pairwise = round(lambda_est2, 4),
   lambda_trio_nz = round(lambda_est_nz, 4)
)
print(weights)


cor_b1 <- cor(b_hat1, X)
cor_b2 <- cor(b_hat2, X)
cor_bnz <- cor(b_hat_nz, X)

cat("\nCorrelation of estimated b_i with each variable (original trios estimator):\n")
print(round(cor_b1, 3))

cat("\nCorrelation of estimated b_i with each variable (pairwise ratio estimator):\n")
print(round(cor_b2, 3))

cat("\nCorrelation of estimated b_i with each variable (nonzero trios estimator):\n")
print(round(cor_bnz, 3))






# ----- Setup -----
# (Assume X, p, lambda_est1/2/nz, b_hat1/2/nz, cor_b1/2/nz are all in workspace)

# ---- 1. SCATTER PLOTS: Each X_p vs estimated b_i for all three estimators ----

par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))  # One var at a time, all estimators stacked
est_names <- c("Trio", "Pairwise", "NZ")
lambda_list <- list(lambda_est1, lambda_est2, lambda_est_nz)
b_list     <- list(b_hat1, b_hat2, b_hat_nz)

for (j in 1:p) {
   for (k in 1:3) {
      plot(
         X[, j], b_list[[k]],
         xlab = colnames(X)[j], ylab = "b_i",
         main = paste0(est_names[k], ": ", colnames(X)[j])
      )
   }
   # Pause between variables (uncomment to step through manually)
   # readline(prompt = "Press [enter] for next variable...")
}

par(mfrow = c(1, 1))

# ---- 2. BARPLOTS: Correlation of b_hat with each variable ----

barplot(as.numeric(cor_b1), names.arg = colnames(X), main = "Correlation of b_i (Trio) with X_p", las = 2)
barplot(as.numeric(cor_b2), names.arg = colnames(X), main = "Correlation of b_i (Pairwise) with X_p", las = 2)
barplot(as.numeric(cor_bnz), names.arg = colnames(X), main = "Correlation of b_i (NZ) with X_p", las = 2)

# ---- 3. HEATMAP: Fitted X (b_i * lambda_p), e.g. Trio estimator ----
library(pheatmap)
# Just the first 50 rows
n_show <- min(50, nrow(X_fit1))
pheatmap(
   X_fit1[1:n_show, ],
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   main = "Fitted X (First 50 subjects)"
)

