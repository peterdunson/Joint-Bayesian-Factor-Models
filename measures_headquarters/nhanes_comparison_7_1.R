# -------------------------------------------------------------------
# NHANES: Trio-Based Factor Diagnostics, Estimators & Decorrelation
# -------------------------------------------------------------------

# 1. Load and standardize your data
X <- scale(dat, center = TRUE, scale = TRUE)   # n x p
n <- nrow(X)
p <- ncol(X)
C <- cor(X)
S <- cov(X)

# -------------------------
# 2. Method of Moments Estimators (Trio-based and alternatives)
# -------------------------

lambda2_est <- numeric(p)
lambda2_all <- vector("list", p)
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
cat("\nEstimated squared loadings (lambda_p^2):\n")
print(round(lambda2_est, 4))

# Alternative estimator (pairwise ratio)
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
cat("\nSecond estimator of squared loadings (lambda_p^2):\n")
print(round(lambda2_est2, 4))

# Trio-based estimator using only nonzero correlations
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

# 3. Estimate b_i for each subject for all estimators
lambda_est1 <- sign(lambda2_est) * sqrt(abs(lambda2_est))
denom1 <- sum(lambda_est1^2)
b_hat1 <- as.numeric(X %*% lambda_est1) / denom1

lambda_est2 <- sign(lambda2_est2) * sqrt(abs(lambda2_est2))
denom2 <- sum(lambda_est2^2)
b_hat2 <- as.numeric(X %*% lambda_est2) / denom2

lambda_est_nz <- sign(lambda2_est_nz) * sqrt(abs(lambda2_est_nz))
denom_nz <- sum(lambda_est_nz^2)
b_hat_nz <- as.numeric(X %*% lambda_est_nz) / denom_nz

# 4. Show summary for each estimator
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

# 5. Fraction of negative trios
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
cat(sprintf("\nFraction of negative trios: %.4f\n", frac_neg_trios))

# 6. Summary statistics for all trio-based lambda^2 estimates
cat("\nSummary statistics for all trio-based lambda^2 estimates (by variable):\n")
for (pp in 1:p) {
   vals <- lambda2_all[[pp]]
   cat(sprintf("%s: mean = %.4f, sd = %.4f, min = %.4f, max = %.4f, n_trios = %d\n",
               colnames(X)[pp], mean(vals, na.rm=TRUE), sd(vals, na.rm=TRUE),
               min(vals, na.rm=TRUE), max(vals, na.rm=TRUE), length(vals)))
}

# ------------------------
# 7. MoM vs MGSP decorrelation comparison (K=1, as above)
# ------------------------

# --- MoM residuals ---
X_fit_mom <- outer(b_hat1, lambda_est1)
resid_mom <- X - X_fit_mom

# --- Bayesian MGSP residuals ---
fit_MGSP <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/direct_paper_fits/storing_fit/MGSP/fit_joint_NHANES1718_k1.rds")
Lambda_hat_mgsp <- fit_MGSP$Lambda_hat
if (nrow(Lambda_hat_mgsp) != ncol(X)) Lambda_hat_mgsp <- t(Lambda_hat_mgsp)
eta_hat_mgsp <- X %*% Lambda_hat_mgsp %*% solve(t(Lambda_hat_mgsp) %*% Lambda_hat_mgsp)
X_fit_mgsp <- eta_hat_mgsp %*% t(Lambda_hat_mgsp)
resid_mgsp <- X - X_fit_mgsp

# --- Overlayed histogram plot ---
get_offdiag <- function(M) M[lower.tri(M)]
r_obs <- get_offdiag(cor(X))
r_mom <- get_offdiag(cor(resid_mom))
r_mgsp <- get_offdiag(cor(resid_mgsp))
breaks <- seq(-1, 1, length.out = 41)

hist(r_obs,  breaks = breaks, col = rgb(0.2, 0.4, 1, 0.4), border = "#003399",
     main = "Off-diagonal Correlations: NHANES",
     xlab = "Correlation", xlim = c(-1, 1))
hist(r_mom,  breaks = breaks, col = rgb(1, 0.4, 0.4, 0.5), add = TRUE, border = "#BB2222")
hist(r_mgsp, breaks = breaks, col = rgb(0.2, 1, 0.6, 0.3), add = TRUE, border = "#008844")
legend("topright",
       legend = c("Original", "MoM K=1", "Bayesian K=1 (MGSP)"),
       fill = c(rgb(0.2,0.4,1,0.4), rgb(1,0.4,0.4,0.5), rgb(0.2,1,0.6,0.3)),
       border = NA)

cat("\nMean absolute off-diagonal correlation:\n")
cat(sprintf("Original:        %.4f\n", mean(abs(r_obs))))
cat(sprintf("MoM K=1:         %.4f\n", mean(abs(r_mom))))
cat(sprintf("Bayesian K=1:    %.4f\n", mean(abs(r_mgsp))))

