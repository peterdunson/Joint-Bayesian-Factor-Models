# -------------------------------------------------------------------
# One-Factor MoM Diagnostics & Identifiability Checks (Ciprian Notes)
# For a single simulation replicate
# -------------------------------------------------------------------

# 1. Load and standardize your data (SIMULATION)
sim_path <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000_1.rds"
sim2 <- readRDS(sim_path)
X <- scale(sim2$Y, center = TRUE, scale = TRUE)   # n x p
n <- nrow(X)
p <- ncol(X)

# 2. Compute Aacovariance and correlation matrices
S <- cov(X)
C <- cor(X)

# 3. Estimate lambda_p^2 for each variable using all valid trios
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

# 4. Check if any lambda_p^2 violate the upper bound: lambda_p^2 > Var(X_{ip}) = 1
violations <- which(lambda2_est > 1)
if (length(violations) > 0) {
   cat("\nWARNING: The following variables have lambda_p^2 > 1:\n")
   print(names(lambda2_est)[violations])
   print(round(lambda2_est[violations], 4))
} else {
   cat("\nNo lambda_p^2 exceed 1 (all <= 1).\n")
}

# 5. Estimate b_i for each subject (using estimated loadings)
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
   xlim = range(b_hat)
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
cat(sprintf("\nFraction of negative trios: %.4f\n", frac_neg_trios))

# 8. Report summary of lambda^2 distributions for each variable
cat("\nSummary statistics for all trio-based lambda^2 estimates (by variable):\n")
for (pp in 1:p) {
   vals <- lambda2_all[[pp]]
   cat(sprintf("%s: mean = %.4f, sd = %.4f, min = %.4f, max = %.4f, n_trios = %d\n",
               colnames(X)[pp], mean(vals, na.rm=TRUE), sd(vals, na.rm=TRUE),
               min(vals, na.rm=TRUE), max(vals, na.rm=TRUE), length(vals)))
}

# 9. Off-diagonal correlation distribution (observed vs permutation null)
get_offdiag <- function(M) M[lower.tri(M)]
R_obs <- cor(X)
r_obs <- get_offdiag(R_obs)
set.seed(123)
X_perm <- apply(X, 2, sample)
R_perm <- cor(X_perm)
r_perm <- get_offdiag(R_perm)

breaks <- seq(-1, 1, length.out = 41)
hist(r_perm, breaks = breaks, col = rgb(0.2, 0.4, 1, 0.5), border = "white",
     main = "Off-diagonal correlations: Observed vs Permutation Null",
     xlab = "Correlation", ylab = "Frequency", xlim = c(-1, 1))
hist(r_obs, breaks = breaks, col = rgb(1, 0.4, 0.4, 0.5), add = TRUE, border = "white")
legend("topright", legend = c("Permutation Null", "Observed"),
       fill = c(rgb(0.2, 0.4, 1, 0.5), rgb(1, 0.4, 0.4, 0.5)), border = NA)

# Done! All outputs and checks are included above.
