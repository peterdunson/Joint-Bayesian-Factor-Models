# evaluate_factor_measures.R

#CIPRIAN SUGGESTIONS:

# R2 and residual correlation reduction
# Correlation reduction between factors (η) and errors (ϵ)
# MSE for zero and nonzero coefficients
# Coverage for individual parameters

# ─── 1) Libraries ─────────────────────────────────────────────────────────
library(mvtnorm)    # for cor()
library(dplyr)

# ─── 2) Load sim truth & fit ──────────────────────────────────────────────
sim   <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Lambda_true <- sim$Lambda           # (p+1) × K
Sigma_true  <- diag(rep(0.2, nrow(Lambda_true)))
Omega_true  <- sim$Omega            # true covariance
Y           <- sim$Y                # original data

obj    <- readRDS("fit_Joint_scen2_scale_all_6k.rds")
fit    <- if (inherits(obj, "stanfit")) obj else obj$fit
post   <- if (is.list(obj) && !inherits(obj,"stanfit")) obj$posterior else extract(fit)

# ─── 3) Estimate posterior‐mean Λ, η and Ω ─────────────────────────────────
Lambda_hat <- apply(post$Lambda, c(2,3), mean)   # [p+1]×K
eta_hat    <- apply(post$eta,    c(2,3), mean)   # [n]×K
Omega_hat  <- Lambda_hat %*% t(Lambda_hat) + Sigma_true

# R² of covariance
R2_cov <- cor(as.vector(Omega_hat), as.vector(Omega_true))^2

# ─── 4) Residual correlation reduction ────────────────────────────────────
# compute raw correlation
corr_raw <- cor(Y)
# reconstructed signal
Y_hat    <- eta_hat %*% t(Lambda_hat)
resid    <- Y - Y_hat
corr_res <- cor(resid)

# average off‐diagonal magnitude
offdiag_mean <- function(M) mean(abs(M[lower.tri(M)]))
raw_off   <- offdiag_mean(corr_raw)
resid_off <- offdiag_mean(corr_res)
reduction <- (raw_off - resid_off) / raw_off

# ─── 5) Factor–error correlation ──────────────────────────────────────────
# compute for posterior‐mean only (you could also do draws)
errs <- resid
fac_err_cor <- max(abs(cor(eta_hat, errs)))

# ─── 6) MSE zero vs. nonzero λ ────────────────────────────────────────────
p1 <- nrow(Lambda_true); K <- ncol(Lambda_true)
# true pattern
pattern <- (Lambda_true != 0)
# posterior‐mean MSE
mse_mat <- (Lambda_hat - Lambda_true)^2
mse_zero    <- mean(mse_mat[pattern == FALSE])
mse_nonzero <- mean(mse_mat[pattern == TRUE])

# ─── 7) 95%‐CI coverage for each λ ────────────────────────────────────────
# extract all draws into [iter, j, k]
Lam_draws <- post$Lambda  # dims: iter × j × k
# build a logical matrix of coverage
cover <- matrix(NA, nrow = p1, ncol = K)
for (j in 1:p1) for (k in 1:K) {
   draws_jk <- Lam_draws[,j,k]
   ci       <- quantile(draws_jk, c(0.025, 0.975))
   cover[j,k] <- (Lambda_true[j,k] >= ci[1] && Lambda_true[j,k] <= ci[2])
}
coverage_rate <- mean(cover)

# ─── 8) Report ────────────────────────────────────────────────────────────
cat(sprintf("R²(Omega)              = %.4f\n",        R2_cov))
cat(sprintf("Corr. off‐diag mean: raw = %.4f, resid = %.4f (reduction = %.2f%%)\n",
            raw_off, resid_off, reduction*100))
cat(sprintf("Max |cor(eta,eps)|     = %.4f\n",        fac_err_cor))
cat(sprintf("MSE λ (zero entries)   = %.6f\n",        mse_zero))
cat(sprintf("MSE λ (nonzero entries)= %.6f\n",        mse_nonzero))
cat(sprintf("95%%‐CI coverage rate   = %.2f%%\n",     coverage_rate*100))



# ─── After your coverage loop in evaluate_factor_measures.R ────────────────

# convert logical cover[ j , k ] → numeric (0/1)
cov_mat <- cover * 1
rownames(cov_mat) <- paste0("V", seq_len(nrow(cov_mat)))
colnames(cov_mat) <- paste0("F", seq_len(ncol(cov_mat)))

# draw a binary heatmap: 1=covered, 0=missed
library(pheatmap)
pheatmap(
   cov_mat,
   color        = c("grey","grey14"),
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   legend       = TRUE,
   main         = "95% CI Coverage for Each λ₍j,k₎",
   display_numbers = TRUE,
   number_color    = "white"
)




# ─── assume post and Lambda_true are already loaded ─────────────────────────

p1 <- nrow(Lambda_true)
K  <- ncol(Lambda_true)

# 1) Preallocate
cover_mat <- matrix(0, nrow = p1, ncol = K)
lower_mat <- matrix(NA, nrow = p1, ncol = K)
upper_mat <- matrix(NA, nrow = p1, ncol = K)

# 2) Fill them
for (j in seq_len(p1)) {
   for (k in seq_len(K)) {
      draws <- post$Lambda[, j, k]
      qs    <- quantile(draws, c(0.025, 0.975))
      lower_mat[j,k] <- qs[1]
      upper_mat[j,k] <- qs[2]
      cover_mat[j,k] <- (Lambda_true[j,k] >= qs[1] && Lambda_true[j,k] <= qs[2])
   }
}

rownames(cover_mat) <- paste0("V", 1:p1)
colnames(cover_mat) <- paste0("F", 1:K)
rownames(lower_mat) <- rownames(upper_mat) <- rownames(cover_mat)
colnames(lower_mat) <- colnames(upper_mat) <- colnames(cover_mat)

# 3) Build label matrix "lower, upper" on one line
label_mat <- matrix(
   paste0(
      sprintf("%.2f", lower_mat),
      ", ",
      sprintf("%.2f", upper_mat)
   ),
   nrow = p1, ncol = K,
   dimnames = dimnames(cover_mat)
)

# 4) Draw heatmap
library(pheatmap)
pheatmap(
   cover_mat,
   color           = c("grey","grey12"),
   cluster_rows    = FALSE,
   cluster_cols    = FALSE,
   display_numbers = label_mat,
   number_color    = "white",
   main            = "95% CI Coverage (color) & Credible Interval Endpoints"
)











#NEWWWWW

# ─── assume `post` (list with posterior samples) and `Lambda_true` are already in your workspace ────────────────

p1 <- nrow(Lambda_true)
K  <- ncol(Lambda_true)

# 1) Preallocate matrices
cover_mat <- matrix(FALSE,  nrow = p1, ncol = K)
lower_mat <- matrix(NA_real_, nrow = p1, ncol = K)
upper_mat <- matrix(NA_real_, nrow = p1, ncol = K)

# 2) Compute 95% CIs and coverage
for (j in seq_len(p1)) {
   for (k in seq_len(K)) {
      draws <- post$Lambda[, j, k]
      ci    <- quantile(draws, c(0.025, 0.975))
      lower_mat[j, k] <- ci[1]
      upper_mat[j, k] <- ci[2]
      cover_mat[j, k] <- (Lambda_true[j, k] >= ci[1] && Lambda_true[j, k] <= ci[2])
   }
}

# 3) Build a label matrix "lower (true) upper"
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

# 4) Draw the heatmap
library(pheatmap)
# turn TRUE/FALSE into 1/0
cover_num <- cover_mat * 1

pheatmap(
   cover_num,
   color           = c("grey80", "grey13"),      # 0→grey, 1→navy
   cluster_rows    = FALSE,
   cluster_cols    = FALSE,
   display_numbers = label_mat,
   number_color    = "white",
   main            = "95% CI Coverage & Interval Endpoints\nlower (true) upper"
)















fit <- if (inherits(obj, "stanfit")) obj else obj$fit
library(shinystan)

launch_shinystan(fit)
