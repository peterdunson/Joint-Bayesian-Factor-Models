# Load fits
fit_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/direct_paper_fits/storing_fit"
fit_MGSP <- readRDS(file.path(fit_dir, "MGSP", "fit_joint_NHANES1718_k1.rds"))
fit_HS   <- readRDS(file.path(fit_dir, "Horseshoe", "fit_HS_NHANES1718_k1.rds"))
fit_SSL  <- readRDS(file.path(fit_dir, "spike_slab", "fit_SSL_NHANES1718_k1.rds"))

# Load and standardize data
X <- scale(dat, center = TRUE, scale = TRUE) # n x p
n <- nrow(X)
p <- ncol(X)

# --- MoM Trios estimator for lambdas ---
C <- cor(X)
lambda2_est_trios <- numeric(p)
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
   lambda2_est_trios[pp] <- mean(vals, na.rm = TRUE)
}
lambda_est_mom_trios <- sign(lambda2_est_trios) * sqrt(abs(lambda2_est_trios))
denom_mom_trios <- sum(lambda_est_mom_trios^2)
b_hat_mom_trios <- as.numeric(X %*% lambda_est_mom_trios) / denom_mom_trios

# --- MoM Eigen (Principal Component) estimator for lambdas ---
S   <- cov(X)
e   <- eigen(S)
d   <- e$values
V   <- e$vectors
sigma2_hat  <- mean(d[-1])
lambda_mag  <- sqrt(max(d[1] - sigma2_hat, 0))
lambda1_hat_eigen <- lambda_mag * V[, 1]
denom_mom_eigen <- sum(lambda1_hat_eigen^2)
b_hat_mom_eigen <- as.numeric(X %*% lambda1_hat_eigen) / denom_mom_eigen

# --- For all fits: MGSP, HS, SSL ---
extract_lambda_b <- function(fit_obj, X) {
   lambda1 <- as.numeric(fit_obj$Lambda_hat[, 1])  # Ensures vector
   denom   <- sum(lambda1^2)
   b_hat   <- as.numeric(X %*% lambda1) / denom    # n-vector
   list(lambda1 = lambda1, b_hat = b_hat)
}
out_MGSP <- extract_lambda_b(fit_MGSP, X)
out_HS   <- extract_lambda_b(fit_HS,   X)
out_SSL  <- extract_lambda_b(fit_SSL,  X)

# --- Compare b_i across all methods (pairwise scatterplots) ---
b_hat_list <- list(
   MoM_trios = b_hat_mom_trios,
   MoM_eigen = b_hat_mom_eigen,
   MGSP = out_MGSP$b_hat,
   HS   = out_HS$b_hat,
   SSL  = out_SSL$b_hat
)
pairs(as.data.frame(b_hat_list), main = "b_i: All Methods")

# --- For the first 3 variables: lambda1_j * b_i for each method, pairwise plots ---
for (j in 1:3) {
   proj_df <- data.frame(
      MoM_trios = lambda_est_mom_trios[j]   * b_hat_mom_trios,
      MoM_eigen = lambda1_hat_eigen[j]      * b_hat_mom_eigen,
      MGSP      = out_MGSP$lambda1[j] * out_MGSP$b_hat,
      HS        = out_HS$lambda1[j]   * out_HS$b_hat,
      SSL       = out_SSL$lambda1[j]  * out_SSL$b_hat
   )
   pairs(proj_df, main = paste("lambda1 * b_i (variable", j, "): All Methods"))
   print(cor(proj_df))
}

# --- lambda1 comparison (MoM_trios vs others) ---
par(mfrow = c(2,2))
plot(lambda_est_mom_trios, lambda1_hat_eigen, main = "MoM_trios vs MoM_eigen: lambda1", xlab = "MoM_trios", ylab = "MoM_eigen")
plot(lambda_est_mom_trios, out_MGSP$lambda1, main = "MoM_trios vs MGSP: lambda1", xlab = "MoM_trios", ylab = "MGSP")
plot(lambda_est_mom_trios, out_HS$lambda1, main = "MoM_trios vs HS: lambda1", xlab = "MoM_trios", ylab = "HS")
plot(lambda_est_mom_trios, out_SSL$lambda1, main = "MoM_trios vs SSL: lambda1", xlab = "MoM_trios", ylab = "SSL")
par(mfrow = c(1,1))

# --- Function for epsilons (residuals for variable j)
epsi <- function(j, X, lambda1, b_hat) {
   X[,j] - lambda1[j] * b_hat
}

# ---- b_i vs epsilon_{i1} (all methods, plot separately) ----
par(mfrow = c(1,1))
plot(b_hat_mom_trios, epsi(1, X, lambda_est_mom_trios, b_hat_mom_trios), main="MoM_trios: b_i vs epsilon_{i1}", xlab="b_i", ylab="epsilon_{i1}")
plot(b_hat_mom_eigen, epsi(1, X, lambda1_hat_eigen, b_hat_mom_eigen), main="MoM_eigen: b_i vs epsilon_{i1}", xlab="b_i", ylab="epsilon_{i1}")
plot(out_MGSP$b_hat, epsi(1, X, out_MGSP$lambda1, out_MGSP$b_hat), main="MGSP: b_i vs epsilon_{i1}", xlab="b_i", ylab="epsilon_{i1}")
plot(out_HS$b_hat,   epsi(1, X, out_HS$lambda1, out_HS$b_hat), main="HS: b_i vs epsilon_{i1}", xlab="b_i", ylab="epsilon_{i1}")
plot(out_SSL$b_hat,  epsi(1, X, out_SSL$lambda1, out_SSL$b_hat), main="SSL: b_i vs epsilon_{i1}", xlab="b_i", ylab="epsilon_{i1}")
par(mfrow = c(1,1))





# --- Residuals after rank-1 removal for all methods ---
resid_mom_trios   <- X - outer(b_hat_mom_trios, lambda_est_mom_trios)
resid_mom_eigen   <- X - outer(b_hat_mom_eigen, lambda1_hat_eigen)
resid_MGSP        <- X - outer(out_MGSP$b_hat, out_MGSP$lambda1)
resid_HS          <- X - outer(out_HS$b_hat,   out_HS$lambda1)
resid_SSL         <- X - outer(out_SSL$b_hat,  out_SSL$lambda1)

# --- Off-diagonal correlations ---
get_offdiag <- function(M) M[lower.tri(M)]
r_obs    <- get_offdiag(cor(X))
r_mom_trios <- get_offdiag(cor(resid_mom_trios))
r_mom_eigen <- get_offdiag(cor(resid_mom_eigen))
r_MGSP   <- get_offdiag(cor(resid_MGSP))
r_HS     <- get_offdiag(cor(resid_HS))
r_SSL    <- get_offdiag(cor(resid_SSL))

# ---- Permutation null ----
set.seed(123)
X_perm  <- apply(X, 2, sample)
r_perm  <- get_offdiag(cor(X_perm))

# common breaks and y-limit
breaks <- seq(-1, 1, length.out = 41)
all_rs <- list(r_perm, r_obs, r_mom_trios, r_mom_eigen, r_MGSP, r_HS, r_SSL)
max_count <- max(sapply(all_rs, function(v) max(table(cut(v, breaks)))))
bin_width <- diff(breaks)[1]

# Plot base histogram (permutation null)
hist(r_perm, breaks = breaks, col = rgb(0.8, 0.8, 0.8, 0.7), border = "white",
     main = "Residual Off-diagonal Correlations", xlab = "Correlation",
     xlim = c(-1, 1), ylim = c(0, max_count + 5))

# Overlay all residuals
hist(r_obs,        breaks = breaks, col = rgb(0.2,0.4,1,0.4), add=TRUE, border="white")
hist(r_mom_trios,  breaks = breaks, col = rgb(1,0.4,0.4,0.5), add=TRUE, border="white")
hist(r_mom_eigen,  breaks = breaks, col = rgb(0.3,0.2,1,0.3), add=TRUE, border="white")
hist(r_MGSP,       breaks = breaks, col = rgb(0.3,1,0.3,0.3), add=TRUE, border="white")
hist(r_HS,         breaks = breaks, col = rgb(1,0.8,0.2,0.3), add=TRUE, border="white")
hist(r_SSL,        breaks = breaks, col = rgb(0.6,0.2,1,0.2), add=TRUE, border="white")

legend("topright",
       legend = c("Perm. Null", "Observed", "MoM_trios", "MoM_eigen", "MGSP", "HS", "SSL"),
       fill   = c(rgb(0.8,0.8,0.8,0.7),
                  rgb(0.2,0.4,1,0.4),
                  rgb(1,0.4,0.4,0.5),
                  rgb(0.3,0.2,1,0.3),
                  rgb(0.3,1,0.3,0.3),
                  rgb(1,0.8,0.2,0.3),
                  rgb(0.6,0.2,1,0.2)),
       border = NA
)

# --- Density plot of off-diagonal correlations ---
plot(density(r_obs), main = "Density: Off-diag Correlations", xlim = c(-1,1), lwd = 2)
lines(density(r_mom_trios), col = "red",    lwd = 2)
lines(density(r_mom_eigen), col = "blue",   lwd = 2)
lines(density(r_MGSP),      col = "green",  lwd = 2)
lines(density(r_HS),        col = "orange", lwd = 2)
lines(density(r_SSL),       col = "purple", lwd = 2)
legend(
   "topright",
   legend = c("Observed", "MoM_trios", "MoM_eigen", "MGSP", "HS", "SSL"),
   col    = c("black", "red", "blue", "green", "orange", "purple"),
   lwd    = 2
)

# --- Histograms for factor scores (all methods) ---
b_hat_all <- list(MoM_trios = b_hat_mom_trios,
                  MoM_eigen = b_hat_mom_eigen,
                  MGSP      = out_MGSP$b_hat,
                  HS        = out_HS$b_hat,
                  SSL       = out_SSL$b_hat)
for (name in names(b_hat_all)) {
   hist(b_hat_all[[name]], breaks = 100, main = paste(name, "factor scores"), xlab = "b_i")
   qqnorm(b_hat_all[[name]], main = paste("QQ:", name, "b_i")); qqline(b_hat_all[[name]], col = 2)
}

# --- Residual histograms for epsilon_{i1} (all methods) ---
par(mfrow = c(3,2))
hist(epsi(1, X, lambda_est_mom_trios, b_hat_mom_trios), breaks = 50, main = "MoM_trios: epsilon_{i1}", xlab = "Residual")
hist(epsi(1, X, lambda1_hat_eigen, b_hat_mom_eigen),    breaks = 50, main = "MoM_eigen: epsilon_{i1}", xlab = "Residual")
hist(epsi(1, X, out_MGSP$lambda1, out_MGSP$b_hat),      breaks = 50, main = "MGSP: epsilon_{i1}", xlab = "Residual")
hist(epsi(1, X, out_HS$lambda1, out_HS$b_hat),          breaks = 50, main = "HS: epsilon_{i1}", xlab = "Residual")
hist(epsi(1, X, out_SSL$lambda1, out_SSL$b_hat),        breaks = 50, main = "SSL: epsilon_{i1}", xlab = "Residual")
par(mfrow = c(1,1))

# --- QQ Plots for off-diagonal correlations (all methods) ---
par(mfrow = c(2, 4))
qqnorm(r_obs,        main = "QQ: Observed",    pch = 20); qqline(r_obs,        col = "black")
qqnorm(r_mom_trios,  main = "QQ: MoM_trios",   pch = 20); qqline(r_mom_trios,  col = "red")
qqnorm(r_mom_eigen,  main = "QQ: MoM_eigen",   pch = 20); qqline(r_mom_eigen,  col = "blue")
qqnorm(r_MGSP,       main = "QQ: MGSP",        pch = 20); qqline(r_MGSP,       col = "green")
qqnorm(r_HS,         main = "QQ: HS",          pch = 20); qqline(r_HS,         col = "orange")
qqnorm(r_SSL,        main = "QQ: SSL",         pch = 20); qqline(r_SSL,        col = "purple")
qqnorm(r_perm,       main = "QQ: Perm Null",   pch = 20); qqline(r_perm,       col = "grey")
par(mfrow = c(1,1))












# --- Focus only on "non-null" (significant) off-diagonal correlations ---

# 1. Get 95% permutation null interval
ci_null <- quantile(r_perm, c(0.025, 0.975))

# 2. Function to keep only values outside the null
outside_null <- function(r, ci) {
   r[(r < ci[1]) | (r > ci[2])]
}

# 3. Apply to all methods
r_obs_out    <- outside_null(r_obs,    ci_null)
r_mom_trios_out  <- outside_null(r_mom_trios,  ci_null)
r_mom_eigen_out  <- outside_null(r_mom_eigen,  ci_null)
r_MGSP_out   <- outside_null(r_MGSP,   ci_null)
r_HS_out     <- outside_null(r_HS,     ci_null)
r_SSL_out    <- outside_null(r_SSL,    ci_null)

# 4. Plot: Histogram for each method
par(mfrow=c(2,3))
hist(r_obs_out, breaks=30, col=rgb(0.2,0.4,1,0.5), main="Observed (outside null)", xlab="Correlation")
hist(r_mom_trios_out, breaks=30, col=rgb(1,0.4,0.4,0.5), main="MoM trios (outside null)", xlab="Correlation")
hist(r_mom_eigen_out, breaks=30, col=rgb(0.3,0.2,1,0.3), main="MoM eigen (outside null)", xlab="Correlation")
hist(r_MGSP_out, breaks=30, col=rgb(0.3,1,0.3,0.5), main="MGSP (outside null)", xlab="Correlation")
hist(r_HS_out, breaks=30, col=rgb(1,0.8,0.2,0.5), main="HS (outside null)", xlab="Correlation")
hist(r_SSL_out, breaks=30, col=rgb(0.6,0.2,1,0.3), main="SSL (outside null)", xlab="Correlation")
par(mfrow=c(1,1))

# 5. Plot: Densities all on one plot
plot(density(r_obs_out), main="Density: Off-diag (outside null)", xlim=c(-1,1), lwd=2, col="blue")
lines(density(r_mom_trios_out), col="red", lwd=2)
lines(density(r_mom_eigen_out), col="purple", lwd=2)
lines(density(r_MGSP_out), col="green", lwd=2)
lines(density(r_HS_out), col="orange", lwd=2)
lines(density(r_SSL_out), col="darkviolet", lwd=2)
legend("topright",
       legend=c("Observed","MoM_trios","MoM_eigen","MGSP","HS","SSL"),
       col=c("blue","red","purple","green","orange","darkviolet"),
       lwd=2)

# 6. Optionally: print counts
cat("Fraction outside null:\n")
cat(sprintf("Observed:   %d / %d (%.1f%%)\n", length(r_obs_out), length(r_obs), 100*length(r_obs_out)/length(r_obs)))
cat(sprintf("MoM_trios:  %d / %d (%.1f%%)\n", length(r_mom_trios_out), length(r_mom_trios), 100*length(r_mom_trios_out)/length(r_mom_trios)))
cat(sprintf("MoM_eigen:  %d / %d (%.1f%%)\n", length(r_mom_eigen_out), length(r_mom_eigen), 100*length(r_mom_eigen_out)/length(r_mom_eigen)))
cat(sprintf("MGSP:       %d / %d (%.1f%%)\n", length(r_MGSP_out), length(r_MGSP), 100*length(r_MGSP_out)/length(r_MGSP)))
cat(sprintf("HS:         %d / %d (%.1f%%)\n", length(r_HS_out), length(r_HS), 100*length(r_HS_out)/length(r_HS)))
cat(sprintf("SSL:        %d / %d (%.1f%%)\n", length(r_SSL_out), length(r_SSL), 100*length(r_SSL_out)/length(r_SSL)))





