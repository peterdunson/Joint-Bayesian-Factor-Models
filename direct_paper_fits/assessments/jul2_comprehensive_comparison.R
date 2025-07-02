
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
lambda_est_mom <- sign(lambda2_est) * sqrt(abs(lambda2_est))
denom_mom <- sum(lambda_est_mom^2)
b_hat_mom <- as.numeric(X %*% lambda_est_mom) / denom_mom


# -- For all fits: MGSP, HS, SSL --

extract_lambda_b <- function(fit_obj, X) {
   # Lambda_hat is always [p x 1] for k=1
   lambda1 <- as.numeric(fit_obj$Lambda_hat[, 1])  # Ensures vector
   denom   <- sum(lambda1^2)
   b_hat   <- as.numeric(X %*% lambda1) / denom    # n-vector
   list(lambda1 = lambda1, b_hat = b_hat)
}

# Extract for each fit (no need for transpose)
out_MGSP <- extract_lambda_b(fit_MGSP, X)
out_HS   <- extract_lambda_b(fit_HS,   X)
out_SSL  <- extract_lambda_b(fit_SSL,  X)



# Compare MoM b_i to each Bayesian model's b_i, one plot at a time
par(mfrow = c(1, 1))
plot(b_hat_mom, out_MGSP$b_hat, main = "MoM vs MGSP", xlab = "MoM b_i", ylab = "MGSP b_i")
plot(b_hat_mom, out_HS$b_hat,   main = "MoM vs HS",   xlab = "MoM b_i", ylab = "HS b_i")
plot(b_hat_mom, out_SSL$b_hat,  main = "MoM vs SSL",  xlab = "MoM b_i", ylab = "SSL b_i")
plot(out_MGSP$b_hat, out_HS$b_hat, main = "MGSP vs HS", xlab = "MGSP b_i", ylab = "HS b_i")
plot(out_MGSP$b_hat, out_SSL$b_hat, main = "MGSP vs SSL", xlab = "MGSP b_i", ylab = "SSL b_i")
plot(out_HS$b_hat, out_SSL$b_hat, main = "HS vs SSL", xlab = "HS b_i", ylab = "SSL b_i")


# For the first 3 variables, plot X_j vs (lambda1_j * b_i) for each method, each plot separately
for (j in 1:3) {
   plot(X[,j], lambda_est_mom[j] * b_hat_mom,
        main = paste("MoM: X_", j, "vs lambda1 b_i"), xlab = paste("X_", j), ylab = "lambda1 b_i")
   plot(X[,j], out_MGSP$lambda1[j] * out_MGSP$b_hat,
        main = paste("MGSP: X_", j, "vs lambda1 b_i"), xlab = paste("X_", j), ylab = "lambda1 b_i")
   plot(X[,j], out_HS$lambda1[j] * out_HS$b_hat,
        main = paste("HS: X_", j, "vs lambda1 b_i"), xlab = paste("X_", j), ylab = "lambda1 b_i")
   plot(X[,j], out_SSL$lambda1[j] * out_SSL$b_hat,
        main = paste("SSL: X_", j, "vs lambda1 b_i"), xlab = paste("X_", j), ylab = "lambda1 b_i")
}



# --- Function for epsilons (residuals for variable j)
epsi <- function(j, X, lambda1, b_hat) {
   X[,j] - lambda1[j] * b_hat
}

# ---- b_i vs epsilon_{i1} (plot separately) ----
plot(b_hat_mom, epsi(1, X, lambda_est_mom, b_hat_mom), main="MoM: b_i vs epsilon_{i1}", xlab="b_i", ylab="epsilon_{i1}")
plot(out_MGSP$b_hat, epsi(1, X, out_MGSP$lambda1, out_MGSP$b_hat), main="MGSP: b_i vs epsilon_{i1}", xlab="b_i", ylab="epsilon_{i1}")
plot(out_HS$b_hat,   epsi(1, X, out_HS$lambda1, out_HS$b_hat), main="HS: b_i vs epsilon_{i1}", xlab="b_i", ylab="epsilon_{i1}")
plot(out_SSL$b_hat,  epsi(1, X, out_SSL$lambda1, out_SSL$b_hat), main="SSL: b_i vs epsilon_{i1}", xlab="b_i", ylab="epsilon_{i1}")

# ---- Residuals after rank-1 removal ----
resid_mom   <- X - outer(b_hat_mom, lambda_est_mom)
resid_MGSP  <- X - outer(out_MGSP$b_hat, out_MGSP$lambda1)
resid_HS    <- X - outer(out_HS$b_hat,   out_HS$lambda1)
resid_SSL   <- X - outer(out_SSL$b_hat,  out_SSL$lambda1)

get_offdiag <- function(M) M[lower.tri(M)]
r_obs    <- get_offdiag(cor(X))
r_mom    <- get_offdiag(cor(resid_mom))
r_MGSP   <- get_offdiag(cor(resid_MGSP))
r_HS     <- get_offdiag(cor(resid_HS))
r_SSL    <- get_offdiag(cor(resid_SSL))

# ---- Permutation null ----
set.seed(123)
X_perm <- apply(X, 2, sample)
r_perm <- get_offdiag(cor(X_perm))

breaks <- seq(-1, 1, length.out = 41)
hist(r_perm, breaks = breaks, col = rgb(0.8, 0.8, 0.8, 0.7), border = "white",
     main = "Residual Off-diagonal Correlations (Perm Null)", xlab = "Correlation", xlim = c(-1, 1),
     ylim = c(0, max(table(cut(r_obs, breaks)))+10))
hist(r_obs,  breaks = breaks, col = rgb(0.2,0.4,1,0.4), add=TRUE, border="white")
hist(r_mom,  breaks = breaks, col = rgb(1,0.4,0.4,0.5), add=TRUE, border="white")
hist(r_MGSP, breaks = breaks, col = rgb(0.3,1,0.3,0.3), add=TRUE, border="white")
hist(r_HS,   breaks = breaks, col = rgb(1,0.8,0.2,0.3), add=TRUE, border="white")
hist(r_SSL,  breaks = breaks, col = rgb(0.6,0.2,1,0.2), add=TRUE, border="white")
legend("topright",
       legend = c("Perm. Null", "Observed", "MoM", "MGSP", "HS", "SSL"),
       fill   = c(rgb(0.8,0.8,0.8,0.7), rgb(0.2,0.4,1,0.4), rgb(1,0.4,0.4,0.5),
                  rgb(0.3,1,0.3,0.3), rgb(1,0.8,0.2,0.3), rgb(0.6,0.2,1,0.2)),
       border = NA)

# ---- Residual histograms for epsilon_{i1} (all separate) ----
hist(epsi(1, X, lambda_est_mom, b_hat_mom), breaks = 50, main = "MoM: epsilon_{i1}", xlab = "Residual")
hist(epsi(1, X, out_MGSP$lambda1, out_MGSP$b_hat), breaks = 50, main = "MGSP: epsilon_{i1}", xlab = "Residual")
hist(epsi(1, X, out_HS$lambda1, out_HS$b_hat), breaks = 50, main = "HS: epsilon_{i1}", xlab = "Residual")
hist(epsi(1, X, out_SSL$lambda1, out_SSL$b_hat), breaks = 50, main = "SSL: epsilon_{i1}", xlab = "Residual")


# ---- Density plot of off-diagonal correlations: MoM vs Observed ----
plot(density(r_obs), main = "Density of Off-diagonal Correlations", xlim = c(-1,1), lwd = 2)
lines(density(r_mom),   col = "red",    lwd = 2)
lines(density(r_MGSP),  col = "green",  lwd = 2)
lines(density(r_HS),    col = "orange", lwd = 2)
lines(density(r_SSL),   col = "purple", lwd = 2)
legend(
   "topright",
   legend = c("Observed", "MoM", "MGSP", "HS", "SSL"),
   col    = c("black", "red", "green", "orange", "purple"),
   lwd    = 2
)


cor(out_HS$b_hat, out_SSL$b_hat)
summary(out_HS$b_hat - (-out_SSL$b_hat))

# QQ Plots for off-diagonal correlations

par(mfrow = c(2, 3))

qqnorm(r_obs,   main = "QQ: Observed",   pch = 20); qqline(r_obs,   col = "blue")
qqnorm(r_mom,   main = "QQ: MoM",        pch = 20); qqline(r_mom,   col = "red")
qqnorm(r_MGSP,  main = "QQ: MGSP",       pch = 20); qqline(r_MGSP,  col = "green")
qqnorm(r_HS,    main = "QQ: HS",         pch = 20); qqline(r_HS,    col = "orange")
qqnorm(r_SSL,   main = "QQ: SSL",        pch = 20); qqline(r_SSL,   col = "purple")

# Permutation null QQ plot
qqnorm(r_perm,  main = "QQ: Perm Null",  pch = 20); qqline(r_perm,  col = "grey")

par(mfrow = c(1, 1))






