# ------------------------------------------------------------
# DSC Analysis for Three Fitted Joint Factor Models
# (Scenarios 1, 2, and 3)
# ------------------------------------------------------------

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment")

# -- Helper functions --
fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
dsc <- function(corrs, mu2 = 0, sd2 = NULL, sk2 = 0, ku2 = 3, n = NULL) {
   mu1 <- mean(corrs)
   sd1 <- sd(corrs)
   if (sd1 < 1e-12) {
      sk1 <- 0; ku1 <- 0
   } else {
      sk1 <- mean((corrs - mu1)^3) / sd1^3
      ku1 <- mean((corrs - mu1)^4) / sd1^4
   }
   if (is.null(sd2)) {
      if (is.null(n)) stop("Need n (sample size) if sd2 is not given")
      sd2 <- 1 / sqrt(n - 3)
   }
   sk_diff <- (sign(sk1) * abs(sk1)^(1/3)) - (sign(sk2) * abs(sk2)^(1/3))
   ku_diff <- (sign(ku1) * abs(ku1)^(1/4)) - (sign(ku2) * abs(ku2)^(1/4))
   D <- sqrt((mu1 - mu2)^2 + (sd1 - sd2)^2 + sk_diff^2 + ku_diff^2)
   return(list(DSC = D, mu1 = mu1, sd1 = sd1, sk1 = sk1, ku1 = ku1))
}

# ----- Paths -----
fit_path1 <- "fit_Joint_scen1_scale_all_randominit_1.rds"
fit_path2 <- "fit_Joint_scen2_scale_all_randominit.rds"
fit_path3 <- "fit_Joint_scen3_scale_all_randominit_3.rds"

sim_path1 <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen1_1000.rds"
sim_path2 <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds"
sim_path3 <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen3_1000.rds"

# ---------- SCENARIO 1 ----------
fit1 <- readRDS(fit_path1)
sim1 <- readRDS(sim_path1)
Y1 <- scale(sim1$Y, center = TRUE, scale = TRUE)

Lambda_hat1 <- fit1$Lambda_hat
n1 <- nrow(Y1)

eta_hat1 <- Y1 %*% Lambda_hat1 %*% solve(t(Lambda_hat1) %*% Lambda_hat1)
Y_hat1 <- eta_hat1 %*% t(Lambda_hat1)
resid1 <- Y1 - Y_hat1

R_obs1 <- cor(Y1)
R_resid1 <- cor(resid1)

z_obs1 <- fisher_z(R_obs1[lower.tri(R_obs1)])
z_resid1 <- fisher_z(R_resid1[lower.tri(R_resid1)])

dsc_obs1 <- dsc(z_obs1, n=n1)
dsc_resid1 <- dsc(z_resid1, n=n1)

cat("=== Results for SCENARIO 1 ===\n")
cat(sprintf(
   "DSC obs: %.3f (mean=%.3f, sd=%.3f, skew=%.3f, kurtosis=%.3f)\n",
   dsc_obs1$DSC, dsc_obs1$mu1, dsc_obs1$sd1, dsc_obs1$sk1, dsc_obs1$ku1))
cat(sprintf(
   "DSC resid: %.3f (mean=%.3f, sd=%.3f, skew=%.3f, kurtosis=%.3f)\n\n",
   dsc_resid1$DSC, dsc_resid1$mu1, dsc_resid1$sd1, dsc_resid1$sk1, dsc_resid1$ku1))

# ---------- SCENARIO 2 ----------
fit2 <- readRDS(fit_path2)
sim2 <- readRDS(sim_path2)
Y2 <- scale(sim2$Y, center = TRUE, scale = TRUE)

Lambda_hat2 <- fit2$Lambda_hat
n2 <- nrow(Y2)

eta_hat2 <- Y2 %*% Lambda_hat2 %*% solve(t(Lambda_hat2) %*% Lambda_hat2)
Y_hat2 <- eta_hat2 %*% t(Lambda_hat2)
resid2 <- Y2 - Y_hat2

R_obs2 <- cor(Y2)
R_resid2 <- cor(resid2)

z_obs2 <- fisher_z(R_obs2[lower.tri(R_obs2)])
z_resid2 <- fisher_z(R_resid2[lower.tri(R_resid2)])

dsc_obs2 <- dsc(z_obs2, n=n2)
dsc_resid2 <- dsc(z_resid2, n=n2)

cat("=== Results for SCENARIO 2 ===\n")
cat(sprintf(
   "DSC obs: %.3f (mean=%.3f, sd=%.3f, skew=%.3f, kurtosis=%.3f)\n",
   dsc_obs2$DSC, dsc_obs2$mu1, dsc_obs2$sd1, dsc_obs2$sk1, dsc_obs2$ku1))
cat(sprintf(
   "DSC resid: %.3f (mean=%.3f, sd=%.3f, skew=%.3f, kurtosis=%.3f)\n\n",
   dsc_resid2$DSC, dsc_resid2$mu1, dsc_resid2$sd1, dsc_resid2$sk1, dsc_resid2$ku1))

# ---------- SCENARIO 3 ----------
fit3 <- readRDS(fit_path3)
sim3 <- readRDS(sim_path3)
Y3 <- scale(sim3$Y, center = TRUE, scale = TRUE)

Lambda_hat3 <- fit3$Lambda_hat
n3 <- nrow(Y3)

eta_hat3 <- Y3 %*% Lambda_hat3 %*% solve(t(Lambda_hat3) %*% Lambda_hat3)
Y_hat3 <- eta_hat3 %*% t(Lambda_hat3)
resid3 <- Y3 - Y_hat3

R_obs3 <- cor(Y3)
R_resid3 <- cor(resid3)

z_obs3 <- fisher_z(R_obs3[lower.tri(R_obs3)])
z_resid3 <- fisher_z(R_resid3[lower.tri(R_resid3)])

dsc_obs3 <- dsc(z_obs3, n=n3)
dsc_resid3 <- dsc(z_resid3, n=n3)

cat("=== Results for SCENARIO 3 ===\n")
cat(sprintf(
   "DSC obs: %.3f (mean=%.3f, sd=%.3f, skew=%.3f, kurtosis=%.3f)\n",
   dsc_obs3$DSC, dsc_obs3$mu1, dsc_obs3$sd1, dsc_obs3$sk1, dsc_obs3$ku1))
cat(sprintf(
   "DSC resid: %.3f (mean=%.3f, sd=%.3f, skew=%.3f, kurtosis=%.3f)\n",
   dsc_resid3$DSC, dsc_resid3$mu1, dsc_resid3$sd1, dsc_resid3$sk1, dsc_resid3$ku1))

