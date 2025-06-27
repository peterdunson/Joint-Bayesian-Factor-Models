# ─── Paths ────────────────────────────────────────────────────────────────────
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment")
sim_path2 <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds"
fit_path2 <- "fit_Joint_scen2_scale_all_randominit.rds"

# ─── 1) Load data ─────────────────────────────────────────────────────────────
sim2      <- readRDS(sim_path2)
Y2        <- sim2$Y                 # your n×p data matrix

fit2      <- readRDS(fit_path2)
Lambda2   <- fit2$Lambda_hat        # or whatever your object names it

# ─── 2) Compute DSC on observed correlations under permutation null ───────────
set.seed(42)                        # for reproducibility of permutations
res_obs   <- dsc_with_permutation_null_obs(Y = Y2, B = 10000)

# Inspect
res_obs$dsc_obs                     # single DSC distance for observed vs null
summary(res_obs$dsc_null)           # null distribution of DSC distances
res_obs$dsc_obs_stats               # mean, sd, skew, kurt of observed Fisher‐Z’s

# ─── 3) Compute DSC on residuals under permutation null ──────────────────────
set.seed(42)
res_resid <- dsc_with_permutation_null_resid(Y = Y2, Lambda_hat = Lambda2, B = 10000)

# Inspect
res_resid$dsc_resid                 # DSC for residual structure vs null
summary(res_resid$dsc_null)         # null distribution for residual DSC
res_resid$dsc_resid_stats   




hist(res_resid$dsc_null,
     breaks = 100,
     main   = "Null Distribution of Residual DSC",
     xlab   = "DSC under Permutation Null",
     col    = "lightgray",
     border = "white")

# 2) Add observed DSC
abline(v = res_resid$dsc_resid,
       col = "red",
       lwd = 2)

# 3) Annotate
legend("topright",
       legend = sprintf("Observed DSC = %.3f", res_resid$dsc_resid),
       col    = "red",
       lwd    = 2,
       bty    = "n")





hist(res_obs$dsc_null,
     breaks = 200,
     main   = "Null Distribution of Observed Correlation DSC",
     xlab   = "DSC under Permutation Null",
     col    = "lightgray",
     border = "white")

# 2) Add your observed DSC
abline(v = res_obs$dsc_obs,
       col = "blue",
       lwd = 2)

# 3) Annotate
legend("topright",
       legend = sprintf("Observed DSC = %.3f", res_obs$dsc_obs),
       col    = "blue",
       lwd    = 2,
       bty    = "n")









# 1) Compute residuals
eta_hat <- Y2 %*% Lambda2 %*% solve(t(Lambda2) %*% Lambda2)
resid   <- Y2 - (eta_hat %*% t(Lambda2))
n       <- nrow(resid)

# permutation null of z’s
B       <- 10000
z_perm  <- replicate(B, {
   Rp   <- apply(resid, 2, sample)
   cor(Rp)[lower.tri(cor(Rp))]
}) %>% as.vector() %>% fisher_z()

# 3) Plot
hist(z_perm,
     prob    = TRUE,
     breaks  = 200,
     main    = "Permutation vs. Theoretical Null of Fisher-Z Correlations",
     xlab    = "Fisher-Z of r",
     col     = "lightgray",
     border  = "white")

# 4) Overlay N(0, 1/(n−3))
curve(dnorm(x, mean = 0, sd = 1/sqrt(n-3)),
      col   = "blue",
      lwd   = 2,
      add   = TRUE)

legend("topright",
       legend = c("Perm. null", sprintf("Theoretical N(0,1/(n-3))")),
       fill   = c("lightgray", NA),
       border = c(NA, NA),
       lty    = c(NA,1),
       col    = c(NA,"blue"),
       bty    = "n")


