# ---- 1. Center and scale if not already done ----
Y <- scale(dat, center = TRUE, scale = TRUE)    # n x p

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/applications")

# ---- 2. Load fit and get Lambda_hat ----
fit_obj <- readRDS("fit_Joint_NHANES1718.rds")
Lambda_hat <- fit_obj$Lambda_hat
if (nrow(Lambda_hat) != ncol(Y)) Lambda_hat <- t(Lambda_hat)

# ---- 3. DSC for observed ----
dsc_obs_result <- dsc_with_permutation_null_obs(Y, B = 500)   # B = 500 or higher for stability

cat("DSC (Observed): ", dsc_obs_result$dsc_obs, "\n")
cat("Permutation null mean: ", mean(dsc_obs_result$dsc_null), "\n")
cat("Permutation null 95% range: ", quantile(dsc_obs_result$dsc_null, c(0.025, 0.975)), "\n")

# ---- 4. DSC for residuals ----
dsc_resid_result <- dsc_with_permutation_null_resid(Y, Lambda_hat, B = 500)

cat("DSC (Residuals): ", dsc_resid_result$dsc_resid, "\n")
cat("Permutation null mean: ", mean(dsc_resid_result$dsc_null), "\n")
cat("Permutation null 95% range: ", quantile(dsc_resid_result$dsc_null, c(0.025, 0.975)), "\n")

# ---- 5. Optional: Plot DSC null distributions ----
par(mfrow = c(1,2))
hist(dsc_obs_result$dsc_null, breaks = 30, main = "DSC Null (Observed)", xlab = "DSC", col = "#AAAAFF")
abline(v = dsc_obs_result$dsc_obs, col = "red", lwd = 2)
hist(dsc_resid_result$dsc_null, breaks = 30, main = "DSC Null (Residuals)", xlab = "DSC", col = "#FFAAAA")
abline(v = dsc_resid_result$dsc_resid, col = "red", lwd = 2)
par(mfrow = c(1,1))
