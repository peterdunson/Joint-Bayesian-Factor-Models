

B <- 500 # Number of permutation draws (increase for final results)

# 1. Compute for original (before, same for all)
dsc_before <- dsc_with_permutation_null_obs(Y, B = B)
dsc_orig     <- dsc_before$dsc_obs
null_mean    <- mean(dsc_before$dsc_null)
null_sd      <- sd(dsc_before$dsc_null)

# 2. Compute for residuals after K=1 removal (each method)
# (Each will have its own null distribution and SD)
dsc_trio    <- dsc_with_permutation_null_obs(Y_resid_trio, B = B)
dsc_mgsp    <- dsc_with_permutation_null_obs(Y_resid_mgsp, B = B)
dsc_hs      <- dsc_with_permutation_null_obs(Y_resid_hs, B = B)
dsc_ssl     <- dsc_with_permutation_null_obs(Y_resid_ssl, B = B)

# Collect observed DSCs and null SDs
method_names <- c("Original Data", "MoM-Trio Residual", "MGSP Residual", "Horseshoe Residual", "Spike-and-Slab Residual")
dsc_obs      <- c(dsc_orig, dsc_trio$dsc_obs, dsc_mgsp$dsc_obs, dsc_hs$dsc_obs, dsc_ssl$dsc_obs)
dsc_null_means <- c(null_mean, mean(dsc_trio$dsc_null), mean(dsc_mgsp$dsc_null), mean(dsc_hs$dsc_null), mean(dsc_ssl$dsc_null))
dsc_null_sds   <- c(null_sd, sd(dsc_trio$dsc_null), sd(dsc_mgsp$dsc_null), sd(dsc_hs$dsc_null), sd(dsc_ssl$dsc_null))

# Print table for easy copy-paste
dsc_table <- data.frame(
   Method      = method_names,
   Observed_DSC = round(dsc_obs, 3),
   Null_Mean   = round(dsc_null_means, 3),
   Null_SD     = round(dsc_null_sds, 3)
)
print(dsc_table)

# ---- BOXPLOT: DSC null distributions ----

# Collect null DSCs into a list for boxplot
dsc_nulls <- list(
   "Original Data"      = dsc_before$dsc_null,
   "MoM-Trio Residual"  = dsc_trio$dsc_null,
   "MGSP Residual"      = dsc_mgsp$dsc_null,
   "Horseshoe Residual" = dsc_hs$dsc_null,
   "Spike-and-Slab Residual" = dsc_ssl$dsc_null
)

par(mar = c(7, 4, 3, 1))
boxplot(
   dsc_nulls,
   col = c("gray50", "firebrick", "royalblue", "darkgreen", "goldenrod"),
   ylab = "DSC Distance (Permutation Null)",
   main = "Permutation Null DSC Distributions (Before/Residuals)",
   las = 2,
   outline = FALSE
)
# Overlay observed DSC as big points
points(1:5, dsc_obs, pch = 19, cex = 1.7, col = "black")
# Optional: Overlay null mean as small dots
points(1:5, dsc_null_means, pch = 21, cex = 1.1, bg = "white", col = "blue")
legend("topright", legend = c("Observed DSC", "Null Mean"), pch = c(19, 21), col = c("black", "blue"), pt.bg = c(NA, "white"), bty = "n")

