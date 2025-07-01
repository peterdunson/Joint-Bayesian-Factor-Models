# ---- Load functions (if not already loaded) ----
source("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/measures_headquarters/dsc_functions.R")

# ---- Load sim2 original data ----
sim2 <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000_1.rds")
Y2 <- scale(sim2$Y, center = TRUE, scale = TRUE)

# ---- List of fits and model names ----
fit_list <- list(
   MGSP   = fit_MGSP,
   HS     = fit_HS,
   SSL    = fit_SSL,
   MASS   = fit_MASS,
   Robust = fit_RobustMGSP
)

B <- 500  # Number of permutation draws (increase to 1000+ for final runs)

dsc_results <- list()

for (model in names(fit_list)) {
   cat(sprintf("\n--- %s ---\n", model))
   fit <- fit_list[[model]]
   Lambda_hat <- fit$Lambda_hat
   
   # DSC for observed data
   dsc_obs <- dsc_with_permutation_null_obs(Y2, B = B)
   
   # DSC for residuals after removing 1 factor
   dsc_resid <- dsc_with_permutation_null_resid(Y2, Lambda_hat, B = B)
   
   dsc_results[[model]] <- list(obs = dsc_obs, resid = dsc_resid)
   
   # Print summary
   cat(sprintf("DSC obs:   %.3f\n", dsc_obs$dsc_obs))
   cat(sprintf("DSC resid: %.3f\n", dsc_resid$dsc_resid))
}

# ---- Visualization ----
par(mfrow = c(2, length(fit_list)), mar = c(4,4,2,1))

for (model in names(fit_list)) {
   # Observed DSC
   hist(dsc_results[[model]]$obs$dsc_null,
        breaks = 40, main = paste(model, "- OBS"), xlab = "DSC", col = "grey80", border = "white")
   abline(v = dsc_results[[model]]$obs$dsc_obs, col = "red", lwd = 2)
   legend("topright", legend = "Observed", col = "red", lwd = 2, bty = "n")
}

for (model in names(fit_list)) {
   # Residual DSC
   hist(dsc_results[[model]]$resid$dsc_null,
        breaks = 40, main = paste(model, "- RESID"), xlab = "DSC", col = "grey80", border = "white")
   abline(v = dsc_results[[model]]$resid$dsc_resid, col = "blue", lwd = 2)
   legend("topright", legend = "Residual", col = "blue", lwd = 2, bty = "n")
}
par(mfrow = c(1,1))
