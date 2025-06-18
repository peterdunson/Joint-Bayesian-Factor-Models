# --------------------------------------------------------------
# Fit VI-MSFA (CAVI/SVI) Bayesian Factor Analysis on Sim Data
# Using blhansen/VI-MSFA package
# --------------------------------------------------------------

# Install VI-MSFA package if needed
# devtools::install_github("blhansen/VI-MSFA")

library(VIMSFA)

# Load for use
library(MatrixCorrelation)

# ---- Load Simulation Data ----
sim_path <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen1_1000.rds"
sim <- readRDS(sim_path)

# sim$Y: [n x (p+1)] matrix, first column is y (outcome), rest are predictors
X <- sim$Y[, -1, drop=FALSE]  # Drop first column (outcome), keep predictors only

n <- nrow(X)
p <- ncol(X)

# ---- Set number of factors (K) ----
J <- 5    # (should match K in simulation: see your K <- 5)

# ---- Fit CAVI ----
cavi_est <- cavi_fa(X, J, scale = FALSE)  # Do NOT rescale (already standardized if needed)
# Extract estimated parameters
cavi_lambda <- cavi_est$mean_lambda   # [p x J] posterior mean factor loadings
cavi_psi    <- cavi_est$mean_psi      # [p] posterior mean residual variances
cavi_sigma  <- tcrossprod(cavi_lambda) + diag(cavi_psi)   # Implied covariance

# ---- Fit SVI (optional) ----
svi_est <- svi_fa(X, J, scale = FALSE)
svi_lambda <- svi_est$mean_lambda
svi_psi <- svi_est$mean_psi
svi_sigma <- tcrossprod(svi_lambda) + diag(svi_psi)

# ---- Compare true to estimated covariance (optional) ----
if ("Omega" %in% names(sim)) {
   # Remove the first row and column (corresponding to outcome y)
   true_sigma <- sim$Omega[-1, -1]
   # Now both are p x p matrices
   print(MatrixCorrelation::RV(cavi_sigma, true_sigma))
   print(MatrixCorrelation::RV(svi_sigma, true_sigma))
   
}

# ---- Save fit if desired ----
saveRDS(cavi_est, file = "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/vimsfa_cavi_fit_scen1_5.rds")
saveRDS(svi_est,  file = "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/vimsfa_svi_fit_scen1_5.rds")


