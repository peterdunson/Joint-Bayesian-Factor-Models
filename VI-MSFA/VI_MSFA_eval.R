# ---------------------------------------------------------------
# Evaluate VI-MSFA (CAVI/FA) Fit on Simulated Data
# ---------------------------------------------------------------

library(VIMSFA)      # For cavi_fa, svi_fa
library(pheatmap)      # For heatmap visualization
library(MatrixCorrelation) # For RV coefficient
library(tidyverse)

# --- USER SETTINGS ---
scenario <- 1
n_train  <- 1000
K        <- 5
sim_dir  <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations"
sim_file <- sprintf("%s/sim_scen%d_%d.rds", sim_dir, scenario, n_train)

# --- Load simulation data ---
if (!file.exists(sim_file)) stop("Simulation file not found: ", sim_file)
sim <- readRDS(sim_file)

yX <- sim$Y      # n x (p+1): first column is y
Lambda_true <- sim$Lambda # (p+1) x K

X <- yX[, -1, drop=FALSE]   # predictors
y <- yX[, 1]                # outcome

# --- Standardize predictors (VI-MSFA default is scale=FALSE, so match that)
X <- scale(X, center=TRUE, scale=FALSE)

# --- Fit VI-MSFA (CAVI and SVI) ---
cavi_est <- cavi_fa(X, K, scale=FALSE)
svi_est  <- svi_fa(X, K, scale=FALSE)

# --- Extract estimated factor loadings and specific variances ---
cavi_lambda <- cavi_est$mean_lambda  # [p, K]
cavi_psi    <- cavi_est$mean_psi    # [p]
svi_lambda  <- svi_est$mean_lambda
svi_psi     <- svi_est$mean_psi

# --- Visualize estimated loadings (CAVI example) ---
pheatmap(cavi_lambda,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Posterior Mean Loadings (CAVI VI-MSFA)",
         color = colorRampPalette(c("blue", "white", "red"))(50))

hist(as.numeric(cavi_lambda), breaks = 50, main = "Histogram of CAVI Loadings")
sparse_frac <- mean(abs(cavi_lambda) < 0.05)
cat("Fraction of near-zero loadings (< 0.05):", round(sparse_frac, 3), "\n")

# --- Compare estimated covariance to true covariance ---
# True Sigma for predictors: sim$Omega[-1, -1]
cavi_sigma <- tcrossprod(cavi_lambda) + diag(cavi_psi)
svi_sigma  <- tcrossprod(svi_lambda) + diag(svi_psi)
true_sigma <- sim$Omega[-1, -1]

cat("RV coefficient (CAVI): ", MatrixCorrelation::RV(cavi_sigma, true_sigma), "\n")
cat("RV coefficient (SVI): ",  MatrixCorrelation::RV(svi_sigma, true_sigma), "\n")

# --- (Optional) Coverage for loadings ---
# If you want CI coverage for loadings, you'll need samples, but VI-MSFA only gives means (no MCMC samples).
# You can check sign/magnitude recovery against truth:
Lambda_true_pred <- Lambda_true[-1, ] # true loadings for predictors only (p x K)
cat("Correlation between true and CAVI-estimated loadings: ",
    cor(as.numeric(Lambda_true_pred), as.numeric(cavi_lambda)), "\n")
cat("Correlation between true and SVI-estimated loadings: ",
    cor(as.numeric(Lambda_true_pred), as.numeric(svi_lambda)), "\n")
