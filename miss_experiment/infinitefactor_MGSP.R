# fit_tebfar_fixed_sigma.R

# ─── 1) Load Dependencies ─────────────────────────────────────────────────
library(Rcpp)
library(mvtnorm)       # for multivariate normals if needed
# no need to load infinitefactor since we source its R files directly

# ─── 2) Source the edited infinitefactor routines ─────────────────────────
# adjust these paths to wherever you unzipped infinitefactor_edits
source("~/Desktop/infinitefactor_edits/R/linearMGSP_fixedSigma1.R")
sourceCpp("~/Desktop/infinitefactor_edits/src/samplerBits.cpp", rebuild = TRUE)

# ─── 3) Load & prepare your Scenario-2 data ───────────────────────────────
sim    <- readRDS("~/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Y      <- sim$Y               # n × (p+1), first column is outcome y
yX     <- scale(Y, center = TRUE, scale = TRUE)

# ─── 4) Set sampler parameters ────────────────────────────────────────────
sigma2_y <- 1.0   # fixed Var(y) after scaling
nrun      <- 5000
burn      <- 2500

# ─── 5) Fit TEB-FAR style model with fixed σ²ᵧ ─────────────────────────────
fit_TEBFAR <- linearMGSP_fixedSigma1(
  X      = yX,
  Sigma1 = sigma2_y,
  nrun   = nrun,
  burn   = burn
)

# ─── 6) Quick look at results ─────────────────────────────────────────────
# Number of factors retained (posterior mode)
cat("Estimated K:", fit_TEBFAR$numFactors, "\n")

# Posterior-mean induced covariance Ω = ΛΛᵀ + Σ
Omega_hat <- apply(fit_TEBFAR$omegaSamps, c(1,2), mean)
print(dim(Omega_hat))   # should be (p+1) × (p+1)

# Example: heatmap of the posterior-mean loadings
library(pheatmap)
Lambda_mean <- apply(fit_TEBFAR$lambdaSamps, c(1,2), mean)
pheatmap(
  Lambda_mean,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Posterior-Mean Loadings (Fixed σ²ᵧ)"
)
