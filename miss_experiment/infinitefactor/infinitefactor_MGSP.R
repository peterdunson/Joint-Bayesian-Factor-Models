# fit_tebfar_fixed_sigma.R

# ─── 1) Load Dependencies ─────────────────────────────────────────────────
library(Rcpp)
library(mvtnorm)       # for multivariate normals if needed
# no need to load infinitefactor since we source its R files directly


# ─── 2) Source the edited infinitefactor routines ─────────────────────────
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
cat("Estimated K:", fit_TEBFAR$numFactors, "\n")

# Posterior-mean induced covariance Ω = ΛΛᵀ + Σ
Omega_hat <- apply(fit_TEBFAR$omegaSamps, c(1,2), mean)
print(dim(Omega_hat))   # should be (p+1) × (p+1)

# 1) See what’s in the result
str(fit_TEBFAR)

# 1) parameters
p <- nrow(sim$Lambda)
K <- ncol(sim$Lambda)      

# 2) build a p×K matrix of medians
Lambda_mean <- matrix(NA_real_, nrow = p, ncol = K)
for (j in 1:p) {
   for (k in 1:K) {
      # collect the (j,k) loading across iterations *only* if that draw had >= k cols
      vals <- vapply(fit_TEBFAR$lambdaSamps, 
                     FUN = function(mat) if (ncol(mat) >= k) mat[j,k] else NA_real_,
                     FUN.VALUE = 0.0)
      Lambda_mean[j,k] <- median(vals, na.rm = TRUE)
   }
}

library(pheatmap)

# 1) Define your custom palette
custom_cols <- colorRampPalette(c("pink","white","grey4"))(100)

# 2) Re‐use the same break logic we had before
maxval <- max(abs(Lambda_mean))
brks   <- seq(-maxval, maxval, length.out = length(custom_cols) + 1)

# 3) Draw the heatmap with your palette
pheatmap(
   Lambda_mean,
   color        = custom_cols,
   breaks       = brks,
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   main         = "Posterior-Median Loadings\n(custom colors)"
)

