# rotate_posteriors.R

# ─── 0) Libraries ─────────────────────────────────────────────────────────────
library(rstan)        # for extract()
library(GPArotation)  # for varimax()
# (no ggplot here—this is purely for rotation + summaries)

# ─── 1) Read in your saved fit object ────────────────────────────────────────
# This file should contain a list with at least a component named "posterior",
# where posterior$Lambda is an array [n_iter × p × K].
dat <- readRDS("fit_Joint_scen2_scale_all.rds")

post <- dat$posterior
Lambda_samps <- post$Lambda    # array of dim [n_iter, p, K]

# ─── 2) Rotate each draw ─────────────────────────────────────────────────────
n_iter <- dim(Lambda_samps)[1]
p      <- dim(Lambda_samps)[2]
K      <- dim(Lambda_samps)[3]

# Pre-allocate rotated draws
Lambda_vm_samps <- array(NA, c(n_iter, p, K))

for (i in seq_len(n_iter)) {
   Ldraw <- Lambda_samps[i, , ]            # p × K matrix
   vm    <- varimax(Ldraw)                 # run varimax
   Lambda_vm_samps[i, , ] <- vm$loadings   # store rotated p×K
}

# ─── 3) Summarize rotated draws ──────────────────────────────────────────────
# Posterior mean of rotated loadings:
Lambda_vm_hat <- apply(Lambda_vm_samps, c(2,3), mean)

# 95% credible intervals for each loading:
ci_lower <- apply(Lambda_vm_samps, c(2,3), quantile, probs = 0.025)
ci_upper <- apply(Lambda_vm_samps, c(2,3), quantile, probs = 0.975)

# ─── 4) Save results ─────────────────────────────────────────────────────────
saveRDS(
   list(
      Lambda_vm_samples = Lambda_vm_samps,
      Lambda_vm_hat     = Lambda_vm_hat,
      ci_lower          = ci_lower,
      ci_upper          = ci_upper
   ),
   file = "fit_Joint_scen2_scale_all_rotated_fullposterior.rds"
)

cat("Done: rotated posterior draws saved to",
    "fit_Joint_scen2_scale_all_rotated_fullposterior.rds\n")
