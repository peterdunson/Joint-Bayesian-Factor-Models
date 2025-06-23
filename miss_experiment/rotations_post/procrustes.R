# rotate_with_procrustes.R

# ─── 0) Libraries ─────────────────────────────────────────────────────────────
library(rstan)        # for extract()
library(vegan)        # for procrustes()

# ─── 1) Read in your saved fit object ────────────────────────────────────────
dat <- readRDS("fit_Joint_scen2_scale_all.rds")
post <- dat$posterior
Lambda_samps <- post$Lambda    # array [n_iter, p, K]

# ─── 2) Compute the posterior‐mean loading matrix as your “target” ───────────
Lambda_hat <- apply(Lambda_samps, c(2,3), mean)   # p × K

# ─── 3) For each draw, Procrustes‐align to the target ────────────────────────
n_iter <- dim(Lambda_samps)[1]
p      <- dim(Lambda_samps)[2]
K      <- dim(Lambda_samps)[3]

# Container for rotated draws
Lambda_pr_samps <- array(NA, c(n_iter, p, K))

for (i in seq_len(n_iter)) {
   Ldraw <- Lambda_samps[i, , ]                # p × K
   
   # run Procrustes: align Ldraw onto Lambda_hat
   pr <- procrustes(
      X      = Lambda_hat,  # target
      Y      = Ldraw,       # to be aligned
      scale  = FALSE        # no scaling, only rotation & translation
   )
   
   # pr$Yrot is the rotated version of Y
   Lambda_pr_samps[i, , ] <- pr$Yrot
}

# ─── 4) Summarize your Procrustes‐rotated draws ──────────────────────────────
Lambda_pr_hat   <- apply(Lambda_pr_samps, c(2,3), mean)
ci_lower_pr     <- apply(Lambda_pr_samps, c(2,3), quantile, probs = 0.025)
ci_upper_pr     <- apply(Lambda_pr_samps, c(2,3), quantile, probs = 0.975)

# ─── 5) Save the aligned results ──────────────────────────────────────────────
saveRDS(
   list(
      Lambda_pr_samples = Lambda_pr_samps,
      Lambda_pr_hat     = Lambda_pr_hat,
      ci_lower          = ci_lower_pr,
      ci_upper          = ci_upper_pr
   ),
   file = "fit_Joint_scen2_scale_all_procrustes_rotated.rds"
)

cat("Done: Procrustes‐aligned posterior draws saved to",
    "fit_Joint_scen2_scale_all_procrustes_rotated.rds\n")
