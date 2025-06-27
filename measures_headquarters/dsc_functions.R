#permutation_null

fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))

# Helper: DSC calculation
dsc <- function(corrs, mu2, sd2, sk2, ku2) {
   mu1 <- mean(corrs)
   sd1 <- sd(corrs)
   sk1 <- if (sd1 < 1e-12) 0 else mean((corrs - mu1)^3) / sd1^3
   ku1 <- if (sd1 < 1e-12) 0 else mean((corrs - mu1)^4) / sd1^4
   sk_diff <- abs(sk1)^(1/3) - abs(sk2)^(1/3)
   ku_diff <- abs(ku1)^(1/4) - abs(ku2)^(1/4)
   sqrt((mu1 - mu2)^2 + (sd1 - sd2)^2 + sk_diff^2 + ku_diff^2)
}

# 1. DSC for observed data, with permutation null
dsc_with_permutation_null_obs <- function(Y, B = 1000) {
   n <- nrow(Y)
   # Observed
   R_obs <- cor(Y)
   z_obs <- fisher_z(R_obs[lower.tri(R_obs)])
   # Permutation null
   dsc_null <- numeric(B)
   mu2s <- sd2s <- sk2s <- ku2s <- numeric(B)
   for (b in 1:B) {
      Y_perm <- apply(Y, 2, sample)
      R_perm <- cor(Y_perm)
      z_perm <- fisher_z(R_perm[lower.tri(R_perm)])
      mu2s[b] <- mean(z_perm)
      sd2s[b] <- sd(z_perm)
      sk2s[b] <- if (sd2s[b] < 1e-12) 0 else mean((z_perm - mu2s[b])^3) / sd2s[b]^3
      ku2s[b] <- if (sd2s[b] < 1e-12) 0 else mean((z_perm - mu2s[b])^4) / sd2s[b]^4
      dsc_null[b] <- dsc(z_obs, mu2s[b], sd2s[b], sk2s[b], ku2s[b])
   }
   # Center null on mean (can also report full null dist)
   dsc_obs <- dsc(z_obs, mean(mu2s), mean(sd2s), mean(sk2s), mean(ku2s))
   list(
      dsc_obs = dsc_obs,
      dsc_null = dsc_null,
      dsc_obs_stats = c(mean = mean(z_obs), sd = sd(z_obs),
                        skew = if (sd(z_obs) < 1e-12) 0 else mean((z_obs - mean(z_obs))^3) / sd(z_obs)^3,
                        kurt = if (sd(z_obs) < 1e-12) 0 else mean((z_obs - mean(z_obs))^4) / sd(z_obs)^4)
   )
}

# 2. DSC for residuals (after removing K factors), with permutation null
dsc_with_permutation_null_resid <- function(Y, Lambda_hat, B = 1000) {
   n <- nrow(Y)
   # Project out factors
   eta_hat <- Y %*% Lambda_hat %*% solve(t(Lambda_hat) %*% Lambda_hat)
   Y_hat <- eta_hat %*% t(Lambda_hat)
   resid <- Y - Y_hat
   R_resid <- cor(resid)
   z_resid <- fisher_z(R_resid[lower.tri(R_resid)])
   # Permutation null
   dsc_null <- numeric(B)
   mu2s <- sd2s <- sk2s <- ku2s <- numeric(B)
   for (b in 1:B) {
      resid_perm <- apply(resid, 2, sample)
      R_perm <- cor(resid_perm)
      z_perm <- fisher_z(R_perm[lower.tri(R_perm)])
      mu2s[b] <- mean(z_perm)
      sd2s[b] <- sd(z_perm)
      sk2s[b] <- if (sd2s[b] < 1e-12) 0 else mean((z_perm - mu2s[b])^3) / sd2s[b]^3
      ku2s[b] <- if (sd2s[b] < 1e-12) 0 else mean((z_perm - mu2s[b])^4) / sd2s[b]^4
      dsc_null[b] <- dsc(z_resid, mu2s[b], sd2s[b], sk2s[b], ku2s[b])
   }
   dsc_resid <- dsc(z_resid, mean(mu2s), mean(sd2s), mean(sk2s), mean(ku2s))
   list(
      dsc_resid = dsc_resid,
      dsc_null = dsc_null,
      dsc_resid_stats = c(mean = mean(z_resid), sd = sd(z_resid),
                          skew = if (sd(z_resid) < 1e-12) 0 else mean((z_resid - mean(z_resid))^3) / sd(z_resid)^3,
                          kurt = if (sd(z_resid) < 1e-12) 0 else mean((z_resid - mean(z_resid))^4) / sd(z_resid)^4)
   )
}
