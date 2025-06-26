# ---------------------------------------------------------
# Empirical Null Distribution for DSC
# Simulate uncorrelated (noise-only) data, compute DSC
# ---------------------------------------------------------

set.seed(1001)
n <- 1000       # number of samples
p <- 20         # number of variables
B <- 1000       # number of Monte Carlo null draws

fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
dsc <- function(corrs, mu2 = 0, sd2 = NULL, sk2 = 0, ku2 = 3, n = NULL) {
   mu1 <- mean(corrs)
   sd1 <- sd(corrs)
   if (sd1 < 1e-12) {
      sk1 <- 0; ku1 <- 0
   } else {
      sk1 <- mean((corrs - mu1)^3) / sd1^3
      ku1 <- mean((corrs - mu1)^4) / sd1^4
   }
   if (is.null(sd2)) {
      if (is.null(n)) stop("Need n (sample size) if sd2 is not given")
      sd2 <- 1 / sqrt(n - 3)
   }
   sk_diff <- (sign(sk1) * abs(sk1)^(1/3)) - (sign(sk2) * abs(sk2)^(1/3))
   ku_diff <- (sign(ku1) * abs(ku1)^(1/4)) - (sign(ku2) * abs(ku2)^(1/4))
   D <- sqrt((mu1 - mu2)^2 + (sd1 - sd2)^2 + sk_diff^2 + ku_diff^2)
   return(D)
}

dsc_null <- numeric(B)
for (b in 1:B) {
   Y_null <- matrix(rnorm(n * p), nrow = n, ncol = p)
   R_null <- cor(Y_null)
   z_null <- fisher_z(R_null[lower.tri(R_null)])
   dsc_null[b] <- dsc(z_null, n = n)
}

hist(dsc_null, breaks = 40, main = "Empirical Null Distribution of DSC",
     xlab = "DSC (noise-only)", col = "gray", border = "white")
abline(v = quantile(dsc_null, c(0.025, 0.975)), col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("95%"), col = "red", lwd = 2, lty = 2)

cat(sprintf("Empirical null DSC: Mean = %.3f, SD = %.3f\n", mean(dsc_null), sd(dsc_null)))
cat(sprintf("95%% null interval: %.3f to %.3f\n",
            quantile(dsc_null, 0.025), quantile(dsc_null, 0.975)))
