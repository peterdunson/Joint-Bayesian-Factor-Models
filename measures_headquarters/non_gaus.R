# ---------------------------------------------------------
# Empirical Null Distribution for DSC via Permutations
# Each column of Y is permuted independently (breaks dependence, keeps marginals)
# ---------------------------------------------------------

set.seed(123)
n <- nrow(Y)
p <- ncol(Y)
B <- 1000  # Number of permutations

# Fisher-z transform
fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))

# DSC distance function
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

dsc_perm <- numeric(B)
for (b in 1:B) {
   Y_perm <- Y
   for (j in 1:p) {
      Y_perm[, j] <- sample(Y[, j], replace = FALSE)  # randomly permute each column
   }
   R_perm <- cor(Y_perm)
   z_perm <- fisher_z(R_perm[lower.tri(R_perm)])
   dsc_perm[b] <- dsc(z_perm, n = n)
}

hist(dsc_perm, breaks = 40, main = "Permutation Null Distribution of DSC",
     xlab = "DSC (column-wise permutation)", col = "gray", border = "white")
abline(v = quantile(dsc_perm, c(0.025, 0.975)), col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("95% envelope"), col = "red", lwd = 2, lty = 2)

cat(sprintf("Permutation null DSC: Mean = %.3f, SD = %.3f\n", mean(dsc_perm), sd(dsc_perm)))
cat(sprintf("95%% permutation interval: %.3f to %.3f\n",
            quantile(dsc_perm, 0.025), quantile(dsc_perm, 0.975)))
