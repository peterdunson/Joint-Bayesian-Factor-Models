# D(F, G): Distance between two distributions based on mean and SD
dsc <- function(corrs, mu2 = 0, sd2 = NULL, n = NULL) {
   mu1 <- mean(corrs)
   sd1 <- sd(corrs)
   if (is.null(sd2)) {
      if (is.null(n)) stop("Need to provide n (sample size) if sd2 is not given")
      sd2 <- 1 / sqrt(n - 3)
   }
   D <- sqrt((mu1 - mu2)^2 + (sd1 - sd2)^2)
   return(list(D = D, mu1 = mu1, sd1 = sd1, mu2 = mu2, sd2 = sd2))
}



dsc_skew <- function(corrs, mu2 = 0, sd2 = NULL, sk2 = 0, n = NULL) {
   mu1 <- mean(corrs)
   sd1 <- sd(corrs)
   sk1 <- mean((corrs - mu1)^3) / sd1^3
   if (is.null(sd2)) {
      if (is.null(n)) stop("Need n (sample size) if sd2 is not given")
      sd2 <- 1 / sqrt(n - 3)
   }
   D <- sqrt((mu1 - mu2)^2 + (sd1 - sd2)^2 + (sign(sk1) * abs(sk1)^(1/3) - sign(sk2) * abs(sk2)^(1/3))^2)
   return(list(DSC = D, mu1 = mu1, sd1 = sd1, sk1 = sk1, mu2 = mu2, sd2 = sd2, sk2 = sk2))
}


dsc_skew_kurt <- function(corrs, mu2 = 0, sd2 = NULL, sk2 = 0, ku2 = 0, n = NULL) {
   mu1 <- mean(corrs)
   sd1 <- sd(corrs)
   # Avoid division by zero in skewness and kurtosis if sd1 is tiny
   if (sd1 < 1e-12) {
      sk1 <- 0
      ku1 <- 0
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
   return(list(
      DSC = D, mu1 = mu1, sd1 = sd1, sk1 = sk1, ku1 = ku1,
      mu2 = mu2, sd2 = sd2, sk2 = sk2, ku2 = ku2
   ))
}