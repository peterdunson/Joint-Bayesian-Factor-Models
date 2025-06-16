set.seed(123)

# Parameters
rep <- 1
n <- 200
N <- n * rep
p <- 30
k <- 5

Lambda <- matrix(0, nrow = p, ncol = k)
numeff <- k + sample(k)  # Each col of Lambda has between k+1 and 2k nonzeros

for (h in 1:k) {
   temp <- sample(1:p, numeff[h])
   Lambda[temp, h] <- rnorm(numeff[h], 0, 1)
}

# True covariance
Ot <- Lambda %*% t(Lambda) + 0.01 * diag(p)

# Generate N x p data matrix, each row ~ N(0, Ot)
library(MASS)
dat <- mvrnorm(n = N, mu = rep(0, p), Sigma = Ot)

# Scale the data
Y <- scale(dat, center = TRUE, scale = TRUE)
n <- nrow(Y)
p <- ncol(Y)
K <- 20   # Upper bound for factors, as in Stan code
