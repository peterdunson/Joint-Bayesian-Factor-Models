# simulation_gen.R
# Fixed λₚ for p = 1…10, zeros for p = 11…20; n = 1000, σ² = 0.2

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim")

set.seed(26031980)

# 1) Parameters
P       <- 20
n       <- 1000
sigmasq <- 0.2

# 2) Fixed loadings vector λₚ
lambdaP <- c(rep(1 / sqrt(10), 10), rep(0, 10))  # length P

# 3) Simulate noise εᵢₚ ∼ N(0, σ²)
Epsilon <- sqrt(sigmasq) * matrix(rnorm(n * P), nrow = n, ncol = P)

# 4) Simulate latent scores bᵢ ∼ N(0,1)
b <- rnorm(n)

# 5) Build data Xᵢₚ = λₚ bᵢ + εᵢₚ
X <- outer(b, lambdaP) + Epsilon

C <- cov(X)
hist(C, breaks = 20, main = "Histogram of Covariances", xlab = "C[p,q]")

sim <- list(
   X       = X,        # n × P data matrix
   lambdaP = lambdaP,  # true loadings
   b       = b,        # latent scores
   Epsilon = Epsilon,  # noise matrix
   C       = C         # empirical covariance
)
saveRDS(sim, "sim_fixed_lambda_k1.rds")
