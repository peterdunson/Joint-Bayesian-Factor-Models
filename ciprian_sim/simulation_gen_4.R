# simulation_gen.R
# Fixed λₚ for p = 1…10, zeros for p = 11…20; n = 1000, heteroskedastic noise

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim")

set.seed(26031980)

# 1) Parameters
P       <- 20
n       <- 1000
sigmasq <- 0.2

# 2) Fixed loadings vector λₚ
lambdaP <- c(rep(1 / sqrt(10), 10), rep(0, 10))  # length P

# 3) Simulate heteroskedastic noise εᵢₚ ∼ N(0, σₚ²)
#    draw P variances around sigmasq
sigma2_p <- runif(P, min = 0.5 * sigmasq, max = 1.5 * sigmasq)
Epsilon  <- matrix(rnorm(n * P), nrow = n, ncol = P) * 
   matrix(sqrt(sigma2_p), n, P, byrow = TRUE)

# 4) Simulate latent scores bᵢ ∼ N(0,1)
b <- rnorm(n)

# 5) Build data Xᵢₚ = λₚ bᵢ + εᵢₚ
X <- outer(b, lambdaP) + Epsilon

# 6) Empirical covariance and quick check
C <- cov(X)
hist(C, breaks = 20, main = "Histogram of Covariances", xlab = "C[p,q]")

# 7) Pack and save
sim <- list(
   X       = X,        # n × P data matrix
   lambdaP = lambdaP,  # true loadings
   b       = b,        # latent scores
   Epsilon = Epsilon,  # noise matrix
   C       = C         # empirical covariance
)
saveRDS(sim, "sim_fixed_lambda_k1_hetnoise.rds")

