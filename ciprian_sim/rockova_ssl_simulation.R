# simulation_gen.R
# Spike-and-Slab LASSO simulation: Block-diagonal correlation, sparse beta0

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim")

set.seed(26031980)

# 1) Parameters
n <- 100
p <- 1000
block_size <- 20
n_blocks <- p / block_size
rho <- 0.9  # within-block correlation

# 2) Build block-diagonal covariance matrix Sigma
block <- matrix(rho, block_size, block_size)
diag(block) <- 1
Sigma <- Matrix::bdiag(replicate(n_blocks, block, simplify = FALSE))
Sigma <- as.matrix(Sigma)

# 3) Generate X ~ MVN(0, Sigma)
library(MASS)
X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

# 4) True beta0 vector: nonzero at 1,51,101,151,201,251
beta0 <- numeric(p)
nz_idx <- c(1, 51, 101, 151, 201, 251)
nz_vals <- c(-2.5, -2, -1.5, 1.5, 2, 2.5) / sqrt(3)
beta0[nz_idx] <- nz_vals

# 5) Generate response y = X %*% beta0 + e
e <- rnorm(n)
y <- as.numeric(X %*% beta0 + e)

# 6) Save the data and parameters
sim <- list(
   X = X,
   y = y,
   beta0 = beta0,
   Sigma = Sigma,
   nz_idx = nz_idx,
   nz_vals = nz_vals
)
saveRDS(sim, "sim_spikeslablasso_block.rds")
