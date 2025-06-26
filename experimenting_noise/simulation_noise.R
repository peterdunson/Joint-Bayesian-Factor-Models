# --- SIMULATE DATA FOR FACTOR MODEL DIAGNOSTICS ---
set.seed(101)
library(MASS)

# Parameters
n <- 500                 # number of samples
p <- 10                  # number of variables
K_true <- 3              # true number of factors
lambda_strengths <- c(4, 2.5, 1.5)   # All strong factors

# Build Lambda (p x K_true)
Lambda <- matrix(0, nrow=p, ncol=K_true)
for (k in 1:K_true) {
   Lambda[,k] <- rnorm(p, mean=0, sd=lambda_strengths[k])
}

# Simulate factor scores (standard normal)
Eta <- matrix(rnorm(n * K_true), nrow=n, ncol=K_true)

# Signal (n x p)
Y_signal <- Eta %*% t(Lambda)   # [n x p]

# Noise levels
noise_levels <- seq(0.1, 2, by=0.4)
all_data <- list()

for (noise in noise_levels) {
   # Add Gaussian noise (independent, homoscedastic)
   Y_noise <- matrix(rnorm(n * p, mean=0, sd=noise), nrow=n, ncol=p)
   Y <- Y_signal + Y_noise
   Y <- scale(Y, center=TRUE, scale=TRUE)   # center & scale for fair comparison
   all_data[[as.character(noise)]] <- Y
}
