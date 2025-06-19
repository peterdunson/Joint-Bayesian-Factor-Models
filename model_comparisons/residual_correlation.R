# ---------------------------------------------------------------
# Residual Correlation Diagnostics: Avg & Max |Corr| vs. # Factors
# ---------------------------------------------------------------

library(rstan)

# --- File paths (adjust as needed) ---
sim_file <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds"
#fit_file <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/mgps_fit_scen2_5.rds"
fit_file <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/stan_tebfar_fit_scen2_5.rds"


sim  <- readRDS(sim_file)
fit  <- readRDS(fit_file)
post <- rstan::extract(fit)


Y <- sim$Y                   # n x (p+1)
X <- Y[, -1, drop = FALSE]   # n x p
n <- nrow(X); p <- ncol(X)
K <- dim(post$Lambda)[3]

# --- Posterior-mean loadings & residual precisions ---
Lambda_full <- apply(post$Lambda, c(2,3), mean)  # (p+1) x K
psi_full    <- apply(post$psi,    2,   mean)     # length p+1

# drop the outcome row
Lambda_mean <- Lambda_full[-1, , drop = FALSE]  # p x K
psi_mean    <- psi_full[-1]                     # length p

# --- Preallocate storage ---
max_k       <- K
avg_abs_cor <- numeric(max_k + 1)
max_abs_cor <- numeric(max_k + 1)

# --- Loop over k = 0..K ---
for (k in 0:max_k) {
   
   if (k == 0) {
      # Residuals = centered X
      Resid <- scale(X, center = TRUE, scale = FALSE)
      
   } else {
      # select first k factors
      Lambda_k <- Lambda_mean[, 1:k, drop = FALSE]  # p x k
      
      # build Psi^{-1}
      Psi_inv <- diag(1 / psi_mean)                 # p x p
      
      # compute A = Psi^{-1} * Lambda_k
      A <- Psi_inv %*% Lambda_k                     # p x k
      
      # compute M = I_k + Lambda_k^T * Psi^{-1} * Lambda_k
      M <- diag(k) + t(Lambda_k) %*% A               # k x k
      
      # factor scores: (X %*% A) %*% inv(M)
      Eta_hat <- (X %*% A) %*% solve(M)              # n x k
      
      # fitted values and residuals
      Fitted <- Eta_hat %*% t(Lambda_k)             # n x p
      Resid  <- X - Fitted                           # n x p
   }
   
   # correlation of residuals
   R <- cor(Resid)
   
   # get off-diagonal entries
   off <- abs(R[upper.tri(R)])
   
   # store avg and max
   avg_abs_cor[k + 1] <- mean(off)
   max_abs_cor[k + 1] <- max(off)
}
K_max <- length(avg_abs_cor) - 1

# 1) average absolute residual correlation
plot(0:K_max, avg_abs_cor, type = "b", pch = 19,
     xlab = "Number of Factors (k)",
     ylab = "Avg |Residual Correlation|",
     main = "Average Residual Correlation vs. k")

# 2) maximum absolute residual correlation
plot(0:K_max, max_abs_cor, type = "b", pch = 17, col = "blue",
     xlab = "Number of Factors (k)",
     ylab = "Max |Residual Correlation|",
     main = "TEB-FAR Sim 2")

