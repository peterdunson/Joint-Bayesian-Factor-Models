# full_mgps_example.R

# ─── 1) Libraries ─────────────────────────────────────────────────────────
library(mvtnorm)       # for rmvnorm
library(magrittr)      # for %>%
library(FactoMineR)    # for RV coefficient
library(tidyr)         # for data manipulation if needed

# ─── 2) Helper functions ──────────────────────────────────────────────────

# Multivariate normal sampler
rmvn <- function(n, mu, Sigma) {
   mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma)
}

# RV coefficient wrapper returning list with $rv
coeffRV <- function(A, B) {
   # A and B must be same dimension matrices
   list(rv = FactoMineR::RVcoef(A, B))
}

# Get mode of a vector
get_mode <- function(x) {
   ux <- unique(x)
   ux[which.max(tabulate(match(x, ux)))]
}

# Metropolis–Hastings step
MH <- function(init = 0, M = 200, num_burn = 100, step_size = 0.1, func, ...) {
   param <- numeric(M)
   param[1] <- init
   for (i in 2:M) {
      cand <- rnorm(1, mean = param[i-1], sd = step_size)
      log_r <- func(cand, ...) - func(param[i-1], ...)
      if (log(runif(1)) < log_r) {
         param[i] <- cand
      } else {
         param[i] <- param[i-1]
      }
   }
   param[M]
}

# ══════════════════════════════════════════════════════════════════════════
# Multiplicative Gamma Process Shrinkage (MGPS) Gibbs sampler
# ══════════════════════════════════════════════════════════════════════════
MGPS <- function(Y,
                 num_iter   = 5000,
                 num_burn   = 2500,
                 thin       = 5,
                 eta0,
                 Lambda0,
                 Sigma0) {
   n <- nrow(Y)
   p <- ncol(Y)
   # initial # factors
   k_ast <- ncol(Lambda0)
   
   # hyperparameters
   a_sigma <- 1;    b_sigma <- 0.3
   a1      <- 2.1;  b1      <- 1
   a2      <- 3.1;  b2      <- 1
   df      <- 3
   alpha0  <- -1;   alpha1  <- -5e-4
   epsilon <- 1e-4
   
   # initialize
   Lambda    <- matrix(0, p, k_ast)
   sigma_inv <- rgamma(p, shape = a_sigma, rate = b_sigma)
   phi       <- matrix(rgamma(p*k_ast, df/2, df/2), nrow = p, ncol = k_ast)
   delta     <- c(rgamma(1, shape = a1, rate = b1),
                  rgamma(k_ast-1, shape = a2, rate = b2))
   tau       <- cumprod(delta)
   
   num_slice    <- (num_iter - num_burn) / thin
   cube_Lambda  <- array(0, dim = c(p, k_ast, num_slice))
   cube_eta     <- array(0, dim = c(n, k_ast, num_slice))
   mult         <- matrix(0, n, p)
   cov_epsilon  <- matrix(0, p, p)
   k_est        <- numeric(num_slice)
   
   slice_idx <- 0
   for (iter in 1:num_iter) {
      # 1) Update eta
      V_eta   <- chol2inv(chol(diag(k_ast) + t(Lambda) %*% diag(sigma_inv) %*% Lambda))
      m_eta   <- V_eta %*% t(Lambda) %*% diag(sigma_inv) %*% t(Y)
      eta     <- t(m_eta) + rmvn(n, rep(0, k_ast), V_eta)
      
      # 2) Update Lambda
      Lambda <- t(sapply(1:p, function(j) {
         D_inv     <- diag(phi[j, ] * tau, k_ast)
         V_lambda  <- chol2inv(chol(D_inv + sigma_inv[j] * t(eta) %*% eta))
         m_lambda  <- as.vector(sigma_inv[j] * V_lambda %*% t(eta) %*% Y[, j])
         rmvn(1, m_lambda, V_lambda)
      }))
      
      # 3) Update phi
      rate_phi <- (df/2) + (Lambda^2) %*% diag(tau)
      phi      <- t(apply(rate_phi, 1, function(rp) {
         rgamma(k_ast, shape = (df+1)/2, rate = rp)
      }))
      
      # 4) Update delta & tau
      phiLambda <- colSums(phi * Lambda^2)
      for (h in 1:k_ast) {
         ad <- if (h==1) a1 + p*k_ast/2 else a2 + p*(k_ast-h+1)/2
         bd <- 1 + 0.5 / delta[h] * sum(tau[h:k_ast] * phiLambda[h:k_ast])
         delta[h] <- rgamma(1, shape = ad, rate = bd)
         tau       <- cumprod(delta)
      }
      
      # 5) Update sigma_inv
      res        <- (Y - eta %*% t(Lambda))^2
      sigma_inv  <- sapply(1:p, function(j) {
         rgamma(1,
                shape = a_sigma + n/2,
                rate  = b_sigma + 0.5 * sum(res[, j]))
      })
      cov_epsilon <- cov_epsilon + diag(1/sigma_inv) / num_slice
      
      # 6) Adapt # factors
      prob      <- exp(alpha0 + alpha1 * iter)
      if (runif(1) < prob) {
         prop_zero <- colMeans(abs(Lambda) < epsilon)
         idx_red   <- which(prop_zero >= 1)
         if (iter>20 && length(idx_red)==0 && all(prop_zero<.995)) {
            # birth
            k_ast <- k_ast + 1
            Lambda    <- cbind(Lambda, rep(0, p))
            eta       <- cbind(eta, rnorm(n))
            phi       <- cbind(phi, rgamma(p, df/2, df/2))
            delta     <- c(delta, rgamma(1, a2, b2))
            tau       <- cumprod(delta)
         } else if (length(idx_red)>0) {
            # death
            k_ast <- max(k_ast - length(idx_red), 1)
            Lambda    <- Lambda[, -idx_red, drop=FALSE]
            eta       <- eta   [, -idx_red, drop=FALSE]
            phi       <- phi   [, -idx_red, drop=FALSE]
            delta     <- delta[-idx_red]
            tau       <- cumprod(delta)
         }
      }
      
      # 7) Store posterior samples
      if (iter > num_burn && (iter - num_burn) %% thin == 0) {
         slice_idx <- slice_idx + 1
         cube_Lambda[ , , slice_idx] <- Lambda
         cube_eta   [ , , slice_idx] <- eta
         mult       <- mult + eta %*% t(Lambda) / num_slice
         k_est[slice_idx] <- k_ast
      }
      
      if (iter %% 1000 == 0) cat("iter =", iter, "\n")
   }
   
   # posterior medians + RV metrics
   Lambda_hat   <- apply(cube_Lambda, c(1,2), median)
   eta_hat      <- apply(cube_eta,    c(1,2), median)
   RV_Lambda    <- coeffRV(Lambda_hat, Lambda0)$rv
   RV_eta       <- coeffRV(eta_hat,    eta0)$rv
   
   list(
      Lambda_hat = Lambda_hat,
      eta_hat    = eta_hat,
      RV_Lambda  = RV_Lambda,
      RV_eta     = RV_eta,
      k_est      = get_mode(k_est),
      mult       = mult,
      cov_epsilon = cov_epsilon
   )
}
# ══════════════════════════════════════════════════════════════════

# ─── 3) Run MGPS on Scenario 2 simulation ────────────────────────────
sim    <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Y      <- sim$Y
eta0   <- sim$eta
Lambda0<- sim$Lambda
Sigma0 <- sim$Omega

res <- MGPS(Y, num_iter=5000, num_burn=2500, thin=5,
            eta0=eta0, Lambda0=Lambda0, Sigma0=Sigma0)

# ─── 4) Print results ────────────────────────────────────────────────
cat("Estimated number of factors:", res$k_est, "\n")
cat("RV( Lambda ): ", res$RV_Lambda, "\n")
cat("RV( eta    ): ", res$RV_eta,    "\n")
