# FULL BASS (Zhao et al., JMLR 2016): Mixture of Dense, Sparse, and Null Factors
# Faithful Gibbs Sampler Implementation (core logic)

library(MASS) # for mvrnorm

# --- Data: Y (n x p), number of factors K ---
# Y: data matrix, centered
# K: number of latent factors

BASS_gibbs <- function(Y, K, n_iter = 3000, burn_in = 1500,
                       a_sigma=1, b_sigma=0.3, a_tau=0.5, b_tau=0.5, verbose=TRUE) {
   n <- nrow(Y)
   p <- ncol(Y)
   # Hyperpriors for mixture weights (Dirichlet prior)
   pi_prior <- c(1, 1, 1)  # prior for (sparse, dense, null)
   # Initialization
   X <- matrix(rnorm(n * K), n, K)
   Lambda <- matrix(0, p, K)
   sigma2 <- rep(1, p)
   tau_col <- rep(1, K)          # for dense
   tau_elem <- matrix(1, p, K)   # for sparse
   z_k <- sample(0:2, K, replace = TRUE)  # 0=sparse, 1=dense, 2=null
   pi_vec <- rep(1/3, 3)         # mixture weights
   
   # Storage
   Lambda_store <- array(NA, dim = c(n_iter - burn_in, p, K))
   z_store <- matrix(NA, n_iter - burn_in, K)
   pi_store <- matrix(NA, n_iter - burn_in, 3)
   iter_keep <- 1
   
   for (iter in 1:n_iter) {
      # --- 1. Update X (factors) ---
      for (i in 1:n) {
         V_inv <- diag(K)
         m <- rep(0, K)
         for (j in 1:p) {
            if (all(Lambda[j, ] == 0)) next
            V_inv <- V_inv + (Lambda[j, ] %*% t(Lambda[j, ])) / sigma2[j]
            m <- m + as.numeric(Y[i, j] * Lambda[j, ] / sigma2[j])
         }
         V <- solve(V_inv)
         mu <- V %*% m
         X[i, ] <- mvrnorm(1, mu, V)
      }
      
      # --- 2. Update z_k for each factor ---
      for (k in 1:K) {
         logp <- rep(-Inf, 3)
         # Sparse (0): elementwise shrinkage
         logp[1] <- sum(dnorm(Lambda[, k], 0, sqrt(tau_elem[, k]), log = TRUE)) + log(pi_vec[1])
         # Dense (1): column shrinkage
         logp[2] <- sum(dnorm(Lambda[, k], 0, sqrt(tau_col[k]), log = TRUE)) + log(pi_vec[2])
         # Null (2): all zero
         logp[3] <- if (all(Lambda[, k] == 0)) 0 else -1e8
         logp[3] <- logp[3] + log(pi_vec[3])
         # Normalize and sample
         maxlp <- max(logp)
         pz <- exp(logp - maxlp)
         pz <- pz / sum(pz)
         z_k[k] <- sample(0:2, 1, prob = pz)
      }
      
      # --- 3. Update Lambda ---
      for (k in 1:K) {
         if (z_k[k] == 2) {
            Lambda[, k] <- 0
         } else {
            for (j in 1:p) {
               res_j <- Y[, j] - X[, -k, drop = FALSE] %*% Lambda[j, -k, drop = FALSE]
               if (z_k[k] == 0) {
                  # Sparse
                  v <- 1 / (sum(X[, k]^2) / sigma2[j] + 1 / tau_elem[j, k])
                  m <- sum(res_j * X[, k]) / sigma2[j] * v
                  Lambda[j, k] <- rnorm(1, m, sqrt(v))
               } else if (z_k[k] == 1) {
                  # Dense
                  v <- 1 / (sum(X[, k]^2) / sigma2[j] + 1 / tau_col[k])
                  m <- sum(res_j * X[, k]) / sigma2[j] * v
                  Lambda[j, k] <- rnorm(1, m, sqrt(v))
               }
            }
         }
      }
      
      # --- 4. Update elementwise and columnwise shrinkage ---
      for (k in 1:K) {
         if (z_k[k] == 0) {
            for (j in 1:p)
               tau_elem[j, k] <- 1 / rgamma(1, a_tau + 0.5, b_tau + 0.5 * Lambda[j, k]^2)
         }
         if (z_k[k] == 1) {
            tau_col[k] <- 1 / rgamma(1, a_tau + 0.5 * p, b_tau + 0.5 * sum(Lambda[, k]^2))
         }
      }
      
      # --- 5. Update mixture weights (Dirichlet) ---
      z_table <- table(factor(z_k, levels = 0:2))
      pi_vec <- as.numeric(rdirichlet(1, pi_prior + z_table))
      
      # --- 6. Update sigma2 (residual variances) ---
      for (j in 1:p) {
         res <- Y[, j] - X %*% Lambda[j, ]
         sigma2[j] <- 1 / rgamma(1, a_sigma + n / 2, b_sigma + 0.5 * sum(res^2))
      }
      
      # --- 7. Store posterior samples ---
      if (iter > burn_in) {
         Lambda_store[iter_keep, , ] <- Lambda
         z_store[iter_keep, ] <- z_k
         pi_store[iter_keep, ] <- pi_vec
         iter_keep <- iter_keep + 1
      }
      if (verbose && iter %% 100 == 0) cat("Iteration", iter, "\n")
   }
   
   return(list(Lambda = Lambda_store, z_k = z_store, pi = pi_store))
}

# ---- Dirichlet sampler (helper function) ----
rdirichlet <- function(n, alpha) {
   x <- matrix(rgamma(length(alpha) * n, shape = alpha), ncol = length(alpha), byrow = TRUE)
   x / rowSums(x)
}
