# ==============================================
# Ciprian‐style K=1 Simulation Analysis (MoM only)
# ==============================================

# 1) Load the fixed‐λ simulation
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim")
sim <- readRDS("sim_fixed_lambda_k1.rds")
X   <- sim$X            # n × P
P   <- ncol(X)

# 2) True λ and λ²
lambdaP  <- sim$lambdaP        # length P
lambdaP2 <- lambdaP^2

# ──────────────────────────────────────────────
# 3) MoM‐Eigen (PC) estimator (k=1)
S            <- cov(X)
e            <- eigen(S)
d            <- e$values
V            <- e$vectors
sigma2hat    <- mean(d[-1])
lambda_eigen <- sqrt(max(d[1] - sigma2hat, 0)) * V[, 1]          # length P
b_eigen      <- as.numeric(X %*% lambda_eigen) / sum(lambda_eigen^2)
pred_eigen   <- lambda_eigen[1] * b_eigen



# ──────────────────────────────────────────────
# 5) MoM‐Trio estimator (k=1)
C        <- cov(X)
lambda2  <- numeric(P)
for (p in 1:P) {
   vals <- c()
   idx  <- setdiff(1:P, p)
   for (i in seq_along(idx)) {
      q <- idx[i]
      for (r in idx[-seq_len(i)]) {
         if (abs(C[q, r]) > 0.02) {
            vals <- c(vals, C[p, q] * C[p, r] / C[q, r])
         }
      }
   }
   lambda2[p] <- mean(vals, na.rm = TRUE)
}
lambda_trio <- sign(lambda2) * sqrt(abs(lambda2))
b_trio      <- as.numeric(X %*% lambda_trio) / sum(lambda_trio^2)
pred_trio   <- lambda_trio[1] * b_trio






#look at this specifcalyy, ALL P Q R THAT ARE DIFFERENT********

vdiff_all <- unlist(lapply(1:P, function(p) {
   idx <- setdiff(1:P,p)
   v <- c()
   for(i in seq_along(idx)) for(r in idx[-seq_len(i)]) {
      if(abs(C[idx[i],r])>0.02)
         v <- c(v, lambdaP2[p]*C[idx[i],r] - C[p,idx[i]]*C[p,r])
   }
   v
}))
hist(vdiff_all, breaks=50,
     main="Distribution of λp²Cqr–CpqCpr",
     xlab="Difference", col=rgb(0.2,0.6,0.2,0.5))
abline(v=0,col="blue",lwd=2)







# 3) Estimate λₚ² by Trio‐MoM: v_pqr = C[p,q]*C[p,r] / C[q,r]
lambda2_est <- numeric(P)
for (p in 1:P) {
   vals <- c()
   others <- setdiff(1:P, p)
   for (i in seq_along(others)) {
      q <- others[i]
      for (r in others[-seq_len(i)]) {
         if (abs(C[q,r]) > 0.02) {
            vals <- c(vals, C[p,q] * C[p,r] / C[q,r])
         }
      }
   }
   lambda2_est[p] <- mean(vals, na.rm = TRUE)
}

# 4) Compute differences: Δ_pqr = λₚ² * C[q,r] - C[p,q]*C[p,r]
#    We’ll just compute summary statistics of these differences for each p
diff_summary <- lapply(1:P, function(p) {
   idx <- setdiff(1:P, p)
   diffs <- c()
   for (i in seq_along(idx)) {
      q <- idx[i]
      for (r in idx[-seq_len(i)]) {
         if (abs(C[q,r]) > 0.02) {
            diffs <- c(diffs, lambdaP2[p] * C[q,r] - C[p,q] * C[p,r])
         }
      }
   }
   # return mean and sd of the differences
   c(mean = mean(diffs, na.rm=TRUE), sd = sd(diffs, na.rm=TRUE))
})

# 5) Display results
# True vs estimated λ²
results <- data.frame(
   p       = 1:P,
   lambda2_true = lambdaP2,
   lambda2_est  = lambda2_est
)
print(results)

# Difference summaries
diff_df <- do.call(rbind, diff_summary)
rownames(diff_df) <- paste0("p", 1:P)
print(diff_df)







# assuming you already have X, lambda_trio, b_trio, lambda_eigen, b_eigen from above…

plot_resid_vs_b <- function(b, lambda, X, method_name) {
   P <- ncol(X)
   par(mfrow = c(3, 4), mar = c(4,4,2,1))  # 12 panels at a time
   for (p in 1:P) {
      eps <- X[, p] - lambda[p] * b
      corr_val <- cor(b, eps)
      plot(b, eps,
           xlab = "b_i", ylab = expression(epsilon[ ip ]),
           main = paste0(method_name, ": p=", p, "\ncor=", round(corr_val,2)),
           pch  = 20, cex = 0.6)
      abline(h=0, col="gray")
   }
   par(mfrow = c(1,1))
}

# MoM‐Trio residual check
plot_resid_vs_b(b_trio,  lambda_trio,  X, "MoM Trio")

# MoM‐Eigen residual check
plot_resid_vs_b(b_eigen, lambda_eigen, X, "MoM Eigen")






