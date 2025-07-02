# -----------------------------------------
# Compare k=2: MoM-Eigen, MoM-Trio, MGSP
# -----------------------------------------

# 1. Load & standardize
# dat should be your n×p data matrix
X <- scale(dat, center = TRUE, scale = TRUE)
n <- nrow(X); p <- ncol(X)

# 2. Load MGSP k=2 fit
fit_dir   <- "/Users/.../storing_fit"
fit_MGSP2 <- readRDS(file.path(fit_dir, "MGSP", "fit_joint_NHANES1718_k2.rds"))

# 3. MoM-Eigen (PC) k=2
S      <- cov(X); e <- eigen(S)
d      <- e$values; V <- e$vectors
sigma2_hat    <- mean(d[-(1:2)])
lambda_mag    <- sqrt(pmax(d[1:2] - sigma2_hat, 0))
Lambda_eigen2 <- V[,1:2] %*% diag(lambda_mag)
B_eigen2      <- X %*% Lambda_eigen2
b_hat_eigen2  <- sweep(B_eigen2, 2, colSums(Lambda_eigen2^2), FUN="/")

# 4. MoM-Trio k=2
mom_trio <- function(C) {
   p <- ncol(C); lam2 <- numeric(p)
   for(i in 1:p){
      others <- setdiff(1:p,i); vals <- c()
      for(j in others) for(k in setdiff(others,j)){
         cij <- C[i,j]; cik <- C[i,k]; cjk <- C[j,k]
         if(!is.na(cjk) && abs(cjk)>1e-8) vals <- c(vals,(cij*cik)/cjk)
      }
      lam2[i] <- mean(vals, na.rm=TRUE)
   }
   s <- sign(lam2)*sqrt(abs(lam2))
   s/sqrt(sum(s^2))
}
C0       <- cor(X)
lambda1  <- mom_trio(C0)
C1       <- C0 - lambda1 %*% t(lambda1)
lambda2r <- mom_trio(C1)
# make λ2 orthogonal to λ1
lambda2  <- (lambda2r - sum(lambda1*lambda2r)*lambda1)
lambda2  <- lambda2/sqrt(sum(lambda2^2))
Lambda_trio2 <- cbind(lambda1, lambda2)
B_trio2      <- X %*% Lambda_trio2
b_hat_trio2  <- sweep(B_trio2, 2, colSums(Lambda_trio2^2), FUN="/")

# 5. MGSP k=2
Lambda_mgsp2 <- as.matrix(fit_MGSP2$Lambda_hat[,1:2])
normalize    <- function(v) v/sqrt(sum(v^2))
Lambda_mgsp2 <- apply(Lambda_mgsp2, 2, normalize)
B_mgsp2      <- X %*% Lambda_mgsp2
b_hat_mgsp2  <- sweep(B_mgsp2, 2, colSums(Lambda_mgsp2^2), FUN="/")

# 6. Pairwise scatterplots of all six score series
library(GGally)
df_scores <- data.frame(
   Trio1   = b_hat_trio2[,1],  Trio2   = b_hat_trio2[,2],
   Eigen1  = b_hat_eigen2[,1], Eigen2  = b_hat_eigen2[,2]#,
   #MGSP1   = b_hat_mgsp2[,1],  MGSP2   = b_hat_mgsp2[,2]
)
ggpairs(df_scores, 
        title="k=2 Factor Scores: Trio vs. Eigen vs. MGSP",
        columns=1:4, axisLabels="show")

# 7. Compute off-diagonal residual correlations for each method
get_offdiag <- function(M) M[lower.tri(M)]
resid_cor <- function(B, Lambda) {
   R <- X - B %*% t(Lambda)      # remove 2-factor fit
   get_offdiag(cor(R))
}
r_trio  <- resid_cor(b_hat_trio2,  Lambda_trio2)
r_eigen <- resid_cor(b_hat_eigen2, Lambda_eigen2)
r_mgsp  <- resid_cor(b_hat_mgsp2,  Lambda_mgsp2)

# --- Compute observed off-diagonal correlations ---
r_obs2 <- get_offdiag(cor(X))

# --- 8. Density overlay including observed ---
plot(density(r_obs2),
     lwd   = 2,
     col   = "grey50",
     main  = "Residual Off-Diag Correlations (k=2)",
     xlab  = "Correlation",
     xlim  = c(-1, 1)
)
lines(density(r_trio),  lwd = 2, col = "black")
lines(density(r_eigen), lwd = 2, col = "blue")
lines(density(r_mgsp),  lwd = 2, col = "darkgreen")
legend("topright",
       legend = c("Observed", "MoM-Trio", "MoM-Eigen", "MGSP"),
       col    = c("grey50",  "black",     "blue",       "darkgreen"),
       lwd    = 2
)


