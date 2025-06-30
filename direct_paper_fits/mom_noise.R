# --- Load pure noise simulation, standardize ---
X <- scale(readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")$Y)

# --- K=1 MoM estimator (top factor) ---
S   <- cov(X)
e   <- eigen(S)
d   <- e$values
V   <- e$vectors

sigma2_hat  <- mean(d[-1])
lambda_mag  <- sqrt(max(d[1] - sigma2_hat, 0))
lambda1_hat <- lambda_mag * V[, 1]
Sigma_hat   <- sigma2_hat * diag(ncol(X))
Sigma_inv   <- solve(Sigma_hat)
denom       <- as.numeric(crossprod(lambda1_hat, Sigma_inv %*% lambda1_hat) + 1)
w_vec       <- (Sigma_inv %*% lambda1_hat) / denom

eta_hat_1   <- as.numeric(X %*% w_vec)   # n-vector

# --- Project out the fitted K=1 factor ---
F1_fit   <- outer(eta_hat_1, lambda1_hat)  # n x p
resid1   <- X - F1_fit

# --- Compare correlation distributions: before vs after ---
get_offdiag <- function(M) M[lower.tri(M)]
cor_before <- cor(X)
cor_after  <- cor(resid1)
r_before <- get_offdiag(cor_before)
r_after  <- get_offdiag(cor_after)

# --- Plot overlayed histograms ---
breaks <- seq(-1, 1, length.out = 41)  # 40 equal-width bins

hist(r_before,
     breaks = breaks,
     col = rgb(0.2, 0.4, 1, 0.5),   # blue, semi-transparent
     main = "Off-diagonal correlations before/after K=1 MoM removal",
     xlab = "Correlation",
     ylab = "Frequency",
     border = "white")

hist(r_after,
     breaks = breaks,
     col = rgb(1, 0.4, 0.4, 0.5),   # red, semi-transparent
     add = TRUE,
     border = "white")

legend("topright",
       legend = c("Before", "After K=1 MoM"),
       fill = c(rgb(0.2, 0.4, 1, 0.5), rgb(1, 0.4, 0.4, 0.5)),
       border = NA)
















# --- Load and standardize ---
X <- scale(readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")$Y)

# --- Function to compute MoM K=1 loadings ---
mom_loading <- function(X) {
   S <- cov(X)
   e <- eigen(S)
   d <- e$values
   V <- e$vectors
   sigma2_hat  <- mean(d[-1])
   lambda_mag  <- sqrt(max(d[1] - sigma2_hat, 0))
   lambda1_hat <- lambda_mag * V[, 1]
   return(lambda1_hat)
}

# --- Bootstrap ---
B <- 1000 # Number of bootstraps
p <- ncol(X)
set.seed(123)
Lambda_boot <- matrix(NA, nrow=B, ncol=p)

for (b in 1:B) {
   idx <- sample(1:nrow(X), nrow(X), replace=TRUE)
   Lambda_boot[b, ] <- mom_loading(X[idx, ])
}

# --- Compute CIs ---
Lambda_hat <- mom_loading(X)
lambda_ci_lower <- apply(Lambda_boot, 2, quantile, 0.025)
lambda_ci_upper <- apply(Lambda_boot, 2, quantile, 0.975)

# --- Summarize ---
summary_table <- data.frame(
   Mean = Lambda_hat,
   CI_2.5 = lambda_ci_lower,
   CI_97.5 = lambda_ci_upper
)
print(summary_table)

# --- Optional: Plot loadings with CIs ---
library(ggplot2)
df_plot <- data.frame(
   Variable = 1:p,
   Mean = Lambda_hat,
   Lower = lambda_ci_lower,
   Upper = lambda_ci_upper
)
ggplot(df_plot, aes(x=Variable, y=Mean)) +
   geom_point() +
   geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.2) +
   labs(title="Bootstrapped 95% CIs for K=1 MoM Loadings", y="Loading") +
   theme_minimal()








#BLUP


# --- Compute BLUP for latent scores using MoM parameters (weights) ---
# Formula: BLUP_eta = w_vec' %*% X[i, ] for each row i

# Recalculate w_vec with the final MoM loadings and noise estimate:
lambda_hat <- Lambda_hat    # Already from your code above
sigma2_hat <- mean(lambda_hat^2) * ((p - 1) / p) # Approx. or use your previous sigma2_hat formula

# This matches your earlier approach:
denom <- as.numeric(crossprod(lambda_hat, lambda_hat) + sigma2_hat)
w_vec_blup <- lambda_hat / denom

# Compute BLUP scores for all observations:
eta_blup <- as.numeric(X %*% w_vec_blup) # n-vector of factor scores

# --- Summarize the BLUPs ---
cat("Summary of BLUP-estimated latent factor scores:\n")
print(summary(eta_blup))
hist(eta_blup, breaks = 40, col = "grey", main = "Histogram of BLUP Factor Scores",
     xlab = "BLUP of Latent Factor", border = "white")

