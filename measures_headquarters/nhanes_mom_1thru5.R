# --- Prepare data ---
X_raw <- as.matrix(phthalates_adults[, phthalate_vars])
X <- scale(X_raw)
get_offdiag <- function(M) M[lower.tri(M)]

maxK <- 5
offdiag_list <- list()
offdiag_list[[1]] <- get_offdiag(cor(X))  # K=0 (original)

resid_mat <- X
for (K in 1:maxK) {
   # Covariance and eigendecomp of current residual
   S   <- cov(resid_mat)
   e   <- eigen(S)
   d   <- e$values
   V   <- e$vectors
   
   sigma2_hat  <- mean(d[-(1:K)])
   lambda_mag  <- sqrt(max(d[1] - sigma2_hat, 0))
   lambda_hat  <- lambda_mag * V[, 1]      # top eigenvector, scaled
   
   Sigma_hat   <- sigma2_hat * diag(ncol(X))
   Sigma_inv   <- solve(Sigma_hat)
   denom       <- as.numeric(crossprod(lambda_hat, Sigma_inv %*% lambda_hat) + 1)
   w_vec       <- (Sigma_inv %*% lambda_hat) / denom
   
   eta_hat     <- as.numeric(resid_mat %*% w_vec)  # n-vector
   F_fit       <- outer(eta_hat, lambda_hat)       # n x p
   resid_mat   <- resid_mat - F_fit                # update residual
   
   # Store off-diagonal correlations
   offdiag_list[[K+1]] <- get_offdiag(cor(resid_mat))
}

# --- Plot overlayed histograms ---
breaks <- seq(-1, 1, length.out = 41)
cols <- c(rgb(0.2, 0.4, 1, 0.5), rgb(1, 0.7, 0.2, 0.4), rgb(0.2, 1, 0.7, 0.4), 
          rgb(1, 0.3, 0.7, 0.4), rgb(1, 0.4, 0.4, 0.4), rgb(0.4, 0.4, 0.4, 0.4))

plot(NULL, xlim = c(-1, 1), ylim = c(0, max(sapply(offdiag_list, function(x) hist(x, breaks=breaks, plot=FALSE)$counts))),
     xlab = "Correlation", ylab = "Frequency",
     main = "Off-diagonal correlations after sequential MoM factor removal")

for (K in 0:maxK) {
   hist(offdiag_list[[K+1]], breaks = breaks, col = cols[K+1], add = K > 0, border = "white")
}
legend("topright",
       legend = c("Before", paste0("After K=", 1:maxK, " MoM")),
       fill = cols[1:(maxK+1)],
       border = NA)

# --- Print summary statistics ---
for (K in 0:maxK) {
   cat(sprintf("Mean abs off-diagonal correlation AFTER K=%d: %.4f\n", K, mean(abs(offdiag_list[[K+1]]))))
}



for (K in 0:maxK) {
   hist(offdiag_list[[K+1]],
        breaks = breaks,
        col = cols[K+1],
        border = "white",
        xlim = c(-1, 1),
        ylim = c(0, max(sapply(offdiag_list, function(x) hist(x, breaks=breaks, plot=FALSE)$counts))),
        main = sprintf("Off-diagonal correlations after K=%d MoM factors removed", K),
        xlab = "Correlation", ylab = "Frequency")
   legend("topright",
          legend = if (K == 0) "Original" else paste("After K=", K),
          fill = cols[K+1],
          border = NA)
   Sys.sleep(1)
}



# --- Assumes offdiag_list (K=0:5), breaks, cols are defined as above ---
mean_before <- mean(abs(offdiag_list[[1]]))
means_after <- sapply(2:6, function(k) mean(abs(offdiag_list[[k]])))

par(mfrow = c(2, 3)) # Arrange 2 rows x 3 columns

for (K in 1:5) {
   hist(offdiag_list[[1]],   # before
        breaks = breaks,
        col = rgb(0.2, 0.4, 1, 0.5),
        main = sprintf("Before vs After K=%d MoM", K),
        xlab = "Correlation",
        ylab = "Frequency",
        border = "white",
        xlim = c(-1, 1),
        ylim = c(0, max(sapply(offdiag_list, function(x) hist(x, breaks=breaks, plot=FALSE)$counts))))
   hist(offdiag_list[[K+1]], # after K
        breaks = breaks,
        col = rgb(1, 0.4, 0.4, 0.5),
        add = TRUE,
        border = "white")
   legend("topright",
          legend = c(
             sprintf("Before (%.3f)", mean_before),
             sprintf("After K=%d (%.3f)", K, means_after[K])
          ),
          fill = c(rgb(0.2, 0.4, 1, 0.5), rgb(1, 0.4, 0.4, 0.5)),
          border = NA,
          cex = 0.85)
   # Optionally, add text annotation in the bottom right
   # text(0.7, 0.9*par("usr")[4],
   #      labels = sprintf("Mean |before|: %.3f\nMean |after|: %.3f", mean_before, means_after[K]),
   #      adj = 0)
}

par(mfrow = c(1, 1))


