# Reconstruct fitted curves from posterior mean
post   <- rstan::extract(fit1)
Theta  <- apply(post$Theta, c(2,3), mean)   # H x M
eta    <- apply(post$eta, c(2,3), mean)     # N x H

# Compute loadings on grid
t_grid <- seq(0, 1, length=50)
B_grid <- bs(t_grid, df=15, intercept=TRUE)
Lambda <- Theta %*% t(B_grid)    # H x 50

# Reconstruct each subject's curve
y_hat <- eta %*% Lambda          # N x 50

# Compare to smoothed curves
matplot(t_grid, t(y_hat[1:10,]), type="l")    # plot for 10 subjects





# Compare estimated Lambda (columns) to PCA of observed curves
# (Assume smoothed_matrix: rows=time, cols=subject)
smoothed_matrix <- matrix(smoothed_long$y, nrow=50, byrow=FALSE)
pca_fit <- prcomp(t(smoothed_matrix), center=TRUE, scale.=FALSE)
matplot(t_grid, pca_fit$rotation[,1:3], type="l", main="PCA loadings")
matplot(t_grid, t(Lambda[1:3,]), type="l", lty=2, add=TRUE)

