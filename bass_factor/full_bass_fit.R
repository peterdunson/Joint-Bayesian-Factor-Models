# Load simulation data as before
scenario <- 1
sim_path <- sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds", scenario)
sim <- readRDS(sim_path)
Y <- sim$Y   # n x p data matrix

# Optionally center the data columns (recommended for factor analysis)
Y <- scale(Y, center = TRUE, scale = FALSE)

# Set number of factors (choose based on your analysis, e.g., K=5)
K <- 5

# Run the full BASS_gibbs and rdirichlet functions from above into your script before running this next part!



# Set sampler parameters
n_iter <- 3000
burn_in <- 1500

# Run BASS
set.seed(1)
out <- BASS_gibbs(Y, K, n_iter = n_iter, burn_in = burn_in)


# Posterior mean of factor loadings
Lambda_post <- apply(out$Lambda, c(2,3), mean)  # (p x K) matrix

# For factor types across posterior (0=sparse, 1=dense, 2=null)
z_k_post <- apply(out$z_k, 2, function(x) as.integer(names(sort(table(x), decreasing=TRUE))[1]))

cat("Most likely factor types (by column):\n")
print(z_k_post)

# Plot the estimated loadings matrix
image(Lambda_post, main = "Posterior Mean of Loadings (BASS)")

# If you want to see which factors are null, dense, sparse:
barplot(table(z_k_post), names.arg = c("Sparse", "Dense", "Null"), main = "Factor Types")


# Lambda_post is p x K (variables x factors)
image(t(Lambda_post), main = "Posterior Mean Factor Loadings", 
      xlab = "Variables", ylab = "Factors", axes = FALSE, col = heat.colors(30))
axis(1, at = seq(0, 1, length.out = ncol(Lambda_post)), labels = 1:ncol(Lambda_post))
axis(2, at = seq(0, 1, length.out = nrow(Lambda_post)), labels = 1:nrow(Lambda_post))
box()



