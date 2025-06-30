# --- Load simulation and fit ---
sim <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000_1.rds")
Y   <- scale(sim$Y, center = TRUE, scale = TRUE)

fit_res <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/direct_paper_fits/storing_fit/fit_joint_scen2_k1.rds")
Lambda_hat <- fit_res$Lambda_hat  # [p x K], here K=1

if (nrow(Lambda_hat) != ncol(Y)) Lambda_hat <- t(Lambda_hat) # [p x 1]

# --- Project out the fitted K=1 factor ---
# (Same math as before)
eta_hat1 <- Y %*% Lambda_hat %*% solve(t(Lambda_hat) %*% Lambda_hat)   # n x 1
Y_hat1   <- eta_hat1 %*% t(Lambda_hat)                                 # n x p
resid1   <- Y - Y_hat1

# --- Compare correlation distributions: before vs after ---
get_offdiag <- function(M) M[lower.tri(M)]
cor_before <- cor(Y)
cor_after  <- cor(resid1)
r_before <- get_offdiag(cor_before)
r_after  <- get_offdiag(cor_after)

# --- Plot overlayed histograms ---
breaks <- seq(-1, 1, length.out = 41)

hist(r_before,
     breaks = breaks,
     col = rgb(0.2, 0.4, 1, 0.5),
     main = "Off-diagonal correlations before/after K=1 MGSP removal",
     xlab = "Correlation",
     ylab = "Frequency",
     border = "white")

hist(r_after,
     breaks = breaks,
     col = rgb(1, 0.4, 0.4, 0.5),
     add = TRUE,
     border = "white")

legend("topright",
       legend = c("Before", "After K=1"),
       fill = c(rgb(0.2, 0.4, 1, 0.5), rgb(1, 0.4, 0.4, 0.5)),
       border = NA)
