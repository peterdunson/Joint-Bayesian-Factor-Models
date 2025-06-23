library(ggplot2)
library(reshape2)
library(gridExtra)

# Load data
sim <- readRDS("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
fit_res <- readRDS(sprintf("fit_Xonly_scen2_scale_all.rds"))

# Extract estimated and true loadings (no alignment)
Lambda_x_hat <- fit_res$Lambda_hat
Lambda_true  <- sim$Lambda[-1, ]   # drop y, keep X only

# Prepare data frames for plotting
df_true <- melt(as.matrix(Lambda_true))
colnames(df_true) <- c("Variable", "Factor", "Loading")
df_true$Type <- "True"

df_est <- melt(as.matrix(Lambda_x_hat))
colnames(df_est) <- c("Variable", "Factor", "Loading")
df_est$Type <- "Estimated"

df_true$Variable <- as.factor(df_true$Variable)
df_est$Variable  <- as.factor(df_est$Variable)

# Make plots
p1 <- ggplot(df_true, aes(x=Factor, y=Variable, fill=Loading)) +
   geom_tile() +
   scale_fill_gradient2() +
   ggtitle("True Loadings (X)") +
   theme_minimal()

p2 <- ggplot(df_est, aes(x=Factor, y=Variable, fill=Loading)) +
   geom_tile() +
   scale_fill_gradient2() +
   ggtitle("Estimated Loadings (No Rotation)") +
   theme_minimal()

gridExtra::grid.arrange(p1, p2, ncol=2)



