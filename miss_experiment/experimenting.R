library(rstan)
library(ggplot2)

# ---- CHOOSE SCENARIO ----
scenario <- 2   # Change this to 1, 2, or 3 as needed
sim_path <- sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds", scenario)
sim <- readRDS(sim_path)

Y <- sim$Y  # (n x p+1) matrix, first column is y
n <- nrow(Y)
p <- ncol(Y)
K <- 5

Y <- scale(Y, center = TRUE, scale = FALSE)
X <- Y[, -1]
p_x <- ncol(X)

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
mod <- stan_model("mgps_factor_model.stan")

# ---- Misspecified Model: Fit on X only ----
stan_data_x <- list(N = n, P = p_x, K = K, Y = X)
fit_fa_x <- sampling(
   object = mod,
   data = stan_data_x,
   chains = 4,
   iter = 4000,
   warmup = 2000,
   seed = 12,
   control = list(adapt_delta = 0.99, max_treedepth = 12)
)
post_x <- rstan::extract(fit_fa_x)
eta_hat_x <- apply(post_x$eta, c(2,3), mean)    # n x K

# Posthoc regression of y on estimated factors
lm_posthoc <- lm(Y[,1] ~ eta_hat_x)
y_pred_posthoc <- predict(lm_posthoc)
mse_posthoc <- mean((Y[,1] - y_pred_posthoc)^2)

# ---- Joint Model: Fit on [y, X] ----
stan_data_joint <- list(N = n, P = p, K = K, Y = Y)
fit_fa_joint <- sampling(
   object = mod,
   data = stan_data_joint,
   chains = 4,
   iter = 4000,
   warmup = 2000,
   seed = 12,
   control = list(adapt_delta = 0.99, max_treedepth = 12)
)
post_joint <- rstan::extract(fit_fa_joint)
eta_hat_joint <- apply(post_joint$eta, c(2,3), mean)            # n x K
Lambda_hat_joint <- apply(post_joint$Lambda, c(2,3), mean)      # p x K

# Predict y directly from joint factor model
y_pred_joint <- eta_hat_joint %*% t(Lambda_hat_joint)[,1]
mse_joint <- mean((Y[,1] - y_pred_joint)^2)

# ---- Visualization: Compare first factor loadings ----
Lambda_x_hat <- apply(post_x$Lambda, c(2,3), mean)   # (p_x x K)
Lambda_joint_hat <- Lambda_hat_joint                 # (p x K)

df_loadings <- data.frame(
   Variable = 1:p,
   Joint_Loading = Lambda_joint_hat[,1],
   Misspecified_Loading = c(NA, Lambda_x_hat[,1]) # NA for outcome var not in X-only
)

ggplot(df_loadings, aes(x = Variable)) +
   geom_line(aes(y = Joint_Loading, color = "Joint (with y)")) +
   geom_line(aes(y = Misspecified_Loading, color = "Misspecified (X only)")) +
   labs(title = "First Factor Loading: Joint vs Misspecified", y = "Loading Value") +
   scale_color_manual(values = c("blue", "red")) +
   theme_minimal()

# ---- Evaluation Output ----
cat(sprintf("MSE (Misspecified/posthoc): %.4f\n", mse_posthoc))
cat(sprintf("MSE (Joint factor model): %.4f\n", mse_joint))

# Optional: Compare factor loading to ground truth if available
if (!is.null(sim$Lambda)) {
   Lambda_true <- sim$Lambda
   cor1 <- cor(Lambda_joint_hat[,1], Lambda_true[,1])
   cat(sprintf("Correlation (first loading: Joint vs Truth): %.3f\n", cor1))
}
