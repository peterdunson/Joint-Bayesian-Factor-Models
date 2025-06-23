library(ggplot2)
library(reshape2)

# ---- CHOOSE SCENARIO ----
scenario <- 2
K        <- 5

# ---- LOAD SIM & FITS ----
sim_path <- sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds", scenario)
sim      <- readRDS(sim_path)

Y        <- scale(sim$Y, center = TRUE, scale = FALSE)
X        <- Y[, -1]
p        <- ncol(Y)
p_x      <- ncol(X)

fit_x_list <- readRDS(sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment/fit_Xonly_scen%d.rds", scenario))
fit_j_list <- readRDS(sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment/fit_Joint_scen%d.rds",  scenario))

post_x       <- fit_x_list$post
Lambda_x_hat <- fit_x_list$Lambda_hat
post_j       <- fit_j_list$post
Lambda_j_hat <- fit_j_list$Lambda_hat

# ---- FIRST FACTOR LINE PLOT ----
df1 <- data.frame(
   Variable = 1:p,
   Joint    = Lambda_j_hat[,1],
   X_only   = c(NA, Lambda_x_hat[,1])
)
df1m <- melt(df1, id.vars="Variable", variable.name="Model", value.name="Loading")

ggplot(df1m, aes(x = Variable, y = Loading, colour = Model)) +
   geom_line() +
   scale_colour_manual(values = c("Joint" = "blue", "X_only" = "red")) +
   labs(title = "First Factor Loadings",
        x     = "Variable",
        y     = "Loading") +
   theme_minimal()

# ---- COVERAGE & MSE HEATMAPS ----

# True loadings
Lambda_true_j <- sim$Lambda           # p × K
Lambda_true_x <- sim$Lambda[-1, ]     # p_x × K

# Posterior draws
draws_j <- post_j$Lambda              # iter × p × K
draws_x <- post_x$Lambda              # iter × p_x × K

# 95% intervals & posterior means
lower_j <- apply(draws_j, c(2,3), quantile, probs = 0.025)
upper_j <- apply(draws_j, c(2,3), quantile, probs = 0.975)
mean_j  <- Lambda_j_hat

lower_x <- apply(draws_x, c(2,3), quantile, probs = 0.025)
upper_x <- apply(draws_x, c(2,3), quantile, probs = 0.975)
mean_x  <- Lambda_x_hat

# Coverage & MSE
cov_j <- (Lambda_true_j >= lower_j & Lambda_true_j <= upper_j)
cov_x <- (Lambda_true_x >= lower_x & Lambda_true_x <= upper_x)
mse_j <- (mean_j - Lambda_true_j)^2
mse_x <- (mean_x - Lambda_true_x)^2

# Pad X-only with NA for the 'y' row
cov_x_full <- rbind(NA, cov_x)
mse_x_full <- rbind(NA, mse_x)

# Melt into long format
df_j_cov <- melt(cov_j, varnames = c("Variable","Factor"), value.name = "Covered")
df_j_mse <- melt(mse_j, varnames = c("Variable","Factor"), value.name = "MSE")
df_j     <- merge(df_j_cov, df_j_mse, by = c("Variable","Factor")); df_j$Model <- "Joint"

df_x_cov <- melt(cov_x_full, varnames = c("Variable","Factor"), value.name = "Covered")
df_x_mse <- melt(mse_x_full, varnames = c("Variable","Factor"), value.name = "MSE")
df_x     <- merge(df_x_cov, df_x_mse, by = c("Variable","Factor")); df_x$Model <- "X-only"

df_all <- rbind(df_j, df_x)
df_all$Covered <- factor(df_all$Covered, levels = c(FALSE,TRUE), labels = c("No","Yes"))

# Coverage heatmap
ggplot(df_all, aes(x = factor(Factor), y = factor(Variable), fill = Covered)) +
   geom_tile(color = "white") +
   facet_wrap(~Model, ncol = 2) +
   scale_fill_manual(values = c("No" = "grey80", "Yes" = "black")) +
   labs(title = "95% Credible Interval Coverage",
        x     = "Factor",
        y     = "Variable") +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# MSE heatmap
ggplot(df_all, aes(x = factor(Factor), y = factor(Variable), fill = MSE)) +
   geom_tile() +
   facet_wrap(~Model, ncol = 2) +
   scale_fill_gradientn(colours = c("white","orange","red")) +
   labs(title = "Per parameter MSE",
        x     = "Factor",
        y     = "Variable") +
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))




# ---- HEATMAP OF FACTOR LOADINGS ----

# Build full X‐only with NA for the 'y' row
Lambda_x_full <- rbind(NA, Lambda_x_hat)  

# Melt into long form
dfL_j <- melt(Lambda_j_hat,    varnames=c("Variable","Factor"), value.name="Loading")
dfL_j$Model <- "Joint"
dfL_x <- melt(Lambda_x_full,   varnames=c("Variable","Factor"), value.name="Loading")
dfL_x$Model <- "X-only"

dfL <- rbind(dfL_j, dfL_x)

# Plot
ggplot(dfL, aes(x=factor(Factor), y=factor(Variable), fill=Loading)) +
   geom_tile() +
   facet_wrap(~Model, ncol=2) +
   scale_fill_gradient2(
      low      = "darkgreen",
      mid      = "white",
      high     = "darkred",
      midpoint = 0,
      name     = "Loading"
   ) +
   labs(title="Heatmap of Factor Loadings",
        x="Factor", y="Variable") +
   theme_minimal() +
   theme(axis.text.x = element_text(angle=90, vjust=0.5),
         panel.grid = element_blank())



#Near 0 and MSE

# ---- Fraction of near-zero loadings (< 0.05) & Train/Test MSE ----
thresh <- 0.05

# pad X-only so dims match (NA for the 'y' row)
Lambda_x_full <- rbind(NA, Lambda_x_hat)

# fraction near zero
frac_joint <- mean(abs(Lambda_j_hat)  < thresh, na.rm=TRUE)
frac_xonly <- mean(abs(Lambda_x_full) < thresh, na.rm=TRUE)

cat("\nFraction of near-zero loadings (< 0.05):\n")
cat(sprintf("  Joint:   %.2f%%\n", 100 * frac_joint))
cat(sprintf("  X-only:  %.2f%%\n", 100 * frac_xonly))

# -- Train/Test split for MSE --
set.seed(123)
n          <- nrow(Y)
train_idx  <- sample(n, floor(0.7 * n))
test_idx   <- setdiff(seq_len(n), train_idx)

# X-only: posthoc regression on scores
eta_x_hat  <- apply(post_x$eta, c(2,3), mean)    # n × K
df_train_x <- data.frame(y = Y[train_idx,1], eta_x_hat[train_idx,])
lm_x       <- lm(y ~ ., data = df_train_x)
pred_train_x <- predict(lm_x, newdata = data.frame(eta_x_hat[train_idx,]))
pred_test_x  <- predict(lm_x, newdata = data.frame(eta_x_hat[test_idx,]))
mse_train_x  <- mean((Y[train_idx,1] - pred_train_x)^2)
mse_test_x   <- mean((Y[test_idx,1]   - pred_test_x )^2)

# Joint model: direct factor prediction
# Joint model: direct factor prediction (fixed)
eta_j_hat   <- apply(post_j$eta, c(2,3), mean)   # n × K
y_hat_j_all <- eta_j_hat %*% matrix(Lambda_j_hat[1, ], ncol=1)  
mse_train_j <- mean((Y[train_idx,1] - y_hat_j_all[train_idx])^2)
mse_test_j  <- mean((Y[test_idx,1]   - y_hat_j_all[test_idx])^2)


cat("\nTrain/Test MSEs:\n")
cat(sprintf("  X-only  - Train: %.4f | Test: %.4f\n", mse_train_x, mse_test_x))
cat(sprintf("  Joint   - Train: %.4f | Test: %.4f\n", mse_train_j, mse_test_j))



