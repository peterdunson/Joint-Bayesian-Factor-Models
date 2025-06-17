# frequentist_coverage_eval.R
# -------------------------------------------
# Evaluate 95% CI coverage for OLS, Lasso, Ridge, PCA+Regression
# using simulation with true beta values
# -------------------------------------------

#IN PROGRESS, NOT RUNNING PROPERLY!!!!


library(glmnet)

# ---- LOAD SIM DATA ----
sim_file <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen1_1000_beta.rds"
sim <- readRDS(sim_file)
X <- sim$Y[, -1]
y <- sim$Y[, 1]
beta_true <- sim$beta

X <- scale(X)
y <- scale(y, center = TRUE, scale = TRUE)

# ---- OLS ----
ols_mod <- lm(y ~ X)
ols_summary <- summary(ols_mod)
coef_ols <- coef(ols_mod)[-1]
se_ols   <- coef(ols_summary)[-1, 2]
ci_ols_lower <- coef_ols - 1.96 * se_ols
ci_ols_upper <- coef_ols + 1.96 * se_ols
coverage_ols <- mean(beta_true >= ci_ols_lower & beta_true <= ci_ols_upper)
cat(sprintf("OLS regression coefficient 95%% coverage: %.2f%%\n", 100 * coverage_ols))

# ---- Ridge ----
set.seed(123)
ridge_fit <- cv.glmnet(X, y, alpha = 0)
beta_ridge <- as.vector(coef(ridge_fit, s = "lambda.min"))[-1]
# Bootstrap CIs
B <- 100
ridge_boot <- matrix(NA, nrow = B, ncol = ncol(X))
for (b in 1:B) {
   idx <- sample(1:nrow(X), replace = TRUE)
   rfit <- cv.glmnet(X[idx, ], y[idx], alpha = 0)
   ridge_boot[b, ] <- as.vector(coef(rfit, s = "lambda.min"))[-1]
}
ridge_se <- apply(ridge_boot, 2, sd)
ci_ridge_lower <- beta_ridge - 1.96 * ridge_se
ci_ridge_upper <- beta_ridge + 1.96 * ridge_se
coverage_ridge <- mean(beta_true >= ci_ridge_lower & beta_true <= ci_ridge_upper)
cat(sprintf("Ridge regression coefficient 95%% coverage: %.2f%%\n", 100 * coverage_ridge))

# ---- Lasso ----
lasso_fit <- cv.glmnet(X, y, alpha = 1)
beta_lasso <- as.vector(coef(lasso_fit, s = "lambda.min"))[-1]
lasso_boot <- matrix(NA, nrow = B, ncol = ncol(X))
for (b in 1:B) {
   idx <- sample(1:nrow(X), replace = TRUE)
   lfit <- cv.glmnet(X[idx, ], y[idx], alpha = 1)
   lasso_boot[b, ] <- as.vector(coef(lfit, s = "lambda.min"))[-1]
}
lasso_se <- apply(lasso_boot, 2, sd)
ci_lasso_lower <- beta_lasso - 1.96 * lasso_se
ci_lasso_upper <- beta_lasso + 1.96 * lasso_se
coverage_lasso <- mean(beta_true >= ci_lasso_lower & beta_true <= ci_lasso_upper)
cat(sprintf("Lasso regression coefficient 95%% coverage: %.2f%%\n", 100 * coverage_lasso))

# ---- PCA + Regression ----
# 1. PCA on X
pca <- prcomp(X, center = FALSE, scale. = FALSE)
expl_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
M <- which(expl_var >= 0.95)[1]      # Number of PCs to explain 95% variance

X_pca <- pca$x[, 1:M, drop = FALSE]  # Scores for regression

# 2. Regression of y on top M PCs
pca_reg <- lm(y ~ X_pca)
pca_coefs <- coef(pca_reg)[-1]         # PC coefficients (not betas on original X)
pca_se    <- summary(pca_reg)$coef[-1, 2]

# To get coefficients on original X, project back:
beta_pca_regression <- pca$rotation[, 1:M] %*% pca_coefs  # [p x 1]

# Approximate SE for each X variable using delta method:
# (SEs are approximate; for precise intervals, use bootstrapping)
boot_pca_betas <- matrix(NA, nrow = B, ncol = ncol(X))
for (b in 1:B) {
   idx <- sample(1:nrow(X), replace = TRUE)
   pca_b <- prcomp(X[idx, ], center = FALSE, scale. = FALSE)
   X_pca_b <- pca_b$x[, 1:M, drop = FALSE]
   reg_b <- lm(y[idx] ~ X_pca_b)
   coefs_b <- coef(reg_b)[-1]
   boot_pca_betas[b, ] <- pca_b$rotation[, 1:M] %*% coefs_b
}
pca_se_betas <- apply(boot_pca_betas, 2, sd)
ci_pca_lower <- as.vector(beta_pca_regression - 1.96 * pca_se_betas)
ci_pca_upper <- as.vector(beta_pca_regression + 1.96 * pca_se_betas)
coverage_pca <- mean(beta_true >= ci_pca_lower & beta_true <= ci_pca_upper)
cat(sprintf("PCA+Regression coefficient 95%% coverage: %.2f%%\n", 100 * coverage_pca))

