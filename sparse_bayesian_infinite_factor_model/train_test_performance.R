# train_test_performance.R
# For joint Bayesian factor models (predicting outcome Y from X via factors)

library(rstan)

# --- USER: Load and preprocess your data ---
# Suppose your data is an n x (p+1) matrix called "XY"
# - Columns 1:p: predictors X
# - Column p+1: outcome Y
# Example: XY <- as.matrix(read.csv("your_data.csv"))

# --- Standardize data ---
XY <- scale(XY, center = TRUE, scale = TRUE)   # Center and scale each column

n <- nrow(XY)
p <- ncol(XY) - 1       # Number of predictors
K <- 20                 # Upper bound for number of factors

# --- Train/test split ---
set.seed(123)
train_idx <- sample(1:n, size = floor(0.75 * n), replace = FALSE)
test_idx <- setdiff(1:n, train_idx)

XY_train <- XY[train_idx, ]
XY_test <- XY[test_idx, ]
n_train <- nrow(XY_train)
n_test <- nrow(XY_test)

# --- Stan data list (joint model: all variables, predictors + outcome) ---
stan_data <- list(
   N = n_train,
   P = p + 1,
   K = K,
   Y = XY_train
)

# --- Fit model to training data ---
fit <- stan(
   file = "mgps_factor_model.stan",
   data = stan_data,
   iter = 2000,
   warmup = 1000,
   chains = 4,
   seed = 42,
   control = list(adapt_delta = 0.95)
)

# --- Posterior means for factor loadings and residuals ---
post <- extract(fit)
Lambda_mean <- apply(post$Lambda, c(2,3), mean)  # (p+1) x K
psi_mean <- apply(post$psi, 2, mean)             # (p+1)

# --- Predict test Y using only test X and trained loadings/psi ---
predict_Y_joint <- function(X_test, Lambda, psi, y_col) {
   # X_test: n_test x p  (predictors only)
   # Lambda: (p+1) x K
   # psi: length p+1
   # y_col: which column in Lambda/psi is the outcome (usually p+1)
   n_test <- nrow(X_test)
   p <- ncol(X_test)
   K <- ncol(Lambda)
   preds <- rep(NA, n_test)
   
   # Partition loadings and psi for X and Y
   Lambda_x <- Lambda[1:p, , drop=FALSE]
   Lambda_y <- Lambda[y_col, , drop=FALSE]
   psi_x <- psi[1:p]
   psi_y <- psi[y_col]
   
   Psi_x_inv <- diag(1 / psi_x)
   Lt_PsiInv <- t(Lambda_x) %*% Psi_x_inv
   V_eta <- solve(diag(K) + Lt_PsiInv %*% Lambda_x)
   
   for (i in 1:n_test) {
      M_eta <- V_eta %*% Lt_PsiInv %*% X_test[i, ]
      preds[i] <- as.numeric(Lambda_y %*% M_eta)   # predicted mean of Y
   }
   preds
}

X_test <- XY_test[, 1:p, drop=FALSE]
Y_test_true <- XY_test[, p+1]
Y_pred <- predict_Y_joint(X_test, Lambda_mean, psi_mean, y_col = p+1)

# --- Compute test set MSE for outcome Y ---
mse <- mean((Y_test_true - Y_pred)^2)
cat("Test set MSE for Y:", mse, "\n")

