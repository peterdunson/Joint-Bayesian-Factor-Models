library(rstan)
library(splines)


# Wrapper to fit the covariate‐adjusted functional MGPS model
fit_func_mgps_covariate <- function(df, Z,
                                    stan_file   = "func_mgps_covariate.stan",
                                    H           = 3,
                                    M           = 15,
                                    iter        = 2000,
                                    warmup      = 1000,
                                    chains      = 4,
                                    adapt_delta = 0.9) {
   # df: data.frame with columns:
   #     subj (integer 1…N), time (numeric), y (numeric)
   # Z:   data.frame or matrix with N rows (subjects) and J covariate columns
   
   # 2) dimensions
   N <- length(unique(df$subj))
   J <- ncol(Z)
   
   # 3) build spline basis matrix
   B <- bs(df$time, df = M, intercept = TRUE)
   
   # 4) assemble data for Stan
   stan_data <- list(
      N     = N,
      n_obs = nrow(df),
      J     = J,
      H     = H,
      M     = M,
      y     = df$y,
      subj  = df$subj,
      Z     = as.matrix(Z),
      B     = as.matrix(B)
   )
   
   # 5) compile & sample
   sm  <- stan_model(stan_file)
   fit <- sampling(sm,
                   data    = stan_data,
                   iter    = iter,
                   warmup  = warmup,
                   chains  = chains,
                   control = list(adapt_delta = adapt_delta))
   return(fit)
}

# Example usage:
# Suppose df_long has subj, time, y
# and Z_mat is a matrix of covariates (one row per subject):
#
# fit2 <- fit_func_mgps_covariate(df_long, Z_mat,
#                                 stan_file = "func_mgps_covariate.stan",
#                                 H = 4, M = 12,
#                                 iter = 1500, warmup = 750, chains = 2)
#
# print(fit2, pars=c("sigma_y","delta"), probs=c(0.1,0.5,0.9))
#
# To reconstruct loading functions:
# post      <- extract(fit2, permuted=TRUE)
# Theta_bar <- apply(post$Theta, c(2,3), mean)  # H×M
# t_grid    <- seq(min(df_long$time), max(df_long$time), length=200)
# B_grid    <- bs(t_grid, df=M, intercept=TRUE)
# Lambda_est <- Theta_bar %*% t(B_grid)         # H×200
