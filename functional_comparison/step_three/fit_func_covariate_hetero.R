library(rstan)
library(splines)


# Wrapper to fit the covariate‐adjusted + heteroskedastic functional MGPS model
fit_func_mgps_cov_het <- function(df, Z,
                                  stan_file   = "func_mgps_cov_het.stan",
                                  H           = 3,
                                  M           = 15,
                                  iter        = 2000,
                                  warmup      = 1000,
                                  chains      = 4,
                                  adapt_delta = 0.9) {
   # df: data.frame with columns subj (1…N), time (numeric), y (numeric)
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
# df_long <- data.frame(subj = rep(1:50, each=100),
#                       time = runif(5000),
#                       y    = rnorm(5000))
# Z_mat    <- matrix(rnorm(50*3), ncol=3)  # e.g. 3 covariates for 50 subjects
#
# fit3 <- fit_func_mgps_cov_het(df_long, Z_mat,
#                               stan_file = "func_mgps_cov_het.stan",
#                               H = 4, M = 12,
#                               iter = 1500, warmup = 750, chains = 2)
#
# print(fit3, pars=c("sigma_beta","sigma_gamma","delta"), probs=c(0.1,0.5,0.9))
