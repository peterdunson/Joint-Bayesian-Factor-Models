library(rstan)
library(splines)

getwd()

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/functional_comparison/step_one")


fit_func_mgps_direct <- function(df,
                                 stan_file       = "func_mgps.stan",
                                 H               = 3,
                                 M               = 15,
                                 iter            = 2000,
                                 warmup          = 1000,
                                 chains          = 4,
                                 adapt_delta     = 0.9,
                                 maxtree_depth   = 15) {
   
   B <- bs(df$time, df = M, intercept = TRUE)
   
   stan_data <- list(
      N     = length(unique(df$subj)),
      n_obs = nrow(df),
      H     = H,
      M     = M,
      y     = df$y,
      subj  = df$subj,
      B     = as.matrix(B)
   )
   
   sm <- stan_model(stan_file)
   fit <- sampling(
      sm,
      data    = stan_data,
      iter    = iter,
      warmup  = warmup,
      chains  = chains,
      control = list(
         adapt_delta    = adapt_delta,
         max_treedepth  = maxtree_depth
      )
   )
   return(fit)
}


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



# fit the direct MGPS functional factor model
fit1 <- fit_func_mgps_direct(
   df            = smoothed_long,
   stan_file     = "func_mgps.stan",
   H             = 3,
   M             = 15,
   iter          = 2000,
   warmup        = 1000,
   chains        = 4,
   adapt_delta   = 0.99,
   maxtree_depth = 15
)

# check convergence and shrinkage
print(fit1,
      pars  = c("sigma_y", paste0("delta[", 1:3, "]")),
      probs = c(0.1, 0.5, 0.9))


# Save the fit object
saveRDS(fit1, file = "fit1_direct_mgps.rds")

# To load it back in later:
# fit1 <- readRDS("fit1_direct_mgps.rds")



post <- rstan::extract(fit1, permuted=TRUE)
Theta_bar <- apply(post$Theta, c(2,3), mean)  # H×M matrix
t_grid    <- seq(0, 1, length = 200)
B_grid    <- bs(t_grid, df=15, intercept=TRUE)
Lambda_est <- Theta_bar %*% t(B_grid)         # H×200 matrix of curves

matplot(t_grid, t(Lambda_est), type="l", lwd=2,
        xlab="Time (normalized)", ylab="Loading", main="Estimated Factors")

