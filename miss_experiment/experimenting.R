# fit_and_save.R
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- CHOOSE SCENARIO ----
scenario  <- 2   # 1, 2 or 3
sim_path <- sprintf("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen%d_1000.rds", scenario)
sim       <- readRDS(sim_path)

Y  <- scale(sim$Y, center = TRUE, scale = TRUE)  # n × (p+1)
X  <- Y[, -1]
n  <- nrow(Y)
p  <- ncol(Y)
p_x<- ncol(X)
K  <- 5

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/sparse_bayesian_infinite_factor_model")
mod <- stan_model("mgps_factor_model.stan")
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment")

# ---- 1) Misspecified (X-only) ----
stan_data_x <- list(N = n, P = p_x, K = K, Y = X)
fit_x       <- sampling(mod, data = stan_data_x,
                        chains=4, iter=6000, warmup=3000,
                        seed=12,
                        control=list(adapt_delta=0.99, max_treedepth=15))

# extract & save
post_x   <- rstan::extract(fit_x)
Lambda_x_hat <- apply(post_x$Lambda, c(2,3), mean)
saveRDS(list(post=post_x,
             Lambda_hat=Lambda_x_hat),
        file = sprintf("fit_Xonly_scen%d_scale.rds", scenario))

# ---- 2) Joint ([y,X]) ----
stan_data_j   <- list(N = n, P = p, K = K, Y = Y)
fit_j         <- sampling(mod, data = stan_data_j,
                          chains=4, iter=10000, warmup=5000,
                          seed=12,
                          control=list(adapt_delta=0.99, max_treedepth=15))

# extract & save
post_j   <- rstan::extract(fit_j)
Lambda_j_hat <- apply(post_j$Lambda, c(2,3), mean)
saveRDS(list(post=post_j,
             Lambda_hat=Lambda_j_hat),
        file = sprintf("fit_Joint_scen%d_scale.rds", scenario))

