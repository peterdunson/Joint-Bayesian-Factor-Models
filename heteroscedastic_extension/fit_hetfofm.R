# file: fit_hetfofm.R
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# 1) Simulate toy data
set.seed(123)
N <- 50; p <- 40; H <- 2
t <- seq(0,1,length.out=p)
# true loadings: sin & cos
Lambda_true <- rbind(sin(2*pi*t), cos(2*pi*t))
eta_true    <- matrix(rnorm(N*H), N, H)
signal      <- eta_true %*% Lambda_true
# hetero noise: bump in middle
tau2_true   <- 0.1 + 0.5 * ( (t>0.4)&(t<0.6) )
noise       <- matrix(rnorm(N*p), N, p) * sqrt(rep(tau2_true, each=N))
X_sim       <- signal + noise

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/heteroscedastic_extension")

# 2) Prepare data for Stan
stan_data <- list(
   N = N,
   p = p,
   H = H,
   X = X_sim
)

# 3) Compile & fit
sm <- stan_model("hetfofm.stan")
fit <- sampling(
   sm,
   data = stan_data,
   iter = 2000,
   chains = 4,
   control = list(adapt_delta=0.95, max_treedepth=20)
)

# suppose you set p <- 40 in your R script
last_idx <- p

print(fit,
      pars = c(
         "sigma_rw",
         "log_tau2[1]",               # first time‐point
         sprintf("log_tau2[%d]", last_idx)  # last time‐point
      ),
      probs = c(0.1, 0.5, 0.9))

