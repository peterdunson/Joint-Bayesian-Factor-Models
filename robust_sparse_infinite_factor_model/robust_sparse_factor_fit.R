# robust_sparse_factor_fit.R
# --------------------------------------------------------
# Fits the robust sparse Bayesian infinite factor model (Student-t errors + MGPS prior) via Stan
# --------------------------------------------------------

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# -------------- USER PARAMETERS -------------------
# Paths
stan_file   <- "robust_sparse_infinite_factor_model.stan"  # Update as needed
sim_data_rds <- "sim_scen1_1000.rds"                      # Simulated data path
output_rds   <- "robust_sparse_factor_fit_scen1.rds"       # Output path

# Model/Hyperparameters
K        <- 5         # Truncation level (upper bound for factors)
nu       <- 3          # Degrees of freedom for Student-t
a1       <- 2          # MGP: delta_1 prior shape
a2       <- 3          # MGP: delta_{l>=2} prior shape
kappa    <- 2          # Local shrinkage
a_sigma  <- 2          # Error variance prior (Inverse Gamma)
b_sigma  <- 2          # Error variance prior (Inverse Gamma)

n_chains <- 4
n_iter   <- 3000
n_warmup <- 1500

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/simulations")

# -------------- LOAD DATA ------------------------
sim <- readRDS(sim_data_rds)
Y   <- sim$Y  # Assumed to be (n x p)
N   <- nrow(Y)
P   <- ncol(Y)

# -------------- PREPARE DATA FOR STAN -------------
stan_data <- list(
   N = N,
   P = P,
   K = K,
   Y = Y,
   nu = nu,
   a1 = a1,
   a2 = a2,
   kappa = kappa,
   a_sigma = a_sigma,
   b_sigma = b_sigma
)

# -------------- FIT MODEL -------------------------
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/robust_sparse_infinite_factor_model")
mod <- stan_model(stan_file)

fit <- sampling(
   object = mod,
   data = stan_data,
   chains = n_chains,
   iter = n_iter,
   warmup = n_warmup,
   seed = 42,
   control = list(adapt_delta = 0.99, max_treedepth = 15)
)

# -------------- SAVE & PRINT SUMMARY --------------
saveRDS(fit, output_rds)
print(fit, pars = c("sigma2", "tau"), probs = c(0.1, 0.5, 0.9))
