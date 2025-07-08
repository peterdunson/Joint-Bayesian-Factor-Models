# ---- 0. Setup ----
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/nhanes")

dat <- readRDS("nhanes_phthalates_adults.rds")

Y <- scale(dat, center = TRUE, scale = TRUE)
n <- nrow(Y)
p <- ncol(Y)
K <- 1

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/method_development")
library(rstan)
mod <- stan_model("mgsp_factor_model_discover.stan") # <-- your new Stan file!
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/ciprian_sim/method_development")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---- 1. Get Top 20% Off-Diagonal Correlations (by abs value) ----
cor_mat <- cor(Y)
off_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
cor_vals <- abs(cor_mat[off_idx])
n_top <- ceiling(length(cor_vals) * 0.05)
top_idx <- order(cor_vals, decreasing = TRUE)[1:n_top]
rows <- off_idx[top_idx, 1]
cols <- off_idx[top_idx, 2]

# ---- 2. Stan Data List ----
stan_data_j <- list(
   N = n,
   P = p,
   K = K,
   Y = Y,
   n_pen = n_top,
   row_idx = rows,
   col_idx = cols,
   pen_weight = 50  # <-- penalty strength (tune as needed)
)

# ---- 3. Run Stan ----
fit_j <- sampling(
   object       = mod,
   data         = stan_data_j,
   chains       = 4,
   iter         = 16000,
   warmup       = 8000,
   seed         = 19,
   init         = "random",
   init_r       = 2,
   control      = list(
      adapt_delta   = 0.99,
      max_treedepth = 15
   )
)

post_j       <- extract(fit_j)
Lambda_j_hat <- apply(post_j$Lambda, c(2,3), mean)

saveRDS(
   list(
      fit        = fit_j,
      posterior  = post_j,
      Lambda_hat = Lambda_j_hat
   ),
   file = "fit_Joint_NHANES1718_penalizedcorr_k1.rds"
)

sum_j      <- summary(fit_j)$summary
max_rhat_j <- max(sum_j[,"Rhat"],   na.rm = TRUE)
min_ess_j  <- min(sum_j[,"n_eff"], na.rm = TRUE)
bfmi_j     <- rstan::get_bfmi(fit_j)

cat("=== Penalized correlations MGSP model diagnostics ===\n")
cat(sprintf("  max R̂    = %.3f\n", max_rhat_j))
cat(sprintf("  min n_eff = %.0f\n",  min_ess_j))
cat(sprintf("  min BFMI  = %.3f\n", min(bfmi_j, na.rm = TRUE)))
