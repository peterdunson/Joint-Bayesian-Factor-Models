library(dplyr)
library(purrr)
library(tidyr)
library(splines)
library(rstan)

# 1. Prepare the smoothed functional data (as before)
df_zlen <- content %>%
   select(id, agedays, zlen, zwei) %>%           # Include zwei here
   rename(subj = id, time = agedays, y = zlen) %>%
   filter(!is.na(y), !is.na(zwei))

t_grid <- seq(min(df_zlen$time), max(df_zlen$time), length = 20)

# Fit per-subject smoothing spline and predict on grid
smoothed_list <- df_zlen %>%
   group_by(subj) %>%
   nest() %>%
   mutate(
      fit    = map(data, ~ smooth.spline(.x$time, .x$y, df = 8)),
      y_pred = map(fit, ~ predict(.x, x = t_grid)$y)
   )

# Create all combinations of subj and time
smoothed_long <- expand.grid(
   subj = smoothed_list$subj,
   time = t_grid
) %>%
   arrange(subj, time)

# Flatten the predicted values (row-wise) to fill y
smoothed_long$y <- unlist(smoothed_list$y_pred)

# Relabel subj to consecutive 1...N
smoothed_long <- smoothed_long %>%
   mutate(subj = as.integer(factor(subj)))

# Normalize time to [0,1]
smoothed_long <- smoothed_long %>%
   mutate(time = (time - min(time)) / (max(time) - min(time)))

# --- Construct covariate matrix: sex, mean zlen, mean zwei for each subject ---
Z_df <- content %>%
   group_by(id) %>%
   summarize(
      sex        = first(ma1fe0),
      mean_zlen  = mean(zlen, na.rm = TRUE),
      mean_zwei  = mean(zwei, na.rm = TRUE),
      .groups    = "drop"
   ) %>%
   mutate(subj = as.integer(factor(id))) %>%
   arrange(subj)

# The covariate matrix for Stan input (no id or subj column)
Z_mat <- as.matrix(Z_df %>% select(sex, mean_zlen, mean_zwei))

# --- Stan wrapper for covariate-adjusted functional factor model ---
fit_func_mgps_covariate <- function(df, Z_mat,
                                    stan_file   = "func_covariate_homo.stan",
                                    H           = 3,
                                    M           = 15,
                                    iter        = 2000,
                                    warmup      = 1000,
                                    chains      = 4,
                                    adapt_delta = 0.99,
                                    maxtree_depth = 15) {
   N <- length(unique(df$subj))
   J <- ncol(Z_mat)
   B <- bs(df$time, df = M, intercept = TRUE)
   stan_data <- list(
      N     = N,
      n_obs = nrow(df),
      J     = J,
      H     = H,
      M     = M,
      y     = df$y,
      subj  = df$subj,
      Z     = as.matrix(Z_mat),
      B     = as.matrix(B)
   )
   sm  <- stan_model(stan_file)
   fit <- sampling(
      sm,
      data    = stan_data,
      iter    = iter,
      warmup  = warmup,
      chains  = chains,
      control = list(adapt_delta = adapt_delta, max_treedepth = maxtree_depth)
   )
   return(fit)
}

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/functional_comparison/step_two")


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# --- Fit the model ---
fit2 <- fit_func_mgps_covariate(
   df            = smoothed_long,
   Z_mat         = Z_mat,
   stan_file     = "func_covariate_homo.stan",
   H             = 3,
   M             = 15,
   iter          = 2000,
   warmup        = 1000,
   chains        = 4,
   adapt_delta   = 0.99,
   maxtree_depth = 15
)

print(fit2, pars = c("sigma_y", paste0("delta[", 1:3, "]")), probs = c(0.1, 0.5, 0.9))

saveRDS(fit2, file = "fit2_covariate_adjusted_mgps_USING_20.rds")

#fit2 <- readRDS("fit2_covariate_adjusted_mgps.rds")


