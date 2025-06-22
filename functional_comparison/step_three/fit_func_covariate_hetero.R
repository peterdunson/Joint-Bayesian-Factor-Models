library(dplyr)
library(purrr)
library(tidyr)
library(splines)
library(rstan)
library(refund)
data(content)

# --- 1. Prepare smoothed functional data ---
df_zlen <- content %>%
   select(id, agedays, zlen, zwei) %>%
   rename(subj = id, time = agedays, y = zlen, zwei = zwei) %>%
   filter(!is.na(y), !is.na(zwei))

t_grid <- seq(min(df_zlen$time), max(df_zlen$time), length = 20)

# Smoothing spline per subject, predict on grid
smoothed_list <- df_zlen %>%
   group_by(subj) %>%
   nest() %>%
   mutate(
      fit    = map(data, ~ smooth.spline(.x$time, .x$y, df = 8)),
      y_pred = map(fit, ~ predict(.x, x = t_grid)$y)
   )

smoothed_long <- expand.grid(
   subj = smoothed_list$subj,
   time = t_grid
) %>%
   arrange(subj, time)

smoothed_long$y <- unlist(smoothed_list$y_pred)
smoothed_long <- smoothed_long %>%
   mutate(subj = as.integer(factor(subj)))
smoothed_long <- smoothed_long %>%
   mutate(time = (time - min(time)) / (max(time) - min(time)))

# --- 2. Covariate matrix ---
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
Z_mat <- as.matrix(Z_df %>% select(sex, mean_zlen, mean_zwei))

# --- 3. Stan wrapper for heteroskedastic model ---
fit_func_mgps_cov_het <- function(df, Z_mat,
                                  stan_file     = "func_covariate_hetero.stan",
                                  H             = 3,
                                  M             = 15,
                                  iter          = 2000,
                                  warmup        = 1000,
                                  chains        = 4,
                                  adapt_delta   = 0.99,
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

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/functional_comparison/step_three")

# --- 4. Fit the model ---
fit3 <- fit_func_mgps_cov_het(
   df            = smoothed_long,
   Z_mat         = Z_mat,
   stan_file     = "func_covariate_hetero.stan",
   H             = 3,
   M             = 15,
   iter          = 2000,
   warmup        = 1000,
   chains        = 4,
   adapt_delta   = 0.99,
   maxtree_depth = 15
)

print(fit3, pars=c("sigma_beta", "sigma_gamma", paste0("delta[", 1:3, "]")), probs=c(0.1,0.5,0.9))

# --- 5. Save fit ---
saveRDS(fit3, file = "fit3_covariate_heteroskedastic_mgps.rds")
# To reload later:
# fit3 <- readRDS("fit3_covariate_heteroskedastic_mgps.rds")
