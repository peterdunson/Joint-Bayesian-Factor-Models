# trace_lambda_eta_product.R

# 1) libs
library(rstan)       # for stanfit methods
library(bayesplot)   # for mcmc_trace()
library(dplyr)

# 2) load your 6 k fit
setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment")

obj <- readRDS("fit_Joint_scen2_scale_all_6k.rds")
fit <- if (inherits(obj, "stanfit")) obj else obj$fit

# 3) turn into an array [iter x chain x param]
post_array <- as.array(fit)  
# dims: (iterations) × (chains) × (parameters)

# 4) trace a few Λ, η, and ψ
pars_to_plot <- c("Lambda[1,1]", "eta[1,1]", "psi[1]")


# 5) build the Λ·η product for a single (j,i), say j=1 (variable), i=1 (obs)
j <- 1
i <- 1
K <- dim(post_array)[3] %>% names() %>%
   str_match_all("^Lambda\\[1,(\\d+)\\]$") %>%
   unlist() %>% as.integer() %>% max()




















library(bayesplot)   # for mcmc_trace
library(stringr)

# 1) grab your stanfit object
#    (replace this with however you load your stanfit of the joint model)
fit_j <- fit

# 2) turn it into an array: [iteration × chain × parameter]
post_array <- as.array(fit_j)


# 3) specify which j’s and which factor k
j_list <- c(5, 8, 20)
k      <- 2

# 4) pre‐allocate an array [iter × chain × length(j_list)]
iters  <- dim(post_array)[1]
chains <- dim(post_array)[2]
mu_arr <- array(
   NA_real_,
   dim = c(iters, chains, length(j_list)),
   dimnames = list(
      iteration = NULL,
      chain     = NULL,
      variable  = paste0("j", j_list)
   )
)

# 5) fill it in: μ_{i=j, j} = Λ[j, k] * η[j, k]
for (m in seq_along(j_list)) {
   j <- j_list[m]
   lam_nm <- sprintf("Lambda[%d,%d]", j, k)
   eta_nm <- sprintf("eta[%d,%d]",    j, k)
   lam_mat <- post_array[ , , lam_nm]
   eta_mat <- post_array[ , , eta_nm]
   mu_arr[ , , m] <- lam_mat * eta_mat
}





# 6) Manual trace‐plot of μ = Λ·η for j = 5,8,20, k = 2
# 6) Manual trace-plot of μ = Λ·η for j = 5,8,20, k = 2

library(reshape2)
library(ggplot2)

# Assume mu_arr is already in your workspace
df_mu <- melt(
   mu_arr,
   varnames   = c("Iteration","Chain","Variable"),
   value.name = "mu"
)

# Fix types
df_mu$Iteration <- as.integer(as.character(df_mu$Iteration))
df_mu$Chain     <- factor(df_mu$Chain, levels = unique(df_mu$Chain))

# Define your four chain colors
chain_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")

ggplot(df_mu, aes(x = Iteration, y = mu, color = Chain)) +
   geom_line(linewidth = 0.4) +
   facet_wrap(~ Variable, ncol = 1, scales = "free_y") +
   scale_color_manual(values = chain_cols) +
   theme_minimal(base_size = 12) +
   labs(
      title    = "Traceplots of μ[j,2] = Λ[j,2]·η[j,2]",
      subtitle = "for j = 5, 8, 20",
      x        = "Iteration",
      y        = expression(mu[j,2]),
      color    = "Chain"
   )









# ─── 7) Trace the raw Lambda[j,2] for j = 5,8,20 ─────────────────────────

# define j’s and factor k
j_list <- c(5, 8, 20)
k      <- 2

# collect the three parameter names
lam_pars <- sprintf("Lambda[%d,%d]", j_list, k)

# extract into an array [iter × chain × variable]
lam_arr <- post_array[,, lam_pars, drop = FALSE]

# melt to long form
df_lam <- reshape2::melt(
   lam_arr,
   varnames   = c("Iteration","Chain","Variable"),
   value.name = "lambda"
)
df_lam$Iteration <- as.integer(as.character(df_lam$Iteration))
df_lam$Chain     <- factor(df_lam$Chain, levels = unique(df_lam$Chain))

# plot
ggplot(df_lam, aes(x = Iteration, y = lambda, color = Chain)) +
   geom_line(linewidth = 0.4) +
   facet_wrap(~ Variable, ncol = 1, scales = "free_y") +
   scale_color_manual(values = chain_cols) +
   theme_minimal(base_size = 12) +
   labs(
      title    = expression(paste("Traces of ", Lambda[j,2], " for j = 5, 8, 20")),
      x        = "Iteration",
      y        = expression(Lambda[j,2]),
      color    = "Chain"
   )









# ─── 8) Trace the residual variances σ² for j = 5,8,20 ──────────────────────

# reuse j_list from above
# fits through stanfit post_array, so we need to extract sigma from post_array
# but in your joint model σ is parameter “psi[j]” on the precision scale,
# so variance is 1/psi[j].  If you saved sigmaSamps directly, adapt accordingly.

# here we get psi samples and invert
psi_pars <- sprintf("psi[%d]", j_list)

# extract into an array [iter × chain × variable]
psi_arr <- post_array[,, psi_pars, drop = FALSE]

# convert to variance = 1/precision
var_arr <- 1 / psi_arr

# melt to long form
df_var <- reshape2::melt(
   var_arr,
   varnames   = c("Iteration","Chain","Variable"),
   value.name = "sigma2"
)
df_var$Iteration <- as.integer(as.character(df_var$Iteration))
df_var$Chain     <- factor(df_var$Chain, levels = unique(df_var$Chain))

# plot
ggplot(df_var, aes(x = Iteration, y = sigma2, color = Chain)) +
   geom_line(linewidth = 0.4) +
   facet_wrap(~ Variable, ncol = 1, scales = "free_y") +
   scale_color_manual(values = chain_cols) +
   theme_minimal(base_size = 12) +
   labs(
      title    = expression(paste("Traceplots of residual variance ", sigma[j]^2, " for j = 5,8,20")),
      x        = "Iteration",
      y        = expression(sigma[j]^2),
      color    = "Chain"
   )
