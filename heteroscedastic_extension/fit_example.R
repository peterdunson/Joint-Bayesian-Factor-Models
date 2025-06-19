# 1) install & load the data

library(fda)

# 2) pull out the daily mean‐temperature curves
data("CanadianWeather", package="fda")
# this is a 365×35 matrix of daily temps
temps    <- CanadianWeather$dailyAv[,,"Temperature.C"]
stations <- colnames(temps)

# 3) reshape to long form
library(dplyr)
library(tidyr)
# 1) Pivot to long and add day
df_long <- as.data.frame(temps) %>%
   mutate(day = 1:365) %>%
   pivot_longer(-day, names_to = "station", values_to = "temp")

# 2) Choose 5 stations to keep
keep_stations <- unique(df_long$station)[1:5]

# 3) Filter down, then add time and a single subj index
df <- df_long %>%
   filter(station %in% keep_stations) %>%
   mutate(
      time = (day - 1) / 364,                                  # normalize to [0,1]
      subj = as.integer(factor(station, levels = keep_stations))  # 1:5
   )


# 5) build a spline basis (e.g. M=15)
library(splines)
M <- 15
B <- bs(df$time, df=M, intercept=TRUE)

# 6) pack stan_data
stan_data_weather <- list(
   N     = length(unique(df$station)),  # 35
   n_obs = nrow(df),                    # 35×365=12 775
   H     = 2,                           # choose how many factors
   M     = M,
   y     = df$temp,
   subj  = df$subj,
   B     = as.matrix(B)
)

str(stan_data_weather)






#fitting




# 1) Install & load rstan (if you haven’t already)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# 2) Compile the model (only needs to happen once)

setwd('/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/heteroscedastic_extension')

sm_weather <- stan_model("nhanes_het_ff.stan")

# 3) Run sampling
fit_weather <- sampling(
   sm_weather,
   data    = stan_data_weather,
   iter    = 2000,        # total iterations per chain
   warmup  = 1000,        # warmup (aka “burn‐in”) iterations
   chains  = 4,           # number of parallel chains
   control = list(adapt_delta = 0.9, max_treedepth=20)  
)

# 4) Quick convergence check & summaries
print(fit_weather, pars = c("sigma_beta","sigma_gamma"), probs = c(0.1,0.5,0.9))
# and to see the noise‐curve basis coeffs:
print(fit_weather, pars = "gamma", probs = c(0.1,0.5,0.9))

# make sure rstan is loaded
library(rstan)

# explicitly call rstan’s extract
post <- rstan::extract(fit_weather, permuted = TRUE)

# then you can do:
Beta_est  <- apply(post$Beta,  c(2,3), mean)   # H×M matrix of posterior means
gamma_est <-     apply(post$gamma, 2,   mean)  # length-M vector
# length‐M vector

# 6) Reconstruct the smooth functions on a fine grid
library(splines)
t_grid <- seq(0,1,length.out=365)
B_grid <- bs(t_grid, df = M, intercept = TRUE)

# Loading functions:
Lambda_est <- Beta_est %*% t(B_grid)   # H × 365
# Noise std‐dev:
tau_est    <- exp(0.5 * (B_grid %*% gamma_est))  # length‐365

# 7) Plot results
matplot(t_grid*365, t(Lambda_est), type="l", lwd=2,
        xlab="Day of year", ylab="Loading", main="Estimated Factors")
plot(t_grid*365, tau_est, type="l", lwd=2,
     xlab="Day of year", ylab="σ(t)", main="Estimated Noise SD")


scores <- apply(post$eta, c(2,3), mean)   # N×H
plot(scores[,1], scores[,2],
     xlab="Score on Factor 1", ylab="Score on Factor 2")
text(scores[,1], scores[,2], labels=keep_stations, pos=3, cex=0.7)




