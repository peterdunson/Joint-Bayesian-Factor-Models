# 1) Load libraries
library(bayesplot)    # install.packages("bayesplot")
library(shinystan)    # install.packages("shinystan")

# 2) Read your saved fit
fit_j <- readRDS("fit_Joint_scen2_scale_all_6k.rds")$fit

# 3) Trace‐plots for a few representative parameters
#    (you can swap these for any pars you care about)
mcmc_trace(as.array(fit_j), pars = c("raw_Lambda[1,1]", "phi[1,1]", "delta[1]"))
mcmc_acf_bar(as.array(fit_j), pars = c("raw_Lambda[1,1]", "raw_Lambda[1,2]"))

# 4) Launch shinyStan for interactive exploration
launch_shinystan(fit_j)
