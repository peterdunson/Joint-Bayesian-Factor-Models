# 1) Load libraries
library(bayesplot)
library(shinystan)    

getwd()

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/miss_experiment")

# 2) Read your saved fit
fit_j <- readRDS("fit_Joint_scen2_scale_all_6k.rds")$fit

# 3) Trace‐plots for a few representative parameters
#    (you can swap these for any pars you care about)
mcmc_trace(as.array(fit_j), pars = c("raw_Lambda[1,1]", "phi[1,1]", "delta[1]"))
mcmc_acf_bar(as.array(fit_j), pars = c("raw_Lambda[1,1]", "raw_Lambda[1,2]"))

# 4) Launch shinyStan for interactive exploration
launch_shinystan(fit_j)



# assume fit_j is your stanfit
library(dplyr)
library(tidyr)
library(stringr)

# get Rhat for every Lambda[p,k]
mon <- rstan::monitor(fit_j, pars="Lambda", print=FALSE)
rhat_df <- as.data.frame(mon[, "Rhat"], row.names=rownames(mon)) %>%
   rownames_to_column("param") %>%
   filter(str_detect(param, "^Lambda\\[")) %>%
   extract(param, into=c("p","k"),
           regex="Lambda\\[(\\d+),(\\d+)\\]", convert=TRUE)

# summarize by factor k
rhat_df %>%
   group_by(k) %>%
   summarize(
      mean_Rhat = mean(Rhat, na.rm=TRUE),
      pct_bad   = mean(Rhat > 1.1, na.rm=TRUE)
   ) %>%
   arrange(desc(mean_Rhat)) %>%
   print()








