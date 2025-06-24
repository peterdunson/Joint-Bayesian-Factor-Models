library(dplyr)
library(tidyr)
library(tibble)
library(rstan)

# 1) Load your 6k Stan fit RDS
obj   <- readRDS("fit_Joint_scen2_scale_all_6k.rds")
fit_j <- if (is.list(obj) && inherits(obj$fit, "stanfit")) obj$fit else obj

# 2) Summarize only the Lambda block
sumry <- summary(fit_j, pars = "Lambda")$summary

# 3) Tidy into a p×K matrix of R̂
rhat_mat <- as.data.frame(sumry) %>%
   rownames_to_column("param") %>%
   filter(grepl("^Lambda\\[", param)) %>%
   extract(param,
           into = c("j","k"),
           regex = "Lambda\\[(\\d+),(\\d+)\\]",
           convert = TRUE) %>%
   select(j, k, Rhat) %>%
   pivot_wider(
      names_from   = k,
      values_from  = Rhat,
      names_prefix = "F"
   ) %>%
   arrange(j) %>%
   column_to_rownames("j")

# 4) Inspect
print(rhat_mat)



library(dplyr)
library(tidyr)
library(tibble)
library(rstan)

# 1) Load your 6k Stan fit RDS
obj   <- readRDS("fit_Joint_scen2_scale_all_attempt.rds")
fit_j <- if (is.list(obj) && inherits(obj$fit, "stanfit")) obj$fit else obj

# 2) Summarize only the Lambda block
sumry <- summary(fit_j, pars = "Lambda")$summary

# 3) Tidy into a p×K matrix of R̂
rhat_mat <- as.data.frame(sumry) %>%
   rownames_to_column("param") %>%
   filter(grepl("^Lambda\\[", param)) %>%
   extract(param,
           into = c("j","k"),
           regex = "Lambda\\[(\\d+),(\\d+)\\]",
           convert = TRUE) %>%
   select(j, k, Rhat) %>%
   pivot_wider(
      names_from   = k,
      values_from  = Rhat,
      names_prefix = "F"
   ) %>%
   arrange(j) %>%
   column_to_rownames("j")

# 4) Inspect
print(rhat_mat)

