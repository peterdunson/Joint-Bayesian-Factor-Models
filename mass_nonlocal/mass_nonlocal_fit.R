# fit_and_save.R
library(rstan)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/mass_nonlocal")