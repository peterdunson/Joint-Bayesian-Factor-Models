library(rstan)

# 1) Divergences + treedepth
rstan::check_hmc_diagnostics(fit_x)

# 2) Energy diagnostic (BFMI)
rstan::check_energy(fit_x)




sum_x <- summary(fit_x)$summary
# pick the top 3 parameters by Rhat
bad_pars <- names(sort(sum_x[,"Rhat"], decreasing=TRUE))[1:3]

library(bayesplot)
posterior_array <- as.array(fit_x)
mcmc_pairs(posterior_array,
           pars = bad_pars,
           off_diag_fun = "hex",
           diag_fun     = "hist")

