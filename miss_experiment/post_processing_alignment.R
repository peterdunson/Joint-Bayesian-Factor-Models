library(GPArotation)
library(ggplot2)
library(bayesplot)

# 1) Load saved objects
obj_x        <- readRDS("fit_Xonly_scen2_scale_all.rds")
Lambda_x_hat <- obj_x$Lambda_hat
fit_x        <- obj_x$fit

obj_j        <- readRDS("fit_Joint_scen2_scale_all.rds")
Lambda_j_hat <- obj_j$Lambda_hat
fit_j        <- obj_j$fit

# 2) Varimax rotate the posterior means
Lambda_x_vm <- varimax(Lambda_x_hat)$loadings
Lambda_j_vm <- varimax(Lambda_j_hat)$loadings

# 3) Threshold small loadings
thresh <- 0.05
Lambda_x_vm[ abs(Lambda_x_vm) < thresh ] <- 0
Lambda_j_vm[ abs(Lambda_j_vm) < thresh ] <- 0

# 4) Quick heatmaps of the rotated + thresholded loadings
df_x <- reshape2::melt(Lambda_x_vm, varnames = c("Variable","Factor"), value.name = "Loading")
df_x$Model <- "X-only"
df_j <- reshape2::melt(Lambda_j_vm, varnames = c("Variable","Factor"), value.name = "Loading")
df_j$Model <- "Joint"
df_all <- rbind(df_x, df_j)

ggplot(df_all, aes(x=factor(Factor), y=factor(Variable), fill=Loading)) +
   geom_tile() +
   facet_wrap(~Model, ncol=2) +
   scale_fill_gradient2(low="darkgreen", mid="white", high="darkred", midpoint=0) +
   theme_minimal() +
   labs(title="Rotated & Thresholded Loadings", x="Factor", y="Variable")

# 5) Post-rotation diagnostics: density overlay of one loading
post_array_x <- as.array(fit_x)
post_array_j <- as.array(fit_j)

mcmc_dens_overlay(post_array_x, pars="Lambda_raw[1,1]") +
   ggtitle("Chain mixing for Lambda_raw[1,1] (X-only)")

mcmc_dens_overlay(post_array_j, pars="Lambda_raw[1,1]") +
   ggtitle("Chain mixing for Lambda_raw[1,1] (Joint)")

# 6) Full Rhat / n_eff tables
sum_x <- summary(fit_x)$summary
sum_j <- summary(fit_j)$summary

rhat_neff_x <- data.frame(Parameter=rownames(sum_x), Rhat=sum_x[,"Rhat"], n_eff=sum_x[,"n_eff"])
rhat_neff_j <- data.frame(Parameter=rownames(sum_j), Rhat=sum_j[,"Rhat"], n_eff=sum_j[,"n_eff"])

print(rhat_neff_x)
print(rhat_neff_j)
