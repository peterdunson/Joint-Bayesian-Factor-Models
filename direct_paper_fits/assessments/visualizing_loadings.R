library(ggplot2)

# ---------- Collect loadings ----------
fit_list <- list(
   MGSP   = fit_MGSP,
   HS     = fit_HS,
   SSL    = fit_SSL,
   MASS   = fit_MASS,
   Robust = fit_RobustMGSP
)

loadings_df <- do.call(rbind, lapply(names(fit_list), function(mod) {
   lambda_hat <- as.numeric(fit_list[[mod]]$Lambda_hat)
   data.frame(
      Variable = seq_along(lambda_hat),
      Loading = lambda_hat,
      Model = mod
   )
}))

# ---------- 1. Faceted plot ----------
library(ggplot2)
p_faceted <- ggplot(loadings_df, aes(x = Variable, y = Loading)) +
   geom_col(fill = "#4287f5", alpha = 0.75) +
   facet_wrap(~Model, scales = "free_y") +
   labs(title = "K=1 Factor Loadings by Model (Faceted)", y = "Loading", x = "Variable") +
   theme_minimal(base_size = 15)

print(p_faceted)

# ---------- 2. Separate plots, one for each model ----------
plot_list <- lapply(unique(loadings_df$Model), function(mod) {
   ggplot(subset(loadings_df, Model == mod), aes(x = Variable, y = Loading)) +
      geom_col(fill = "#4287f5", alpha = 0.8) +
      labs(title = paste("K=1 Factor Loadings:", mod),
           y = "Loading",
           x = "Variable") +
      theme_minimal(base_size = 15)
})

# Show each plot (can use for saving as well)
for (p in plot_list) print(p)


