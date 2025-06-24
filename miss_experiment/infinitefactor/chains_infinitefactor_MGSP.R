# run_mgps_chains_and_trace.R

# ─── 1) Libraries ─────────────────────────────────────────────────────────
library(Rcpp)
library(mvtnorm)
library(infinitefactor)  # for msf()
library(reshape2)
library(ggplot2)

# ─── 2) Source infinitefactor routines ────────────────────────────────────
source("~/Desktop/infinitefactor_edits/R/linearMGSP_fixedSigma1.R")
sourceCpp("~/Desktop/infinitefactor_edits/src/samplerBits.cpp", rebuild = TRUE)

# ─── 3) Load & prepare data ───────────────────────────────────────────────
sim  <- readRDS("~/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
Y    <- sim$Y                # n × (p+1)
yX   <- scale(Y, center = TRUE, scale = TRUE)

# ─── 4) Sampling settings ─────────────────────────────────────────────────
sigma2_y <- 1.0
nrun     <- 16000
burn     <- 8000
seeds    <- c(123, 456, 789, 101112)  # four different RNG seeds

# ─── 5) Run 4 independent chains ─────────────────────────────────────────
fits <- lapply(seeds, function(s) {
   set.seed(s)
   linearMGSP_fixedSigma1(
      X      = yX,
      Sigma1 = sigma2_y,
      nrun   = nrun,
      burn   = burn
   )
})

# ─── 6) Stack each chain’s lambdaSamps into an array & label chain ──────
chain_dfs <- lapply(seq_along(fits), function(i) {
   lambda_list <- fits[[i]]$lambdaSamps
   S <- length(lambda_list)
   p <- nrow(lambda_list[[1]])
   K <- ncol(lambda_list[[1]])
   arr <- array(unlist(lambda_list), dim = c(p, K, S),
                dimnames = list(Variable = paste0("V",1:p),
                                Factor   = paste0("F",1:K),
                                Iteration= as.character(1:S)))
   df <- melt(arr, varnames = c("Variable","Factor","Iteration"), value.name = "Loading")
   df$Iteration <- as.integer(df$Iteration)
   df$Chain     <- paste0("chain", i)
   df
})
df_all <- do.call(rbind, chain_dfs)

# ─── 7) Trace‐overlay for a few (j,k) pairs ───────────────────────────────
# choose three pairs:
pairs_to_plot <- subset(
   df_all,
   (Variable=="V5"&Factor=="F2") |
      (Variable=="V8"&Factor=="F3") |
      (Variable=="V1"&Factor=="F1")
)

p1 <- ggplot(pairs_to_plot, aes(x=Iteration, y=Loading, color=Chain)) +
   geom_line(alpha=1) +
   facet_wrap(~ Variable + Factor, scales="free_y", ncol=1) +
   theme_minimal() +
   labs(
      title = "Overlayed Chains for Selected Loadings",
      x     = "Iteration",
      y     = "Loading"
   )
print(p1)





# ─── Compute R̂ for MGSP loadings across 4 chains ────────────────────────
library(coda)

# ─── Full p×K matrix of R̂ ────────────────────────────────────────────────
sim    <- readRDS("~/Desktop/Joint-Bayesian-Factor-Models/simulations/sim_scen2_1000.rds")
p      <- nrow(sim$Lambda)
K      <- ncol(sim$Lambda)
rhat_mat <- matrix(NA_real_, nrow = p, ncol = K)

for (j in 1:p) {
   for (k in 1:K) {
      chains_jk <- lapply(fits, function(f) {
         mcmc( sapply(f$lambdaSamps, function(mat) if(ncol(mat)>=k) mat[j,k] else NA) )
      })
      rhat_mat[j,k] <- gelman.diag(mcmc.list(chains_jk), autoburnin = FALSE)$psrf[,"Point est."]
   }
}

print(rhat_mat)


library(pheatmap)

pheatmap(
   rhat_mat,
   cluster_rows = FALSE,
   cluster_cols = FALSE,
   main         = expression(hat(R) ~ "for each " ~ Lambda[j,k])
)



