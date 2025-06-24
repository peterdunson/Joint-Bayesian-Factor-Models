# ─── 1) Libraries ─────────────────────────────────────────────────────────
library(ggplot2)
library(reshape2)

# ─── 2) Pull out lambdaSamps from your infinitefactor fit ────────────────
# (assumes you've already run `fit_TEBFAR <- linearMGSP_fixedSigma1(...)`)
lambda_list <- fit_TEBFAR$lambdaSamps  # list of S draws, each a p×K matrix

# ─── 3) Stack into an array [p × K × S] ───────────────────────────────────
S <- length(lambda_list)
p <- nrow(lambda_list[[1]])
K <- ncol(lambda_list[[1]])
arr <- array(
   unlist(lambda_list),
   dim = c(p, K, S),
   dimnames = list(
      Variable  = paste0("V", 1:p),
      Factor    = paste0("F", 1:K),
      Iteration = as.character(1:S)
   )
)

# ─── 4) Melt to long form ─────────────────────────────────────────────────
df_traces <- melt(
   arr,
   varnames   = c("Variable","Factor","Iteration"),
   value.name = "Loading"
)
df_traces$Iteration <- as.integer(df_traces$Iteration)

# ─── 5) (Optional) subset to a few J,K pairs ──────────────────────────────
# e.g. only V5–F2, V8–F3, V1–F1; otherwise comment out
# df_traces <- subset(df_traces,
#   (Variable=="V5" & Factor=="F2") |
#   (Variable=="V8" & Factor=="F3") |
#   (Variable=="V1" & Factor=="F1") 
# )

# ─── 6) Plot traces in a grid ─────────────────────────────────────────────
ggplot(df_traces, aes(x = Iteration, y = Loading)) +
   geom_line(alpha = 0.4) +
   facet_grid(Variable ~ Factor, scales = "free_y") +
   theme_minimal(base_size = 11) +
   labs(
      title = "MCMC Traces of λ₍j,k₎ from infinitefactor::linearMGSP",
      x     = "Iteration",
      y     = "Loading"
   )
