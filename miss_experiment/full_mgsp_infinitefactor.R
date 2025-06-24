# ─── Load & source infinitefactor edits ───────────────────
library(Rcpp)
source("~/Desktop/infinitefactor_edits/R/linearMGSP.R")
sourceCpp("~/Desktop/infinitefactor_edits/src/samplerBits.cpp", rebuild=TRUE)

# ─── Prepare your data ────────────────────────────────────
sim <- readRDS("~/Desktop/.../sim_scen2_1000.rds")
yX  <- scale(sim$Y, center=TRUE, scale=TRUE)  
# here yX[ ,1] is your outcome y, and its variance will be estimated

# ─── Run the full joint sampler ───────────────────────────
nrun <- 5000
burn <- 2500
fit_JBFM <- linearMGSP(
  X    = yX,
  nrun = nrun,
  burn = burn
)

# ─── What it does ─────────────────────────────────────────
# • Treats all P+1 columns (including y) identically, with their own residual σ²ₚ  
# • Samples that σ²_y in the same Gibbs/MH steps as the other ψ’s  
# • Returns omegaSamps = posterior samples of Ω = ΛΛᵀ + Σ (including Σ₁₁ = σ²_y)

# ─── Quick checks ─────────────────────────────────────────
table(fit_JBFM$numFactors)             # inferred K over iterations
Omega_hat <- apply(fit_JBFM$omegaSamps, c(1,2), mean)
