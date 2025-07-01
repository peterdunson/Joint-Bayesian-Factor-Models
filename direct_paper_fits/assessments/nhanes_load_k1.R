# load_all_fits_nhanes.R

fit_dir <- "/Users/peterdunson/Desktop/Joint-Bayesian-Factor-Models/direct_paper_fits/storing_fit"

fit_MGSP        <- readRDS(file.path(fit_dir, "MGSP",         "fit_MGSP_NHANES1718_k1.rds"))
fit_HS          <- readRDS(file.path(fit_dir, "Horseshoe",    "fit_HS_NHANES1718_k1.rds"))
fit_SSL         <- readRDS(file.path(fit_dir, "spike_slab",   "fit_spikeslab_NHANES1718_k1.rds"))
fit_MASS        <- readRDS(file.path(fit_dir, "mass_nonlocal","fit_mass_nonlocal_NHANES1718_k1.rds"))
fit_RobustMGSP  <- readRDS(file.path(fit_dir, "robust",       "fit_robust_NHANES1718_k1.rds"))


