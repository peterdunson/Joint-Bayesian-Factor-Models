# Bayesian Factor Models

**UPDATE COMING SOON**













This repository contains simulation scripts and Stan implementations for **sparse Bayesian infinite factor models**, **targeted empirical bayes for supervised joint factor analysis (TEB-FAR)**, **robust and heavy-tailed factor models**, **horseshoe factor models**, and **continuous spike-and-slab lasso (SSL) factor models**. It includes tools for model comparison and simulation-based benchmarking.

---

## Folder and File Structure

**Root Directory**

- **Joint-Bayesian-Factor-Models.Rproj**
  - RStudio project file for managing the repository.

---

**simulations/**

- Scripts and data for generating simulated datasets under different scenarios.
    - **simulation_scen1.R**
      - Simulates **Scenario 1**: Outcome (*y*) generated from a single latent factor.
    - **simulation_scen2.R**
      - Simulates **Scenario 2**: Outcome (*y*) generated similarly to predictors (not just one factor).
    - **simulation_scen3.R**
      - Simulates **Scenario 3**: Predictors generated from a factor model; outcome weakly related (not from factors).
    - **sim_scenX_*.rds**
      - Simulated data files for use in model fitting and evaluation.

---

**sparse_bayesian_infinite_factor_model/**

- Baseline sparse Bayesian infinite factor model (MGPS prior).
    - **fit_baseline_factor_model.R**
      - R script to fit the baseline factor model to simulated data.
    - **mgps_factor_model.stan**
      - Stan code for the MGPS prior factor model.
    - **mgps_factor_model.rds**
      - Example Stan fit output (may be large).

---

**teb_far/**

- Targeted Empirical Bayes for Supervised Joint Factor Analysis (TEB-FAR) model.
    - **fit_tebfar_factor_model.R**
      - R script to fit the TEB-FAR model to simulated data.
    - **tebfar_factor_model.stan**
      - Stan code for TEB-FAR (factor model with fixed outcome variance).
    - **tebfar_factor_model.rds**
      - Example Stan fit output (may be large).

---

**robust_sparse_infinite_factor_model/**

- Robust sparse Bayesian infinite factor model (heavy-tailed errors).
    - **robust_sparse_factor_fit.R**
      - R script to fit the robust factor model to data.
    - **robust_sparse_infinite_factor_model.stan**
      - Stan code for the robust model (t-distributed errors).

---

**horseshoe_estimator/**

- **Horseshoe prior Bayesian factor analysis model.**
    - **fit_horseshoe_factor_model.R**
      - R script to fit the horseshoe factor model.
    - **horseshoe_factor_model.stan**
      - Stan code for the horseshoe prior factor model.
    - **horseshoe_factor_model.rds**
      - Example Stan fit output.

---

**continuous_spike_slab/**

- **Continuous spike-and-slab lasso (SSL) prior Bayesian factor model.**
    - **spike_slab_factor.R**
      - R script to fit the SSL factor model.
    - **spike_slab_factor_model.stan**
      - Stan code for the continuous spike-and-slab prior factor model.
    - **spike_slab_factor_model.rds**
      - Example Stan fit output.

---

**generalized-infinite-factorization-models/**

- Scripts for other infinite factorization models (from Schiavon et al.).
    - **Schiavon_SISGaussian.R**
      - Structured infinite factorization for Gaussian data.
    - **Schiavon_SIScovariatesRegression.R**
      - Structured infinite factorization for covariate regression.

---

**model_comparisons/**

- Scripts for comparing models and post-processing results.
    - **model_comparison_summary.R**
      - Loads simulation data and model fits, compares test-set MSE for:
        - TEB-FAR (Stan)
        - Baseline factor model (MGPS, Stan)
        - Horseshoe factor model (Stan)
        - Continuous spike-and-slab (SSL) factor model (Stan)
        - Lasso, Ridge, OLS, PCA+regression
    - **coverage_eval.R**
      - Evaluates 95% credible interval coverage for factor loadings (Bayesian models).
    - **frequentist_coverage_eval.R**
      - Evaluates 95% CI coverage for frequentist methods (OLS, Ridge, Lasso, PCA+regression).
    - **evaluate_lambda_coverage.R**
      - Additional coverage and symmetry analysis for loadings.
    - **post_process_results.R**
      - Utilities for post-processing fit outputs.
    - **horseshoe_eval.R**
      - Diagnostics and evaluation scripts for the horseshoe factor model.
    - **ssl_eval.R**
      - Diagnostics and evaluation scripts for the continuous spike-and-slab factor model.

---

## Usage

- Generate data: Run any script in **simulations/** for a given scenario.
- Fit models: Use scripts in **sparse_bayesian_infinite_factor_model/**, **teb_far/**, **horseshoe_estimator/**, **continuous_spike_slab/**, or **robust_sparse_infinite_factor_model/**.
- Compare models: Run scripts in **model_comparisons/** to benchmark and analyze results.
- Output/results are saved as `.rds` files for downstream analysis.

---

**Feel free to use or adapt these scripts for your own simulation studies or model comparison tasks.**

---

## Credits

- **Sparse Bayesian Infinite Factor Model:** Bhattacharya & Dunson (2011, Biometrika) [[paper link]](https://academic.oup.com/biomet/article/98/2/291/1745561)
- **TEB-FAR:** Palmer & Dunson (2025, arXiv) [[arXiv link]](https://arxiv.org/abs/2505.11351)
- **Robust Sparse Bayesian Infinite Factor Model:** Lee, Jo, Lee (2022, Comput Stat) [[paper link]](https://link.springer.com/article/10.1007/s00180-022-01208-5)
- **Horseshoe Factor Model:** Carvalho, Polson, Scott (2010, Biometrika) [[paper link]](https://www.jstor.org/stable/25734098)
- **Continuous Spike-and-Slab Lasso Factor Model:** Rocková (2018, Annals of Statistics) [[paper link]](https://projecteuclid.org/journals/annals-of-statistics/volume-46/issue-1/Bayesian-estimation-of-sparse-signals-with-a-continuous-spike-and/10.1214/17-AOS1554.full)

---
