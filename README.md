# Joint Bayesian Factor Models

This repository contains simulation scripts and Stan implementations for sparse Bayesian infinite factor models, targeted empirical Bayes (TEB-FAR) regression, and robust factor models. It includes tools for model comparison and simulation-based benchmarking.

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

- Targeted Empirical Bayes Factor Analysis Regression (TEB-FAR) model.
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
        - Baseline factor model (Stan)
        - Lasso, Ridge, OLS, PCA+regression
    - **coverage_eval.R**
      - Evaluates 95% credible interval coverage for factor loadings (Bayesian models).
    - **frequentist_coverage_eval.R**
      - Evaluates 95% CI coverage for frequentist methods (OLS, Ridge, Lasso, PCA+regression).
    - **evaluate_lambda_coverage.R**
      - Additional coverage and symmetry analysis.
    - **post_process_results.R**
      - Utilities for post-processing fit outputs.

---

## Usage

- Generate data: Run any script in **simulations/** for a given scenario.
- Fit models: Use scripts in **sparse_bayesian_infinite_factor_model/**, **teb_far/**, or **robust_sparse_infinite_factor_model/**.
- Compare models: Run scripts in **model_comparisons/** to benchmark and analyze results.
- Output/results are saved as `.rds` files for downstream analysis.

---

**Feel free to use or adapt these scripts for your own simulation studies or model comparison tasks.**

---

## Credits

- **Sparse Bayesian Infinite Factor Model**: Bhattacharya & Dunson (2011, Biometrika) [link](https://github.com/david-dunson/sparse_bayesian_infinite_factor_models)
- **TEB-FAR**: Palmer & Dunson (2025, arXiv) [link](https://github.com/glennpalmer/TEB-FAR)
- **Robust Sparse Bayesian Infinite Factor Model**: Lee, Jo, Lee (2022, Comput Stat) [link](https://github.com/lee-jaejoon/robust-sparse-bayesian-infinite-factor-models/tree/main)


