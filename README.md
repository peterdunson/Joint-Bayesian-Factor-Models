# Joint-Bayesian-Factor-Models

This repository contains simulation scripts and Stan implementations for sparse Bayesian infinite factor models and targeted empirical Bayes (TEB-FAR) factor regression models.

---

## Folder and File Structure

**Root Directory**

- **Joint-Bayesian-Factor-Models.Rproj**
  - **RStudio project file** for managing the repository.
- **README.md**
  - **This documentation file**.

---

**simulations/**

- Contains scripts for generating simulated datasets for different model testing scenarios.
    - **simulation_scen1.R**
      - Simulates **Scenario 1**: Outcome (*y*) is generated from a single latent factor.
      - Saves data as `.rds` for model input.
    - **simulation_scen2.R**
      - Simulates **Scenario 2**: Outcome (*y*) is generated similarly to the predictors, not just one factor.
      - Saves data as `.rds` for model input.
    - **simulation_scen3.R**
      - Simulates **Scenario 3**: Predictors are generated from a factor model, but outcome is only weakly related to predictors (not from factors).
      - Saves data as `.rds` for model input.

---

**sparse_bayesian_infinite_factor_model/**

- Contains scripts and Stan code for fitting the **baseline sparse Bayesian infinite factor model** (MGPS prior).
    - **fit_baseline_factor_model.R**
      - **R script for fitting the baseline factor model** to simulated data.
      - Loads data, preprocesses, runs Stan, and saves results.
    - **mgps_factor_model.stan**
      - **Stan model** for the baseline sparse Bayesian infinite factor model using the multiplicative gamma process shrinkage prior.

---

**teb_far/**

- Contains scripts and Stan code for fitting the **TEB-FAR model** (Targeted Empirical Bayes Factor Analysis Regression).
    - **fit_tebfar_factor_model.R**
      - **R script for fitting the TEB-FAR model** to simulated data.
      - Loads data, preprocesses, runs Stan, and saves results.
    - **tebfar_factor_model.stan**
      - **Stan model** for TEB-FAR, which fits a factor model with a fixed outcome variance.

---

## Usage

- Run any script in the **simulations/** folder to generate data for a given scenario.
- Use the relevant R script in **sparse_bayesian_infinite_factor_model/** or **teb_far/** to fit the corresponding model to your data.
- Output and results are saved in `.rds` format for downstream analysis.

---

**Feel free to use or adapt these scripts for your own simulation studies or model comparison tasks.**


---

## Credits

- **Sparse Bayesian Infinite Factor Model**: Based on the method from Bhattacharya & Dunson (2011), *Bayesian Analysis*.
- **Targeted Empirical Bayes Factor Regression (TEB-FAR)**: Based on the method from Palmer & Dunson (2025, arXiv preprint).
