# Joint Bayesian Factor Models
 A comprehensive repository with simple code for using Joint Bayesian Factor Models


**Folders**:

simulations/

This folder contains scripts and resources for running simulation studies to evaluate the performance and robustness of the Joint Bayesian Factor Models. Typical uses include generating synthetic datasets, applying the factor models under different scenarios, and benchmarking statistical properties.
Key files:

run_simulation.R: Main driver script for running simulation experiments across varying setups.
simulation_utils.R: Helper functions for data generation, results aggregation, and plotting.
results/: Contains output files and figures generated from simulation runs.
sparse_bayesian_infinite_factor_model/

Implements methods and code for sparse Bayesian infinite factor models, which allow for automatic determination of the number of latent factors and encourage sparsity in factor loadings.
Key files:

sbifm_model.stan: Stan model specification for the sparse Bayesian infinite factor model.
fit_sbifm.R: R script for fitting the model to data using Stan.
sbifm_utils.R: Utilities for preprocessing, postprocessing, and summarizing factor analysis results.
teb_far/

This directory provides code for the Two-Exchangeable Blocks Factor Analysis Regression (TEB-FAR) model, suitable for settings with structured latent factors and grouped data.
Key files:

teb_far_model.stan: Stan model file for TEB-FAR.
fit_teb_far.R: Main script to fit the TEB-FAR model to provided datasets.
teb_far_postprocess.R: Functions and scripts for extracting, visualizing, and interpreting TEB-FAR model results.
