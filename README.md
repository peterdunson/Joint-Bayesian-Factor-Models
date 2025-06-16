# Joint Bayesian Factor Models
 A comprehensive repository with simple code for using Joint Bayesian Factor Models


**Folders**

**sparse_bayesian_infinite_factor_model** - Includes code for “Sparse Bayesian Infinite Factor Models” (Bhattacharya & Dunson, 2011)

Files:

- mgps_factor_model.stan — Stan code for the Bayesian infinite factor model
- simulation_data.R — R script to simulate data with a sparse factor structure
- implementation.R — R script to fit the Stan model to data (minimal version, just fitting and saving)
- train_test_performance.R — R script for evaluating predictive performance using a train/test split; written in a general format for joint Bayesian factor models (expects real datasets with predictors and outcomes, not the simulation data here)
