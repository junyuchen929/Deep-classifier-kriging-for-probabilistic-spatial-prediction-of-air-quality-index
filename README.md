# Deep Classifier Kriging for Probabilistic Spatial Prediction of Air Quality Index

## Overview

This repository implements **Deep Classifier Kriging (DCK)** for probabilistic spatial prediction of the Air Quality Index (AQI), including:

- Univariate simulation studies  
- Bivariate simulation studies  
- Real data analysis for the entire United States  
- Regional (subset) real data analysis  


---

## Univariate Simulation

Code for univariate simulation is located in the `dck_uni` folder, covering the following scenarios:

- Gaussian (1600 samples)  
- Gaussian (3600 samples)  
- Non-Gaussian (1600 samples)  
- Non-Gaussian (3600 samples)  

### Script Order

- `data_generation.R`  
  Generates univariate spatial simulation data.

- `classification.R`  
  Discretizes the response variable into multiple classes.

- `dck_g1600.ipynb`  
  Trains the Deep Classifier Kriging model.

- `Gaussian_Kriging.R`  
  Fits the classical Kriging model.

- `comparison.R`  
  Compares classical Kriging and Deep Classifier Kriging.

- `dk_g1600.ipynb`  
  Runs the DeepKriging model for comparison.

---

## Bivariate Simulation

Code for bivariate simulation is located in the `dck_bi` folder, covering the following scenarios:

- Gaussian (1600 samples)  
- Gaussian (3600 samples)  
- Non-Gaussian (1600 samples)  
- Non-Gaussian (3600 samples)  

### Script Order

- `data_generation.R`  
  Generates bivariate spatial simulation data.

- `2d-to-1d_Projection.R`  
  Performs 2D-to-1D projection and discretization.

- `dck.ipynb`  
  Trains the Deep Classifier Kriging model.

- `Gaussian_Kriging.R`  
  Fits the classical Kriging model.

- `Compute-Comparison_results.R`  
  Compares classical Kriging and Deep Classifier Kriging.

---

## Real Data Analysis: Whole United States

Code for real data analysis on the full U.S. domain is located in the `realdata_whole` folder.

### Script Order

- `EDA.R`  
  Conducts exploratory data analysis of U.S. AQI data.

- `data_clean.ipynb`  
  Cleans and preprocesses AQI data.

- `2d-to-1d_Projection.R`  
  Performs 2D-to-1D projection and discretization.

- `embedding_matrix.py`  
  Generates spatial embedding matrices.

- `dck_train.py`  
  Trains the Deep Classifier Kriging model.

- `Gaussian_Kriging.R`  
  Fits the classical Kriging model.

- `Compute-Comparison_results.R`  
  Compares classical Kriging and Deep Classifier Kriging.

---

## Real Data Analysis: Regional Subsets

Code for regional (subset) analysis is located in the `realdata_subset` folder.

### Script Order

- `data_clean.ipynb`  
  Cleans and preprocesses AQI data for regional subsets.

- `projection_subset.R`  
  Performs 2D-to-1D projection and discretization for each subset.

- `embedding_matrix_subset.py`  
  Generates spatial embedding matrices for subsets.

- `generate_heatmaps.ipynb`  
  Trains Deep Classifier Kriging models for each subset and saves predicted class probability outputs.

- `Gmixture_subset.ipynb`  
  Generates samples from conditional Gaussian mixture models.

- `map_subset.R`  
  Visualizes predictive quantiles alongside observed data, together with maps of extreme AQI exceedance probabilities.
