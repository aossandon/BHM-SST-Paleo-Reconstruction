# BHM-SST-Paleo-Reconstruction
This repository provides the dataset, R scripts, and STAN scripts for implementing the Bayesian Hierarchical Model (BHM) for estimating Paleo-SST reconstruction in the Equatorial Pacific, as proposed in _Ossandón et al. (2024)_.
<div align="center">
  <img src="Fig01.png" alt="" height="400"/>
  <p><em>Description: Schematic of the Bayesian Hierarchical Model (BHM) consisting of four layers: the Data Layer represents observed spatio-temporal SST data Y(s,t) and its mean μ(s,t) and covariance; Process Layer 1 models the temporal variation of μ; Process Layer 2 handles spatial modeling of covariance parameters and temporal regression coefficients; Process Layer 3 models spatially varying regression coefficients using Gaussian kernels; and the Prior Layer defines the prior distribution of hyperparameters.</em></p>
</div>

## Dataset
This dataset is used as input for implementing a space-time Bayesian hierarchical model (BHM) to reconstruct annual Sea Surface Temperature (SST) during the Holocene over a large domain based on SST at limited proxy locations in the equatorial (10°N-10°S) Pacific. The dataset consists of 2° gridded annual mean SST between 1854 and 2014, the paleo proxy SST data from 28 locations, locations of 15 spatial knots, and the location of points (150) used to calibrate the model.

## Scripts
### R
- __Library.R__: Install and load all the packages required. Also, it contains all the personal functions created for this implementation. 
- __Bayesian_model_15knots.R__: Calibrate the BHM for the contemporary SST anomaly data.
- Get_forecast_2021.R: Generate the ensemble streamflow forecast for the peak monsoon season 2021. Before running this script, the script Par_bayesian_model_gamma_best_mod_kLT.rds must be run for 1 to 5-days lead time. 
### STAN
- bayesian_model_gamma_2step_diff_st_2_best_kdlt.STAN: Contain stan code for the model structure of the best candidate BHMC for each lead time. Where k indicates the lead time.
