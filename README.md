# Estimation and convergence rates in the distributional single index model

This repository contains code and replication material for the preprint

Fadoua Balabdaoui, Alexander Henzi, and Lukas Looser. "Estimation and convergence rates in the distributional single index model" (2023).

Functions implementing our methods are in `functions`, and data for replicating the simulations and the case study are in `data`.

For the simulations, `convergence_rates_experiments.R` generates a single run of the simulations, where the parameters can be chosen with the variable `id`. The collected simulation results are included as `simulation_results.zip` in the `data` folder, and can be processed with `simulation_results.R` after unzipping.

Code and results for the real data application are in the folder `application`. Results for the single index quantile regression are obtained with the code from [https://doi.org/10.1080/07350015.2021.2013245](https://doi.org/10.1080/07350015.2021.2013245).
