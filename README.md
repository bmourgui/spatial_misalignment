# Area-to-point spatial misalignment in SDMs

This repository contains the R-scripts to reproduce the analysis presented in the paper:

> B. Mourguiart, M. Chevalier, M. Marzloff, N. Caill-Milly, K. Mengersen, B. Liquet. Dealing with area-to-point spatial misalignment in species distribution models. _Ecography_, 2024(5), e07104, doi: 10.1111/ecog.07104

The study uses a virtual species approach to investigate the effects of area-to-point spatial misalignment (i.e. a resolution or a grain mismatch between coarse-grain covariates and fine-grain response variables) on the explanatory and predictive performance of three species distribution models (SDMs).

Running the "make.R" file allows you to run the entire workflow from simulating the data to producing the results presented in the manuscript.

Description of the files in the workflow:

  * 01_simulate_data.R generates virtual environmental and species presence-absence data under different scenarios with area-to-point spatial misalignment (i.e., environmental data with coarser resolution than presence-absence data);
  * 02_run_xxx.R fits the three models: the GLM (generalised linear model) that does not account for spatial misalignment, the spGLM (spatial GLM) that estimates residual spatial autocorrelation not captured by the environmental covariates, and the BEM (Berkson error model) that accounts for area-to-point spatial misalignment;
  * 03_compute_perfmetrics.R calculates performance metrics (e.g., AUC) to assess the explanatory and predictive performance of the three models over different levels of misalignment and spatial heterogeneity;
  * 04_analyze_perfmetrics.R computes summary statistics of the performance metrics;
  * 05_make_outputs.R generates the figures presented in the manuscript and supplementary materials.
