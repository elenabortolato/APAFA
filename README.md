# APAFA

This repository contains the code to reproduce simulations, figures, and real data analysis results from the paper "Adaptive Partition Factor Analysis" by Elena Bortolato and Antonio Canale (https://arxiv.org/html/2410.18939v1)

## Organization 

- Real_data contains the code to reproduce the analysis and figures on the real data of Section 4 and the code for the out-of-sample predictive experiment of Section 4.2.1 (APAFA and TETRIS)
- Supplementary contains the material to reproduce some extra simulations that are reported in the Supplementary Materials 1:
  - sensitivity.R
  - effect_standardize_factors.R
  - convergence_and_uncertainty.R
  
- largep1000.R contains one example of code used to run the simulation studies with APAFA
- sampler.R contains the code for the Gibbs sampler for Gaussian data and binary data (N.B. relies on .cpp dependencies)
- simulation.zip contains the simulation results (as .RDS files) and the R code to reproduce Figure 4 in the manuscript.



 
