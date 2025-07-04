# APAFA

This repository contains the code to reproduce simulations, figures, and real data analysis results from the paper "Adaptive Partition Factor Analysis" by Elena Bortolato and Antonio Canale (https://arxiv.org/html/2410.18939v2)

## Organization 

- largep1000.R contains one example of code used to run the simulation studies with APAFA
- sampler.R contains the code for the Gibbs sampler for Gaussian data and binary data (N.B. relies on .cpp dependencies)
- simulation.zip contains the simulation results (as .RDS and .RData files) and the R code to reproduce Figure 4 in the manuscript.
- Real_data contains the code and data to reproduce the analysis and figures on the real data examples of Section 4 and the code for the out-of-sample predictive experiment reported in the Supplementary materials
- Supplementary contains short examples of workflow and the material to reproduce some extra simulations that are reported in the Supplementary Materials: (the effect of the prior hyperparameters choice, the effect of standardizing the factors, extra visualizations and uncertainty quantification summaries)

- In each of the subfolders there's a decription of the files contained


 
