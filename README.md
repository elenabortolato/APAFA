# APAFA
**By E. Bortolato & A. Canale**
This repository contains the code to reproduce simulations, figures, and real data analysis results from the paper "Adaptive Partition Factor Analysis" by Elena Bortolato and Antonio Canale (https://arxiv.org/html/2410.18939v2)

The method is designed to estimate factor models, accommodating study-specific and shared components, finding adaptively the specific structure and latent dimensions.

---

## 📁 Repository Structure

APAFA/
├── sampler.R # Main Gibbs sampler / algorithm implementation
├── apa­fa_utils.R # Utility functions (data prep, post-processing, diagnostics)
├── examples/ # Example scripts & datasets
│ ├── simulate_data.R # contains one example of code used to run the simulation studies with APAFA
│ └── run_apa­fa_example.R
├── README.md # This file
└── Supplementary Materials/ # code to reproduce supplementary materials

- sampler.R contains the code for the Gibbs sampler for Gaussian data and binary data (N.B. relies on .cpp dependencies)
- simulation.zip contains the simulation results (as .RDS and .RData files) and the R code to reproduce Figure 4  in the manuscript.
- Real_data contains the code and data to reproduce the analysis and figures on the real data examples of Section 4 and the code for the out-of-sample predictive experiment reported in the Supplementary materials
- Supplementary contains short examples of workflow and the material to reproduce some extra simulations that are reported in the Supplementary Materials: (the effect of the prior hyperparameters choice, the effect of standardizing the factors, extra visualizations and uncertainty quantification summaries)

- In each of the subfolders there's a decription of the files contained


 
