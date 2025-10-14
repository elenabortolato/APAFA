# APAFA
**By E. Bortolato & A. Canale**
This repository contains the code to reproduce simulations, figures, and real data analysis results from the paper "Adaptive Partition Factor Analysis" by Elena Bortolato and Antonio Canale (https://arxiv.org/html/2410.18939v2)

The method is designed to estimate factor models, accommodating study-specific and shared components, finding adaptively the specific structure and latent dimensions.

---

## 📁  Repository
```
APAFA/
├── sampler.R # Core Gibbs sampler implementing the APAFA model
├── simulate_data.R # Example script that generates toy data and runs a small simulation using APAFA
├── R/ # Reusable R functions (utilities, plotting, postprocessing)
│ ├── utils.R
│ ├── postprocessing.R
│ └── data_prep.R
│
├── simulations/ # Code + results to reproduce simulation studies (Section 2)
│ ├── run_simulations.R # Driver script for running the full set of simulation scenarios
│ ├── results/ # Cached posterior draws and summary tables (large files omitted from repo)
│ └── figures/ # Figures and tables generated from simulation results
│
├── Realdata/
│ ├── birds/ # Bird community study (Section 4.1)
│ │ ├── birds.R # Script performing the full real-data analysis for the birds example
│ │ ├── data.csv # Count data: species × locations/time (CSV)
│ │ ├── Ctree.tre # Phylogenetic tree (Newick/tre format)
│ │ ├── traits.csv # Species traits (e.g., mass, habitat)
│ │ ├── grid1000.csv # Location-level covariates for 1,000 locations (habitat, temperature, ...)
│ │ └── grid10000.csv # Larger grid (10,000 locations) used for some robustness checks / predictions
│ │
│ └── immune/ # Immune response dataset (other real-data example)
│ ├── run_immune.R # Script performing the immune response analysis
│ └── results/ # Results & figures for the immune analysis
│
└── Supplementary/ # Scripts reproducing supplementary experiments and robustness checks
  ├── identifiability_test.R
  └── sensitivity_analysis.R
```


- **`sampler.R`** → Core Gibbs sampler implementing the APAFA model.  
- **`simulate_data.R`** →  # contains one example of code used to run the simulation studies with APAFA
- **`simulations/`** → Code and results to reproduce simulation studies (Section 2 in the paper).  
- **`Realdata/birds/`** → Scripts for the bird community dataset analysis.
- **`Realdata/birds/birds.R`** → the code for performing all tha analysis in the  real data example of the paper of Section 4.1.
- **`Realdata/birds/data.csv`** →  count data of bird species observed in 200 locations over time
- **`Realdata/birds/Ctree.tre`** → phylogenetic tree of the species observed
- **`Realdata/birds/traits.csv`** → information on the species observed (mass, habitat)
- **`Realdata/birds/grid1000.csv`** → (and grid10000.csv) information about the locations (habitat, temperature...)
- **`Realdata/immune/`** → Scripts for the immune response dataset analysis.  
- **`Supplementary/`** → Reproduces additional experiments and robustness checksthat are reported in the Supplementary Materials of the paper.

---

## ⚙️ Software Requirements

- **R version:** ≥ 4.4.1  
- **Operating System:** Tested on macOS and Windows

### Required R packages

Install all required packages with:

```r
install.packages(c(
  "Rcpp", "RcppEigen", "mvtnorm", "matrixStats", 
  "ggplot2", "cowplot", "coda", "tidyverse"
))
```
 Optional packages (for plotting and diagnostics):
```r
install.packages(c(
  "ggpubr", "reshape2", "patchwork"
))
```


APAFA/
├── sampler.R # Main Gibbs sampler / algorithm implementation
├── simulate_data.R # contains one example of code used to run the simulation studies with APAFA
├── apa­fa_utils.R # Utility functions (data prep, post-processing, diagnostics)
├── examples/ # Example scripts & datasets
│ 
│ └── run_apa­fa_example.R
├── README.md # This file
└── Supplementary Materials/ # code to reproduce supplementary materials

- sampler.R contains the code for the Gibbs sampler for Gaussian data and binary data (N.B. relies on .cpp dependencies)
- simulation.zip contains the simulation results (as .RDS and .RData files) and the R code to reproduce Figure 4  in the manuscript.
- Real_data contains the code and data to reproduce the analysis and figures on the real data examples of Section 4 and the code for the out-of-sample predictive experiment reported in the Supplementary materials
- Supplementary contains short examples of workflow and the material to reproduce some extra simulations that are reported in the Supplementary Materials: (the effect of the prior hyperparameters choice, the effect of standardizing the factors, extra visualizations and uncertainty quantification summaries)



 
