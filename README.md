# APAFA
**By E. Bortolato & A. Canale**

This repository contains the code to reproduce simulations, figures, and real data analysis results from the paper "Adaptive Partition Factor Analysis" by Elena Bortolato and Antonio Canale (https://arxiv.org/html/2410.18939v2)

The method is designed to fit Bayesian (multi-study) factor models, accommodating study-specific and shared components, finding adaptively group-specific structure and latent dimensions - acknowledging potential model misspecification— and does so flexibly to accommodate a wide range of underlying patterns.

To reproduce the results of Section 3 of the paper, go to the [Simulations](Simulations/Simulation_results.md) workflow. 
To reproduce the results of the real data analysis of Section 4 of the paper, go to the [Birds workflow](Real_Data/birdsmd.md) or the [Immune workflow](Real_Data/immune_workflow.md) in the [Real_Data](Real_data) folder. 
Results contained in the Supplementary Materials are obtained by running the scritps in the [Supplementary](Supplementary) folder.

 

---


## 📁  Repository map
```
APAFA/
├── sampler.R # Core Gibbs sampler implementing the APAFA model
├── vignette.Rmd # explains in detail how to use the method
├── vignette.md # knitted file from  vignette.Rmd
│
├── Realdata/
│ ├── birds/ # Bird community study (Section 4.1)
│ │ ├── 0birds_folder_descriprion.md # description of the folder
│ │ ├── birdsmd.Rmd # Performs the full real-data analysis for the birds example
│ │ ├── birdsmd.md # output of the full real-data analysis for the birds example, Figure 6
│ │ └── data/ # data
│ │   ├── data.csv # Count data: species × locations/time 
│ │   ├── Ctree.tre # Phylogenetic tree of bird species (Newick/tre format)
│ │   ├── traits.csv # Species traits (e.g., mass, habitat)
│ │   ├── grid1000.csv # Location-level covariates for 1,000 locations (habitat, temperature, ...)
│ │   ├── grid10000.csv # Larger grid (10,000 locations) used for some robustness checks / predictions
│ │   ├── Ctree.tre # Phylogenetic tree of bird species (Newick/tre format) grid10000.csv 
│ │   └── birds_results.RData # results of the analysis.
│ │ 
│ └── immune/ # Immune response dataset (Section 4.2)
│   ├── 0immune_descripton_folder.md # description of the folder
│   ├── immune_data/immune_data.rda # dataset
|   ├── imputation_APAFA.R # code for the prediction exercise (figure SM20)
|   ├── imputation_TETRIS.R # code for the prediction exercise (figure SM21)
|   ├── immune_workflow.md # output of the analysis 
│   └── immune_workflow.Rmd # workflow to reproduce the analysis (Figures 7, 9)
│
├── Simulations/ # Code + results to reproduce simulation studies (Section 2)
│ ├── Simulation_results.Rmd # Reproduces the results of the simuation studies (Tables 1 and 2 and Figure 4 of the manuscript)
│ ├── simulate_APAFA.R # Example script that  runs a  simulation study (APAFA)
│ ├── simulate_TETRIS.R # Example script that  runs a  simulation study (TETRIS)
│ ├── example_out.R # Example script for producing the files in the "long" and "large" folders
│ │     # Figures 3 and 5, Figures SM9 and SM10 in the supplementary materials
│ ├── large/ # Results of the simulation studies for data in the "large" format (n<p)
│ ├── long/ # Results of the simulation studies for data in the "long" format (n>p)
│ └── Simulation_results_files/ # figures for the Simulation_results.Rmd
│
└── Supplementary/ # Scripts reproducing additional experiments and  checks reported in the Supplementary Materials of the paper.
  ├── beta_plot.R # reproduces Figure SM11 SM12 and SM13 in the supplementary Materials
  ├── effect_standardize_factors.R # reproduces Figure SM1 and SM2 in the supplementary Materials
  ├── convergence_uncertainty.R # reproduces Figure SM3 and SM14 in the supplementary Materials
  ├── covariance_uncertainty.R # reproduces Figure SM16, SM17, SM18 and SM19 in the supplementary Materials
  ├── sensitivity.R # reproduces Figure SM4 and SM5 in the supplementary Materials
  └── identifiability_test.R  # reproduces Figure SM6, SM7 and SM8 in the supplementary Materials
  
```


---

## ⚙️ Software Requirements

- **R version:** ≥ 4.4.1  
- **Operating System:** Tested on macOS and Windows

### Required R packages

Install all required packages with:

```r
install.packages(c(
  "Rcpp", "RcppEigen", "RcppArmadillo", "mvtnorm", "matrixStats", 
  "ggplot2", "cowplot", "coda", "tidyverse"
))

# ---- Optional packages ---- for reproducing figures in the supplementary materials
install.packages(c(
  "ggpubr", "reshape2", "patchwork", "rgl"
))

# ---- Version checks ----
required_versions <- list(
  Rcpp        = "1.0.14",
  RcppEigen = "0.3.4.0.2",
  RcppArmadillo = "15.0.2.2",
  mvtnorm     = "1.3.1",
  MCMCpack    = "1.7.1",
  calculus    = "1.0.1",
  pgdraw      = "1.1",
  ggplot2     = "3.5.1"
  matrixStats = "1.4.1",
  ggplot2=       "4.0.0",
  cowplot=     "1.1.3",  
  coda =  "0.19.4.1",
  tidyverse =   "2.0.0" 
)

check_versions <- function(pkgs) {
  for (pkg in names(pkgs)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(sprintf("Package '%s' is not installed.", pkg))
    } else {
      current <- as.character(utils::packageVersion(pkg))
      required <- pkgs[[pkg]]
      if (utils::compareVersion(current, required) < 0) {
        warning(sprintf(
          "Package '%s' version %s found, but %s or higher is required.",
          pkg, current, required
        ))
      } else {
        message(sprintf("✔ %s (version %s)", pkg, current))
      }
    }
  }
}

# Run version check
check_versions(required_versions)
```
 
