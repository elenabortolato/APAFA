# APAFA
**By E. Bortolato & A. Canale**


This repository contains the code to reproduce simulations, figures, and real data analysis results from the paper "Adaptive Partition Factor Analysis" by Elena Bortolato and Antonio Canale (https://arxiv.org/html/2410.18939v2)

The method is designed to estimate factor models, accommodating study-specific and shared components, finding adaptively the specific structure and latent dimensions.

---

## 📁  Repository
```
APAFA/
├── sampler.R # Core Gibbs sampler implementing the APAFA model
├── simulate_data.R # Example script that  runs a  simulation study
├── workflow.Rmd # explains in detail how to use the method
├── workflow.md # 
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
│ ├──immune_data.RDS
│ ├──imputation_APAFA.R
│ ├──imputation_TETRIS.R
│ ├──results/ # Results & figures for the immune analysis
│ └── results/ # Results & figures for the immune analysis
│
└── Supplementary/ # Scripts reproducing additional experiments and  checkst reported in the Supplementary Materials of the paper.
  ├── effect_standardize_factors.R # code for obtaining Figure SM1 and SM2 in the supplementary Materials
  └── convergence_uncertainty.R # code for obtaining Figure SM3 and SM14 in the supplementary Materials
  └── sensitivity.R # code for obtaining Figure SM4 and SM5 in the supplementary Materials
  └── identifiability_test.R  # code for obtaining Figure SM6 in the supplementary Materials
  
```


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

# ---- Optional packages ----
install.packages(c(
  "ggpubr", "reshape2", "patchwork"
))

# ---- Version checks ----
required_versions <- list(
  Rcpp        = "1.0.14",
  mvtnorm     = "1.3.1",
  MCMCpack    = "1.7.1",
  calculus    = "1.0.1",
  unbiasedmcmc = "0.3.0",
  pgdraw      = "1.1",
  ggplot2     = "3.5.1"
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
 
