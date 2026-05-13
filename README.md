# swaoBN: Survey-Weighted, Covariate-Adjusted Ordinal Bayesian Network Learning

## Survey-Weighted Ordinal Causal Discovery: Application to Veteran Social Needs and Health Indicators

`swaoBN` is an R package for causal discovery from complex survey data with ordinal variables. The package implements a two-step scoring approach that separates skeleton discovery (using penalized pseudo-BIC) from edge orientation (using unpenalized pseudo-likelihood ratio), ensuring DAG selection consistency even when the number of ordinal levels varies across nodes.

## Authors

**Soumik Purkayastha** [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--3619--2804-brightgreen)](https://orcid.org/0000-0002-3619-2804)
Department of Biostatistics and Health Data Science, University of Pittsburgh School of Public Health;
Center for Healthcare Evaluation, Research, and Promotion, VA Pittsburgh Healthcare System
✉️ soumik@pitt.edu

**Lauren E. Russell** [![ORCID](https://img.shields.io/badge/ORCID-0000--0001--8703--9271-brightgreen)](https://orcid.org/0000-0001-8703-9271)
Office of Health Equity, Veterans Health Administration
✉️ Lauren.Russell3@va.gov

**David A. Frank** [![ORCID](https://img.shields.io/badge/ORCID-0009--0003--0282--787X-brightgreen)](https://orcid.org/0009-0003-0282-787X)
Center for Healthcare Evaluation, Research, and Promotion, VA Pittsburgh Healthcare System
✉️ David.Frank2@va.gov

**Leslie R. M. Hausmann** [![ORCID](https://img.shields.io/badge/ORCID-0000--0002--3909--6009-brightgreen)](https://orcid.org/0000-0002-3909-6009)
Center for Healthcare Evaluation, Research, and Promotion, VA Pittsburgh Healthcare System;
Department of Medicine, University of Pittsburgh School of Medicine
✉️ Leslie.Hausmann@va.gov

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("soumikp/swaoBN")
```

## Quick Start

```r
library(swaoBN)

# Load the synthetic toy dataset bundled with the package
data("swaoBN_toy")
toy_data <- swaoBN_toy$data
toy_weights <- swaoBN_toy$weights

# 1. Point estimate (Greedy hill-climbing search)
fit <- swa_obn(
  y = toy_data, 
  weights = toy_weights, 
  nstart = 10
)

# Inspect adjacency matrix (gam[i,j] = 1 means j -> i)
print(fit$gam)

# 2. Bootstrap edge inclusion probabilities
fit_boot <- swa_obn(
  y = toy_data, 
  weights = toy_weights, 
  nstart = 5,
  boot = 100 # Recommend 500+ for real applications
)

# Probabilities of edges
print(fit_boot$edge_probs)
```

## Features

- **Covariate Adjustment**: Adjust for continuous/categorical confounders (`Z` argument) that are not part of the DAG.
- **Survey Weights**: Appropriately scale log-likelihoods using Kish's effective sample size to account for informative sampling.
- **Expert Constraints**: Use `blacklist` and `whitelist` matrices to specify forbidden or required edges.
- **Two-Step Scoring**: Robust orientation of edges without finite-sample bias artifacts due to heterogeneous ordinal levels.

## Repository Structure

For reproducibility, this repository also houses the data analysis and simulation scripts from the accompanying manuscript. These materials can be found in the `inst/` subdirectory:

- `inst/analysis/`: R scripts for the main empirical application on the SHEP data.
- `inst/simulations/`: Code for the synthetic simulation studies (laptop and cluster versions).
- `inst/submission/`: Draft materials and cover letters.
