# swaoBN 0.1.0

## Initial Release

* Implements `swa_obn()`: survey-weighted, covariate-adjusted ordinal Bayesian
  network learning with expert constraints and two-step scoring.
* Greedy hill-climbing search (`search = "greedy"`) with multi-start random
  restarts.
* Exhaustive search (`search = "exhaust"`) for q ≤ 3 variables.
* Bootstrap uncertainty quantification for edge inclusion probabilities.
* Supports blacklists (forbidden edges) and whitelists (required edges).
* BIC and AIC information criteria with design-effect correction.
* Probit and logistic link functions.
