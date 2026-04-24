# ==============================================================================
# File: svy_fci.R
# Description: Survey-Weighted Fast Causal Inference (FCI) Algorithm
# ==============================================================================
library(pcalg)
library(survey)

#' Custom Survey-Weighted Conditional Independence Test
#' Uses design-adjusted Wald tests to account for complex survey weights
svy_ci_test <- function(x, y, S, suffStat) {
  dat <- suffStat$data
  wts <- suffStat$weights
  
  # Create the survey design object
  des <- svydesign(ids = ~1, weights = ~wts, data = dat)
  
  var_x <- colnames(dat)[x]
  var_y <- colnames(dat)[y]
  
  # Build the regression formula: X ~ Y + S
  if (length(S) == 0) {
    form <- as.formula(paste(var_x, "~", var_y))
  } else {
    var_S <- colnames(dat)[S]
    form <- as.formula(paste(var_x, "~", var_y, "+", paste(var_S, collapse = " + ")))
  }
  
  # Check if variable X is binary to use logistic regression, otherwise linear
  is_binary <- length(unique(na.omit(dat[[var_x]]))) == 2
  
  if (is_binary) {
    fit <- try(svyglm(form, design = des, family = quasibinomial()), silent = TRUE)
  } else {
    fit <- try(svyglm(form, design = des), silent = TRUE)
  }
  
  if (inherits(fit, "try-error")) return(1) # Return non-significant p-value if model fails
  
  # Extract the p-value for the coefficient of Y
  coefs <- summary(fit)$coefficients
  if (var_y %in% rownames(coefs)) {
    return(coefs[var_y, "Pr(>|t|)"])
  } else {
    return(1) # Y was dropped due to collinearity
  }
}

#' Survey-Weighted FCI Wrapper
#' @param data Dataframe of survey responses
#' @param weights Numeric vector of survey weights
#' @param alpha Significance level for CI tests
svy_fci <- function(data, weights, alpha = 0.05, ...) {
  # Package data and weights into the sufficient statistic
  suffStat_svy <- list(data = data, weights = weights)
  
  # Run the standard FCI algorithm using the survey-aware CI test
  res <- fci(suffStat = suffStat_svy, 
             indepTest = svy_ci_test, 
             alpha = alpha, 
             labels = colnames(data), 
             p = ncol(data), 
             ...)
  return(res)
}