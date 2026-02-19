library(survey)
library(pcalg)

#' Custom Survey-Weighted Conditional Independence Test
#' @param x, y, S Column indices for the variables to test (X _|_ Y | S)
#' @param suffStat A list containing the 'data' and 'weights'
svy_ci_test <- function(x, y, S, suffStat) {
  dat <- suffStat$data
  wts <- suffStat$weights
  
  # Create the survey design object
  des <- svydesign(ids = ~1, weights = ~wts, data = dat)
  
  # Get variable names
  var_x <- colnames(dat)[x]
  var_y <- colnames(dat)[y]
  
  # Build the regression formula: X ~ Y + S
  if (length(S) == 0) {
    form <- as.formula(paste(var_x, "~", var_y))
  } else {
    var_S <- colnames(dat)[S]
    form <- as.formula(paste(var_x, "~", var_y, "+", paste(var_S, collapse = " + ")))
  }
  
  # Run the survey-weighted regression
  # Note: If X is continuous, use gaussian. If binary, use quasibinomial.
  # For mixed data, you can dynamically check the class of dat[[var_x]] here.
  if (is.factor(dat[[var_x]]) && nlevels(dat[[var_x]]) == 2) {
    fit <- try(svyglm(form, design = des, family = quasibinomial()), silent = TRUE)
  } else {
    fit <- try(svyglm(form, design = des), silent = TRUE)
  }
  
  if (inherits(fit, "try-error")) return(1) # Return non-significant p-value if model fails
  
  # Extract the p-value for the coefficient of Y using a design-adjusted Wald test
  summ <- summary(fit)
  coefs <- summ$coefficients
  
  if (var_y %in% rownames(coefs)) {
    pval <- coefs[var_y, "Pr(>|t|)"]
    return(pval)
  } else {
    return(1) # Y was dropped (e.g., collinearity)
  }
}

#' Survey-Weighted Fast Causal Inference (FCI)
#' @param data A dataframe of the survey responses
#' @param weights A numeric vector of survey weights
#' @param alpha Significance level for CI tests
svy_fci <- function(data, weights, alpha = 0.05, ...) {
  # Package the data and weights into the suffStat object
  suffStat <- list(data = data, weights = weights)
  
  # Call your provided fci() or pcalg::fci() using our custom survey CI test
  res <- fci(suffStat = suffStat, 
             indepTest = svy_ci_test, 
             alpha = alpha, 
             labels = colnames(data), 
             p = ncol(data), 
             ...)
  return(res)
}