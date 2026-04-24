# 1. Load Required Libraries
library(MASS)
library(survey)

set.seed(2026)

# 2. Generate the TRUE Population (N = 50,000)
# True Causal Chain: X -> Y -> Z
N <- 50000
X_cont <- rnorm(N, 0, 1)
Y_cont <- 1.5 * X_cont + rnorm(N, 0, 1)
Z_cont <- 1.5 * Y_cont + rnorm(N, 0, 1)

# Discretize into 3 ordinal categories (e.g., 1=Low, 2=Medium, 3=High)
X <- as.numeric(cut(X_cont, breaks = quantile(X_cont, c(0, 0.33, 0.66, 1)), include.lowest = TRUE))
Y <- as.numeric(cut(Y_cont, breaks = quantile(Y_cont, c(0, 0.33, 0.66, 1)), include.lowest = TRUE))
Z <- as.numeric(cut(Z_cont, breaks = quantile(Z_cont, c(0, 0.33, 0.66, 1)), include.lowest = TRUE))

pop_data <- data.frame(X = as.factor(X), Y = as.factor(Y), Z = as.factor(Z))

# 3. Impose Extreme Stratified Survey Sampling (Selection Bias)
# Oversample individuals with High Demographics (X=3) AND High Wellbeing (Z=3)
target_group <- (pop_data$X == 3 & pop_data$Z == 3)
prob_selection <- ifelse(target_group, 1.0, 0.05)

sampled_indices <- rbinom(N, 1, prob_selection) == 1
sample_data <- pop_data[sampled_indices, ]
sample_weights <- 1 / prob_selection[sampled_indices]

cat("Sample Size:", nrow(sample_data), "\n")

# 4. Define the Naive Scorer (Unweighted Standard BIC)
naive_score <- function(y_target, x_predictors) {
  if (is.null(x_predictors)) {
    fit <- polr(y_target ~ 1, method = "probit")
  } else {
    dat <- data.frame(Y = y_target, x_predictors)
    fit <- polr(Y ~ ., data = dat, method = "probit")
  }
  return(BIC(fit))
}

# 5. Define the Survey-Weighted Scorer (Pseudo-BIC)
svy_score <- function(y_target, x_predictors, weights) {
  if (is.null(x_predictors)) {
    dat <- data.frame(Y = y_target)
    form <- Y ~ 1
  } else {
    dat <- data.frame(Y = y_target, x_predictors)
    form <- as.formula(paste("Y ~", paste(colnames(x_predictors), collapse = " + ")))
  }
  
  des <- svydesign(ids = ~1, weights = ~weights, data = dat)
  fit <- svyolr(form, design = des, method = "probit")
  
  # Calculate Pseudo-BIC
  k <- length(fit$coefficients) + length(fit$zeta)
  n_eff <- sum(weights)
  pseudo_bic <- -2 * fit$deviance + k * log(n_eff)
  return(pseudo_bic)
}

# 6. Evaluate Graph 1: The TRUE Graph (X -> Y -> Z)
# Score = Score(X) + Score(Y|X) + Score(Z|Y)
cat("\nEvaluating Graph 1: X -> Y -> Z (True Graph, NO X->Z edge)\n")

score1_naive <- naive_score(sample_data$X, NULL) +
  naive_score(sample_data$Y, data.frame(X=sample_data$X)) +
  naive_score(sample_data$Z, data.frame(Y=sample_data$Y))

score1_svy <- svy_score(sample_data$X, NULL, sample_weights) +
  svy_score(sample_data$Y, data.frame(X=sample_data$X), sample_weights) +
  svy_score(sample_data$Z, data.frame(Y=sample_data$Y), sample_weights)

# 7. Evaluate Graph 2: The SPURIOUS Graph (X -> Y -> Z  AND  X -> Z)
# Score = Score(X) + Score(Y|X) + Score(Z | X, Y)
cat("Evaluating Graph 2: Spurious Graph (Has extra X->Z edge)\n")

score2_naive <- naive_score(sample_data$X, NULL) +
  naive_score(sample_data$Y, data.frame(X=sample_data$X)) +
  naive_score(sample_data$Z, data.frame(X=sample_data$X, Y=sample_data$Y))

score2_svy <- svy_score(sample_data$X, NULL, sample_weights) +
  svy_score(sample_data$Y, data.frame(X=sample_data$X), sample_weights) +
  svy_score(sample_data$Z, data.frame(X=sample_data$X, Y=sample_data$Y), sample_weights)

# 8. Compare the Results (Remember: LOWER BIC IS BETTER)
cat("\n--- NAIVE oBN RESULTS (Unweighted) ---\n")
cat("Score True Graph:    ", score1_naive, "\n")
cat("Score Spurious Graph:", score2_naive, "\n")
if(score2_naive < score1_naive) cat("=> FAILURE: Naive oBN prefers the Spurious Graph.\n")

cat("\n--- SURVEY-WEIGHTED oBN RESULTS ---\n")
cat("Score True Graph:    ", score1_svy, "\n")
cat("Score Spurious Graph:", score2_svy, "\n")
if(score1_svy < score2_svy) cat("=> SUCCESS: Survey-Weighted oBN prefers the True Graph.\n")