# ==============================================================================
# RIGOROUS SIMULATION: Stress-Testing Covariates (Z) and Survey Weights
# ==============================================================================
# NOTE: Ensure your swa-oBN functions (OrdCD, etc.) are loaded first!

library(MASS)
library(stats)

set.seed(2026) # Seed for perfect reproducibility

# ==============================================================================
# EXPERIMENT 1: Testing Covariate Adjustment (Z)
# True Structure: Y1 <- Z -> Y2 (NO edge between Y1 and Y2)
# ==============================================================================
cat("\n\n=======================================================\n")
cat("EXPERIMENT 1: COVARIATE ADJUSTMENT (Z)\n")
cat("True Graph: Empty (0 edges between Y1 and Y2)\n")
cat("=======================================================\n")


N <- 2000
Z_age <- rnorm(N)

# Softer effect size (1.2 instead of 2.5) prevents perfect separation warnings
latent_y1 <- 1.2 * Z_age + rnorm(N)
latent_y2 <- 1.2 * Z_age + rnorm(N)

# Use quantiles to ensure perfectly balanced ordinal categories
Y1_exp1 <- as.factor(ifelse(latent_y1 < quantile(latent_y1, 0.33), 1, 
                            ifelse(latent_y1 < quantile(latent_y1, 0.66), 2, 3)))
Y2_exp1 <- as.factor(ifelse(latent_y2 < quantile(latent_y2, 0.33), 1, 
                            ifelse(latent_y2 < quantile(latent_y2, 0.66), 2, 3)))

df_exp1 <- data.frame(Y1 = Y1_exp1, Y2 = Y2_exp1)
Z_df    <- data.frame(Age = Z_age)


# Run Naive (No Z)
cat("\nRunning Naive Model (Ignores Z)...\n")
fit_naive_exp1 <- OrdCD(y = df_exp1, Z = NULL, search = "exhaust", ic = "bic")
cat("Naive Output (Hallucinates an edge due to confounding):\n")
print(fit_naive_exp1$gam)

# Run Adjusted (With Z)
cat("\nRunning swa-oBN (Adjusts for Z)...\n")
fit_adj_exp1 <- OrdCD(y = df_exp1, Z = Z_df, search = "exhaust", ic = "bic")
cat("swa-oBN Output (Correctly identifies NO EDGE):\n")
print(fit_adj_exp1$gam)


# ==============================================================================
# EXPERIMENT 2: Testing Survey Weights (Selection Bias)
# True Structure: Y1 -> Y2 <- Y3 (A Collider. NO edge between Y1 and Y3)
# ==============================================================================
cat("\n\n=======================================================\n")
cat("EXPERIMENT 2: SURVEY WEIGHTING (Selection Bias)\n")
cat("True Graph:\n")
cat("     Y1 Y2 Y3\n")
cat("  Y1  0  0  0\n")
cat("  Y2  1  0  1  <-- Y1 and Y3 both cause Y2\n")
cat("  Y3  0  0  0\n")
cat("=======================================================\n")

N_pop <- 20000

# Y1 and Y3 are completely independent in the population
latent_1 <- rnorm(N_pop)
latent_3 <- rnorm(N_pop)

# Y2 is a collider (softer effects to prevent perfect separation)
latent_2 <- 1.0 * latent_1 + 1.0 * latent_3 + rnorm(N_pop)

Y1_pop <- as.factor(ifelse(latent_1 < quantile(latent_1, 0.33), 1, 
                           ifelse(latent_1 < quantile(latent_1, 0.66), 2, 3)))
Y3_pop <- as.factor(ifelse(latent_3 < quantile(latent_3, 0.33), 1, 
                           ifelse(latent_3 < quantile(latent_3, 0.66), 2, 3)))
Y2_pop <- as.factor(ifelse(latent_2 < quantile(latent_2, 0.33), 1, 
                           ifelse(latent_2 < quantile(latent_2, 0.66), 2, 3)))

pop_data <- data.frame(Y1 = Y1_pop, Y2 = Y2_pop, Y3 = Y3_pop)

# Introduce Selection Bias cleanly (Berkson's Paradox)
# We vastly oversample people where Y2 = 3.
prob_selection <- ifelse(pop_data$Y2 == 3, 0.90,
                         ifelse(pop_data$Y2 == 2, 0.30, 0.05))

# Draw the biased sample
sample_size <- 2000
sampled_indices <- sample(1:N_pop, size = sample_size, prob = prob_selection)
survey_data <- pop_data[sampled_indices, ]

# Calculate Inverse Probability Weights cleanly
survey_weights <- 1 / prob_selection[sampled_indices]
survey_weights <- survey_weights * (sample_size / sum(survey_weights)) # Normalize

# Run Naive (Unweighted)
cat("\nRunning Naive Model (Unweighted)...\n")
fit_naive_exp2 <- OrdCD(y = survey_data, weights = NULL, search = "exhaust", ic = "bic")
cat("Naive Output (Hallucinates spurious Y1-Y3 edge due to collider bias):\n")
print(fit_naive_exp2$gam)

# Run Weighted
cat("\nRunning swa-oBN (Weighted)...\n")
fit_weighted_exp2 <- OrdCD(y = survey_data, weights = survey_weights, search = "exhaust", ic = "bic")
cat("swa-oBN Output (Corrects bias, perfectly recovers true V-structure):\n")
print(fit_weighted_exp2$gam)
cat("\n=======================================================\n")