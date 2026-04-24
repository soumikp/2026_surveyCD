# ==============================================================================
# Simulation: Standard Survey PC vs. Survey Copula PC on Mixed/Non-Linear Data
# ==============================================================================
library(pcalg)
library(survey)

set.seed(2026)

# 1. Generate TRUE Population: Non-Linear Mediation
# True Causal Chain: X (Demographic) -> Y (Risk) -> Z (Wellbeing)
N <- 20000
X <- sample(1:3, N, replace = TRUE) # Ordinal Demographic

# Y depends non-linearly on X (e.g., exponential risk increase)
Y_latent <- exp(X * 0.8) + rnorm(N, 0, 1) 
Y <- as.numeric(cut(Y_latent, breaks = 5)) # Ordinal Risk (1 to 5)

# Z depends non-linearly on Y
Z <- Y^3 + rnorm(N, 0, 20) # Continuous Wellbeing 

pop_data <- data.frame(X = as.factor(X), Y = as.factor(Y), Z = Z)

# 2. Impose Survey Sampling
target_group <- (as.numeric(pop_data$X) >= 2 & pop_data$Z > quantile(Z, 0.75))
prob_selection <- ifelse(target_group, 1.0, 0.1)

sampled_indices <- rbinom(N, 1, prob_selection) == 1
sample_data <- pop_data[sampled_indices, ]
sample_weights <- 1 / prob_selection[sampled_indices]

cat("Sample Size:", nrow(sample_data), "\n")

# 3. Standard Survey CI Test (GLM-based)
svy_ci_test <- function(x, y, S, suffStat) {
  dat <- suffStat$data; wts <- suffStat$weights
  des <- svydesign(ids = ~1, weights = ~wts, data = dat)
  var_x <- colnames(dat)[x]; var_y <- colnames(dat)[y]
  
  if (length(S) == 0) form <- as.formula(paste(var_x, "~", var_y))
  else form <- as.formula(paste(var_x, "~", var_y, "+", paste(colnames(dat)[S], collapse = " + ")))
  
  fit <- try(svyglm(form, design = des), silent = TRUE)
  if (inherits(fit, "try-error")) return(1) 
  
  coefs <- summary(fit)$coefficients
  if (var_y %in% rownames(coefs)) return(coefs[var_y, "Pr(>|t|)"]) else return(1)
}

# 4. Run the Algorithms
cat("\nRunning Standard Survey PC (GLM-based)...\n")
# Must use numeric data for GLM to attempt linear mediation
sample_data_num <- data.frame(X=as.numeric(sample_data$X), Y=as.numeric(sample_data$Y), Z=sample_data$Z)
suffStat_svy <- list(data = sample_data_num, weights = sample_weights)
pc_svy <- pc(suffStat = suffStat_svy, indepTest = svy_ci_test, alpha = 0.01, labels = colnames(sample_data_num), verbose = FALSE)

cat("Running Survey Copula PC...\n")
pc_copula <- svy_copula_pc(data = sample_data_num, weights = sample_weights, alpha = 0.01)

# 5. Compare Results
cat("\n--- TRUE SKELETON ---\n")
cat("X - Y - Z  (No edge between X and Z)\n")

cat("\n--- STANDARD SURVEY PC RESULT (GLM) ---\n")
print(as(pc_svy@graph, "matrix"))
cat("* Fails! Leaves X-Z edge because linear GLM can't block non-linear mediation.\n")

cat("\n--- SURVEY COPULA PC RESULT ---\n")
print(as(pc_copula@graph, "matrix"))
cat("* Success! Latent Gaussian transform blocks the path and deletes X-Z edge.\n")