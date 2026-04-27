# 1. Load Required Libraries
library(pcalg)
library(survey)

set.seed(2026)

# 2. Generate the TRUE Population (N = 50,000)
N <- 50000
X <- rnorm(N, 0, 1)
Y <- 0.8 * X + rnorm(N, 0, 1)
Z <- 0.8 * Y + rnorm(N, 0, 1)

pop_data <- data.frame(X, Y, Z)

# 3. Impose Extreme Stratified Survey Sampling
# We heavily oversample a specific cross-section (e.g., High Demographic + High Wellbeing).
# 100% chance of selection if in the target group, 5% chance if not.
target_group <- (X > 0.5 & Z > 0.5)
prob_selection <- ifelse(target_group, 1.0, 0.05) 

sampled_indices <- rbinom(N, 1, prob_selection) == 1
sample_data <- pop_data[sampled_indices, ]

# Calculate survey weights (Inverse Probability of Selection: 1.0 or 20.0)
sample_weights <- 1 / prob_selection[sampled_indices]

cat("Sample Size:", nrow(sample_data), "\n")

# 4. Define the Survey-Weighted CI Test
svy_ci_test <- function(x, y, S, suffStat) {
  dat <- suffStat$data
  wts <- suffStat$weights
  des <- svydesign(ids = ~1, weights = ~wts, data = dat)
  
  var_x <- colnames(dat)[x]
  var_y <- colnames(dat)[y]
  
  if (length(S) == 0) {
    form <- as.formula(paste(var_x, "~", var_y))
  } else {
    var_S <- colnames(dat)[S]
    form <- as.formula(paste(var_x, "~", var_y, "+", paste(var_S, collapse = " + ")))
  }
  
  fit <- try(svyglm(form, design = des), silent = TRUE)
  if (inherits(fit, "try-error")) return(1) 
  
  coefs <- summary(fit)$coefficients
  if (var_y %in% rownames(coefs)) {
    return(coefs[var_y, "Pr(>|t|)"])
  } else {
    return(1)
  }
}

# 5. Run Standard FCI (Naive, Unweighted)
cat("\nRunning Standard FCI (Unweighted)...\n")
suffStat_naive <- list(C = cor(sample_data), n = nrow(sample_data))
fci_naive <- fci(suffStat = suffStat_naive, 
                 indepTest = gaussCItest, 
                 alpha = 0.05,  # Relaxed slightly to standard 0.05
                 labels = colnames(sample_data), 
                 p = ncol(sample_data),
                 verbose = FALSE)

# 6. Run Survey-Weighted FCI
cat("Running Survey-Weighted FCI...\n")
suffStat_svy <- list(data = sample_data, weights = sample_weights)
fci_svy <- fci(suffStat = suffStat_svy, 
               indepTest = svy_ci_test,   
               alpha = 0.05, 
               labels = colnames(sample_data), 
               p = ncol(sample_data),
               verbose = FALSE)

# 7. Compare the Adjacency Matrices
cat("\n--- TRUE GRAPH SKELETON ---\n")
cat("X - Y - Z  (No edge between X and Z)\n")

cat("\n--- NAIVE FCI RESULT (UNWEIGHTED) ---\n")
print(fci_naive@amat)

cat("\n--- SURVEY-WEIGHTED FCI RESULT ---\n")
print(fci_svy@amat)