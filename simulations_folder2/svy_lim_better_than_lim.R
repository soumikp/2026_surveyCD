set.seed(2026)
pacman::p_load(fastICA, clue)
# 2. Generate the TRUE Population (N = 20,000)
# True Causal Chain: X -> Y -> Z 
# IMPORTANT for LiNGAM: Data MUST be non-Gaussian (we use Uniform noise)
N <- 20000
X <- runif(N, -2, 2)
Y <- 1.5 * X + runif(N, -2, 2)
Z <- 1.5 * Y + runif(N, -2, 2)

pop_data <- data.frame(X, Y, Z)

# 3. Impose Extreme Stratified Survey Sampling (Selection Bias)
# We oversample people where both X and Z are highly positive.
# This artificially distorts the joint non-Gaussian distribution.
target_group <- (pop_data$X > 0.5 & pop_data$Z > 0.5)
prob_selection <- ifelse(target_group, 1.0, 0.05) 

sampled_indices <- rbinom(N, 1, prob_selection) == 1
sample_data <- pop_data[sampled_indices, ]

# Calculate survey weights
sample_weights <- 1 / prob_selection[sampled_indices]

cat("Sample Size:", nrow(sample_data), "\n")

# 4. Define the svy_lim wrapper (from the package we built)
svy_lim <- function(X, weights, B = 50, verbose = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  edge_frequencies <- matrix(0, p, p)
  prob_wts <- weights / sum(weights)
  successful_runs <- 0
  
  for(b in 1:B) {
    # Bootstrap a Pseudo-Population using the survey weights
    idx <- sample(1:n, size = n, replace = TRUE, prob = prob_wts)
    pseudo_X <- X[idx, ]
    
    # Run naive LiNGAM on the corrected pseudo-population
    res <- try(lingam(as.matrix(pseudo_X), verbose = verbose), silent = TRUE)
    
    if(!inherits(res, "try-error")) {
      adj_matrix <- t(res$Bpruned != 0) * 1 
      edge_frequencies <- edge_frequencies + adj_matrix
      successful_runs <- successful_runs + 1
    }
  }
  
  inclusion_probs <- edge_frequencies / successful_runs
  consensus_adj <- t((inclusion_probs > 0.5) * 1) # Transposed to match [Cause, Effect]
  rownames(consensus_adj) <- colnames(consensus_adj) <- colnames(X)
  
  return(consensus_adj)
}

# 5. Run Standard LiNGAM (Naive, Unweighted)
cat("\nRunning Standard LiNGAM (Unweighted)...\n")
lingam_naive <- lingam(as.matrix(sample_data), verbose = FALSE)
# Extract adjacency matrix (transpose so rows are causes, cols are effects)
naive_adj <- t(lingam_naive$Bpruned != 0) * 1
rownames(naive_adj) <- colnames(naive_adj) <- colnames(sample_data)

# 6. Run Survey-Weighted LiNGAM
cat("Running Survey-Weighted LiNGAM...\n")
svy_adj <- svy_lim(sample_data, sample_weights, B = 50, verbose = FALSE)

# 7. Compare the Results
cat("\n--- TRUE GRAPH SKELETON ---\n")
cat("X -> Y -> Z  (No edge between X and Z)\n")

cat("\n--- NAIVE LiNGAM RESULT (UNWEIGHTED) ---\n")
print(naive_adj)

cat("\n--- SURVEY-WEIGHTED LiNGAM RESULT ---\n")
print(svy_adj)