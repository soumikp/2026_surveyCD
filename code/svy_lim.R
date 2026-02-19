#' Survey-Weighted LiNGAM via Resampled Pseudo-Population
#' @param X A dataframe of continuous/mixed data
#' @param weights A numeric vector of survey weights
#' @param B Number of bootstrapped pseudo-populations to run (default 100 for stability)
svy_lim <- function(X, weights, B = 100, verbose = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Matrix to store the frequency of discovered edges across bootstraps
  edge_frequencies <- matrix(0, p, p)
  
  for(b in 1:B) {
    # 1. Create Pseudo-Population by sampling with replacement proportional to survey weights
    idx <- sample(1:n, size = n, replace = TRUE, prob = weights)
    pseudo_X <- X[idx, ]
    
    # 2. Run the provided lingam() function on the i.i.d. pseudo-population
    res <- lingam(as.matrix(pseudo_X), verbose = verbose)
    
    # 3. Record discovered edges (binary adjacency matrix)
    adj_matrix <- t(res$Bpruned != 0) * 1 
    edge_frequencies <- edge_frequencies + adj_matrix
  }
  
  # Average the edges to get inclusion probabilities (Consensus Graph approach)
  inclusion_probs <- edge_frequencies / B
  
  # Return final consensus graph (e.g., edges appearing in > 50% of resamples)
  consensus_adj <- (inclusion_probs > 0.5) * 1
  rownames(consensus_adj) <- colnames(consensus_adj) <- colnames(X)
  
  return(list(inclusion_probs = inclusion_probs, consensus_graph = consensus_adj))
}