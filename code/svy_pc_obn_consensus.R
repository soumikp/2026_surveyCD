#' Survey-Weighted PC + oBN Hybrid Consensus
svy_pc_obn_consensus <- function(data, weights, alpha = 0.05, link = "probit") {
  
  # STEP 1: Get the Skeleton / CPDAG using Survey-Weighted PC
  suffStat <- list(data = data, weights = weights)
  pc_res <- pc(suffStat = suffStat, 
               indepTest = svy_ci_test, 
               alpha = alpha, 
               labels = colnames(data), 
               p = ncol(data))
  
  # Extract the adjacency matrix from the PC result
  cpdag_adj <- as(pc_res@graph, "matrix")
  
  # STEP 2: Find Undetermined Edges (Markov Equivalence)
  # An edge is undetermined if cpdag_adj[i,j] == 1 AND cpdag_adj[j,i] == 1
  undirected_edges <- which(cpdag_adj == 1 & t(cpdag_adj) == 1, arr.ind = TRUE)
  
  # Keep only unique pairs (upper triangle)
  undirected_edges <- undirected_edges[undirected_edges[,1] < undirected_edges[,2], , drop = FALSE]
  
  # STEP 3: Orient edges using Survey-Weighted oBN
  final_adj <- cpdag_adj
  
  if (nrow(undirected_edges) > 0) {
    for (i in 1:nrow(undirected_edges)) {
      u <- undirected_edges[i, 1]
      v <- undirected_edges[i, 2]
      
      # Subset data for just these two variables
      sub_data <- data[, c(u, v)]
      
      # Run the exact survey-weighted bivariate oBN exhaust search
      # (Assuming svy_obn wrapper handles bivariate exhaustion)
      obn_res <- svy_obn(sub_data, weights = weights, search = "exhaust", link = link)
      
      # Update the final adjacency matrix based on oBN orientation
      biv_adj <- obn_res$gam
      final_adj[u, v] <- biv_adj[1, 2]
      final_adj[v, u] <- biv_adj[2, 1]
    }
  }
  
  return(list(initial_cpdag = cpdag_adj, final_dag = final_adj))
}