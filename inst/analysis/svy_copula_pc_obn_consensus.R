# ==============================================================================
# File: svy_copula_pc.R
# Description: Survey-Weighted Copula PC Algorithm for Mixed Data
# ==============================================================================
library(pcalg)

#' Weighted Empirical Gaussian Copula Transform
#' Maps any continuous/ordinal variable to a standard normal distribution 
#' using its weighted empirical CDF.
svy_gaussian_copula_transform <- function(x, weights) {
  n <- length(x)
  
  # Get sorting order
  ord <- order(x)
  x_sorted <- x[ord]
  w_sorted <- weights[ord]
  
  # Calculate weighted CDF
  w_cdf <- cumsum(w_sorted) / sum(w_sorted)
  
  # Shrink slightly away from 1 to avoid qnorm(1) = Inf
  w_cdf <- w_cdf * (n / (n + 1))
  
  # Map back to original order
  cdf_orig <- numeric(n)
  cdf_orig[ord] <- w_cdf
  
  # Transform to Standard Normal
  z <- qnorm(cdf_orig)
  return(z)
}

#' Survey-Weighted Copula PC Algorithm
svy_copula_pc <- function(data, weights, alpha = 0.05, verbose = FALSE) {
  p <- ncol(data)
  n <- nrow(data)
  
  # 1. Transform all mixed variables to latent Gaussian space
  Z_matrix <- matrix(NA, nrow = n, ncol = p)
  for (j in 1:p) {
    # Convert factors/ordinals to numeric ranks for the CDF
    num_x <- as.numeric(data[[j]])
    Z_matrix[, j] <- svy_gaussian_copula_transform(num_x, weights)
  }
  
  # 2. Compute the Survey-Weighted Covariance/Correlation Matrix
  # Using the latent Gaussian variables
  wt_cov <- cov.wt(Z_matrix, wt = weights, cor = TRUE)
  Z_cor <- wt_cov$cor
  
  # 3. Run PC Algorithm using the exact Gaussian CI test on the latent space
  suffStat_copula <- list(C = Z_cor, n = n)
  
  pc_res <- pc(suffStat = suffStat_copula, 
               indepTest = gaussCItest, 
               alpha = alpha, 
               labels = colnames(data), 
               p = p,
               verbose = verbose)
  
  return(pc_res)
}


# ==============================================================================
# File: svy_copula_hybrid.R
# Description: Survey Copula PC + oBN Hybrid Consensus 
# ==============================================================================
library(pcalg)

# Ensure dependencies are available
# source("svy_copula_pc.R")
# source("svy_obn.R")

#' Survey-Weighted Copula PC + oBN Hybrid Consensus
svy_copula_hybrid <- function(data, weights, alpha = 0.05, link = "probit") {
  
  # STEP 1: Get the Skeleton using Survey-Weighted Copula PC
  pc_res <- svy_copula_pc(data = data, weights = weights, alpha = alpha)
  cpdag_adj <- as(pc_res@graph, "matrix")
  final_adj <- cpdag_adj
  
  # STEP 2: Find Undetermined Edges (Markov Equivalence)
  undirected_edges <- which(cpdag_adj == 1 & t(cpdag_adj) == 1, arr.ind = TRUE)
  undirected_edges <- undirected_edges[undirected_edges[,1] < undirected_edges[,2], , drop = FALSE]
  
  # STEP 3: Orient ambiguous edges using Survey-Weighted oBN
  if (nrow(undirected_edges) > 0) {
    for (i in 1:nrow(undirected_edges)) {
      u <- undirected_edges[i, 1]
      v <- undirected_edges[i, 2]
      
      sub_data <- data[, c(u, v)]
      
      # Run survey-weighted bivariate oBN exhaust search
      obn_res <- svy_obn(sub_data, weights = weights, ic = "bic", link = link, nstart = 1, maxit = 5)
      
      biv_adj <- obn_res$gam
      
      # Lock in the direction found by the ordinal scorer
      final_adj[u, v] <- biv_adj[1, 2]
      final_adj[v, u] <- biv_adj[2, 1]
    }
  }
  
  return(list(initial_cpdag = cpdag_adj, final_dag = final_adj))
}


# ==============================================================================
# File: svy_copula_hybrid.R
# Description: Survey Copula PC + oBN Hybrid Consensus 
# ==============================================================================
library(pcalg)

# Ensure dependencies are available
# source("svy_copula_pc.R")
# source("svy_obn.R")

#' Survey-Weighted Copula PC + oBN Hybrid Consensus
svy_copula_hybrid <- function(data, weights, alpha = 0.05, link = "probit") {
  
  # STEP 1: Get the Skeleton using Survey-Weighted Copula PC
  pc_res <- svy_copula_pc(data = data, weights = weights, alpha = alpha)
  cpdag_adj <- as(pc_res@graph, "matrix")
  final_adj <- cpdag_adj
  
  # STEP 2: Find Undetermined Edges (Markov Equivalence)
  undirected_edges <- which(cpdag_adj == 1 & t(cpdag_adj) == 1, arr.ind = TRUE)
  undirected_edges <- undirected_edges[undirected_edges[,1] < undirected_edges[,2], , drop = FALSE]
  
  # STEP 3: Orient ambiguous edges using Survey-Weighted oBN
  if (nrow(undirected_edges) > 0) {
    for (i in 1:nrow(undirected_edges)) {
      u <- undirected_edges[i, 1]
      v <- undirected_edges[i, 2]
      
      sub_data <- data[, c(u, v)]
      
      # Run survey-weighted bivariate oBN exhaust search
      obn_res <- svy_obn(sub_data, weights = weights, ic = "bic", link = link, nstart = 1, maxit = 5)
      
      biv_adj <- obn_res$gam
      
      # Lock in the direction found by the ordinal scorer
      final_adj[u, v] <- biv_adj[1, 2]
      final_adj[v, u] <- biv_adj[2, 1]
    }
  }
  
  return(list(initial_cpdag = cpdag_adj, final_dag = final_adj))
}