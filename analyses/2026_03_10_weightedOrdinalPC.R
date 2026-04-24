################################################################################
# Weighted PC-Stable Algorithm for Ordinal Data
#
# Modifications from bnlearn::pc.stable:
#   1. Input is a matrix/data.frame of ordinal (integer-coded) columns
#   2. A vector of survey weights is incorporated into the CI tests
#
# Approach: Use nonparanormal weighted correlation matrix, then use partial
# correlations derived from it for conditional independence testing.
################################################################################

# ------------------------------------------------------------------------------
# 1. Weighted Gaussian copula correlation matrix
#    (Semiparametric normal-scores approach)
# ------------------------------------------------------------------------------
weighted_copula_cormat <- function(data, weights) {
  # data: n x p matrix of ordinal integer-coded columns
  # weights: length-n vector of survey weights
  #
  # For each column, compute weighted ECDF -> qnorm to get normal scores,
  # then compute a weighted Pearson correlation on the scores.
  
  n <- nrow(data)
  p <- ncol(data)
  varnames <- colnames(data)
  if (is.null(varnames)) varnames <- paste0("V", seq_len(p))
  
  # Step 1: Transform each column to normal scores via weighted ECDF
  Z <- matrix(0, n, p)
  colnames(Z) <- varnames
  
  for (j in seq_len(p)) {
    Z[, j] <- weighted_normal_scores(data[, j], weights)
  }
  
  # Step 2: Weighted Pearson correlation on normal scores
  R <- weighted_pearson_matrix(Z, weights)
  
  # Project to nearest positive-definite matrix if needed
  eig <- eigen(R, symmetric = TRUE)
  if (any(eig$values < 0)) {
    eig$values[eig$values < 0] <- 1e-8
    R <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    d <- sqrt(diag(R))
    R <- R / (d %o% d)
  }
  
  rownames(R) <- colnames(R) <- varnames
  R
}

# ------------------------------------------------------------------------------
# 2. Weighted normal scores for a single ordinal vector
#    Weighted ECDF -> Winsorize -> qnorm
# ------------------------------------------------------------------------------
weighted_normal_scores <- function(x, w) {
  n <- length(x)
  sw <- sum(w)
  
  # For each observation, compute weighted ECDF: F(x_i) = sum(w[x <= x_i]) / sum(w)
  # Use midpoint version for ties: F_mid(x_i) = (F(x_i) + F(x_i^-)) / 2
  vals <- sort(unique(x))
  cum_w <- numeric(length(vals))
  for (k in seq_along(vals)) {
    cum_w[k] <- sum(w[x <= vals[k]])
  }
  cum_w <- cum_w / sw  # F(vals[k])
  
  # F(vals[k]^-) = F(vals[k]) - P(X = vals[k])
  point_w <- numeric(length(vals))
  for (k in seq_along(vals)) {
    point_w[k] <- sum(w[x == vals[k]]) / sw
  }
  
  # Map each observation to its midpoint ECDF value
  u <- numeric(n)
  for (i in seq_len(n)) {
    k <- match(x[i], vals)
    u[i] <- cum_w[k] - point_w[k] / 2  # midpoint
  }
  
  # Winsorize to avoid qnorm(0) = -Inf and qnorm(1) = Inf
  eps <- 1 / (2 * sw)
  u <- pmax(eps, pmin(1 - eps, u))
  
  qnorm(u)
}

# ------------------------------------------------------------------------------
# 3. Weighted Pearson correlation matrix
# ------------------------------------------------------------------------------
weighted_pearson_matrix <- function(Z, w) {
  p <- ncol(Z)
  sw <- sum(w)
  
  # Weighted means
  mu <- colSums(w * Z) / sw
  
  # Center
  Zc <- sweep(Z, 2, mu)
  
  # Weighted covariance: (1/sw) * t(Zc) %*% diag(w) %*% Zc
  wZc <- Zc * sqrt(w)  # scale rows by sqrt(w)
  S <- crossprod(wZc) / sw
  
  # Convert to correlation
  d <- sqrt(diag(S))
  d[d < 1e-15] <- 1e-15
  R <- S / (d %o% d)
  
  # Ensure exact 1s on diagonal
  diag(R) <- 1
  R
}

# ------------------------------------------------------------------------------
# 4. Partial correlation from a correlation matrix
#    cor(X, Y | Z) via matrix inversion
# ------------------------------------------------------------------------------
partial_cor <- function(R, i, j, k = integer(0)) {
  # R: correlation matrix
  # i, j: indices of the two variables
  # k: integer vector of conditioning set indices
  
  if (length(k) == 0) {
    return(R[i, j])
  }
  
  idx <- c(i, j, k)
  S <- R[idx, idx]
  
  P <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
  
  # partial cor = -P[1,2] / sqrt(P[1,1] * P[2,2])
  denom <- sqrt(abs(P[1, 1] * P[2, 2]))
  if (denom < 1e-15) return(0)
  
  -P[1, 2] / denom
}

# ------------------------------------------------------------------------------
# 5. CI test using partial correlation + Fisher's z
#    Effective sample size accounts for survey weights
# ------------------------------------------------------------------------------
ci_test_weighted <- function(R, i, j, k, n_eff, alpha) {
  # Returns TRUE if X_i _||_ X_j | X_k (i.e., independent)
  
  rho <- partial_cor(R, i, j, k)
  q <- length(k)
  
  # Fisher's z transform
  z <- 0.5 * log((1 + rho) / (1 - rho + 1e-15))
  
  # Standard error uses effective sample size
  se <- 1 / sqrt(n_eff - q - 3)
  
  # Two-sided test
  stat <- abs(z) / se
  pval <- 2 * pnorm(stat, lower.tail = FALSE)
  
  list(independent = (pval > alpha), pval = pval, stat = stat)
}

# ------------------------------------------------------------------------------
# 6. Effective sample size (Kish's approximation)
# ------------------------------------------------------------------------------
kish_eff_n <- function(weights) {
  (sum(weights))^2 / sum(weights^2)
}

# ------------------------------------------------------------------------------
# 7. PC-Stable Algorithm (main function)
# ------------------------------------------------------------------------------
pc_stable_weighted <- function(data, weights, alpha = 0.05, max_ord = Inf,
                               verbose = FALSE) {
  # data:    n x p matrix of ordinal integer-coded columns
  # weights: length-n survey weight vector
  # alpha:   significance level for CI tests
  # max_ord: maximum conditioning set size
  # verbose: print progress
  #
  
  # Returns: list with
  #   $amat   - p x p adjacency matrix (CPDAG)
  #   $sepset - separation sets
  
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  varnames <- colnames(data)
  if (is.null(varnames)) varnames <- paste0("V", seq_len(p))
  
  stopifnot(length(weights) == n)
  weights <- weights / sum(weights) * n  # normalize so sum = n
  
  n_eff <- kish_eff_n(weights)
  if (verbose) cat("Effective sample size (Kish):", round(n_eff, 1), "\n")
  
  if (n_eff < p + 4) warning("Effective sample size is very small relative to p.")
  
  # Step 0: Compute weighted Gaussian copula correlation matrix
  if (verbose) cat("Computing weighted Gaussian copula correlation matrix...\n")
  R <- weighted_copula_cormat(data, weights)
  
  # Initialize complete undirected graph
  G <- matrix(TRUE, p, p)
  diag(G) <- FALSE
  rownames(G) <- colnames(G) <- varnames
  
  # Separation sets: sepset[[i]][[j]] stores the separating set for (i, j)
  sepset <- lapply(seq_len(p), function(i) {
    lapply(seq_len(p), function(j) NULL)
  })
  
  ord <- 0  # current conditioning set size
  
  # PC-Stable: skeleton discovery
  while (ord <= max_ord) {
    if (verbose) cat("Order:", ord, "\n")
    
    # STABLE modification: fix adjacency at start of each order
    G_stable <- G
    any_removed <- FALSE
    
    for (i in seq_len(p)) {
      # Neighbors of i in the STABLE (frozen) graph
      adj_i <- which(G_stable[i, ])
      
      for (j in adj_i) {
        if (i >= j) next  # process each pair once
        if (!G[i, j]) next  # already removed in this order
        
        # Possible conditioning sets: adj(i)\{j} from the frozen graph
        candidates <- setdiff(adj_i, j)
        
        if (length(candidates) < ord) next
        
        # Enumerate all subsets of size ord
        if (ord == 0) {
          combos <- list(integer(0))
        } else {
          combos <- combn(candidates, ord, simplify = FALSE)
        }
        
        for (S in combos) {
          test <- ci_test_weighted(R, i, j, S, n_eff, alpha)
          
          if (test$independent) {
            if (verbose) {
              cat(sprintf("  %s _||_ %s | {%s}  p=%.4f\n",
                          varnames[i], varnames[j],
                          paste(varnames[S], collapse = ", "),
                          test$pval))
            }
            G[i, j] <- G[j, i] <- FALSE
            sepset[[i]][[j]] <- sepset[[j]][[i]] <- S
            any_removed <- TRUE
            break
          }
        }
      }
    }
    
    # Check if any edge has enough neighbors for the next order
    max_adj <- 0
    for (i in seq_len(p)) {
      max_adj <- max(max_adj, sum(G[i, ]) - 1)
    }
    
    if (max_adj < ord + 1) break
    ord <- ord + 1
  }
  
  if (verbose) cat("Skeleton complete. Orienting edges...\n")
  
  # Step 1: Orient v-structures
  # Convert to a PDAG: 0 = no edge, 1 = tail, 2 = arrowhead
  pdag <- matrix(0L, p, p)
  rownames(pdag) <- colnames(pdag) <- varnames
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (G[i, j]) pdag[i, j] <- 1L  # undirected: i -- j stored as 1 in both
    }
  }
  
  # Find unshielded triples i - k - j where i and j are not adjacent
  for (k in seq_len(p)) {
    adj_k <- which(G[k, ])
    if (length(adj_k) < 2) next
    
    pairs <- combn(adj_k, 2)
    for (col in seq_len(ncol(pairs))) {
      i <- pairs[1, col]
      j <- pairs[2, col]
      
      # Unshielded: i and j not adjacent
      if (G[i, j]) next
      
      # V-structure: i -> k <- j iff k not in sepset(i, j)
      S <- sepset[[i]][[j]]
      if (is.null(S)) S <- sepset[[j]][[i]]
      
      if (!(k %in% S)) {
        # Orient i -> k <- j
        pdag[i, k] <- 2L  # arrowhead at k
        pdag[k, i] <- 3L  # tail at i
        pdag[j, k] <- 2L
        pdag[k, j] <- 3L
      }
    }
  }
  
  # Step 2: Apply Meek's orientation rules (R1-R3) until no changes
  pdag <- apply_meek_rules(pdag, p, verbose)
  
  # Convert pdag to a simpler adjacency matrix for output
  # amat[i,j] = 1 means i -> j; amat[i,j] = amat[j,i] = 1 means i -- j
  amat <- matrix(0L, p, p)
  rownames(amat) <- colnames(amat) <- varnames
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (pdag[i, j] == 2L) {
        # i -> j (arrowhead at j, tail at i)
        amat[i, j] <- 1L
      } else if (pdag[i, j] == 1L) {
        # undirected
        amat[i, j] <- 1L
      }
    }
  }
  
  list(
    amat = amat,
    sepset = sepset,
    cormat = R,
    n_eff = n_eff,
    alpha = alpha,
    varnames = varnames
  )
}

# ------------------------------------------------------------------------------
# 8. Meek's orientation rules
# ------------------------------------------------------------------------------
apply_meek_rules <- function(pdag, p, verbose = FALSE) {
  # pdag coding: 0 = no edge
  #   undirected i--j:  pdag[i,j]=1, pdag[j,i]=1
  #   directed i->j:    pdag[i,j]=2 (arrowhead), pdag[j,i]=3 (tail)
  
  is_directed <- function(a, b) pdag[a, b] == 2L && pdag[b, a] == 3L
  is_undirected <- function(a, b) pdag[a, b] == 1L && pdag[b, a] == 1L
  has_edge <- function(a, b) pdag[a, b] > 0L
  
  orient <- function(a, b) {
    pdag[a, b] <<- 2L
    pdag[b, a] <<- 3L
  }
  
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (i == j) next
        if (!is_undirected(i, j)) next
        
        # R1: i -> k -- j, and i and j not adjacent => orient k -> j
        # Here we check: exists k such that k -> i undirected i--j, no edge k--j
        # Actually R1: a -> b -- c, no edge a-c => orient b -> c
        # So for undirected i--j, find k: k -> i and no edge k--j
        # Wait, let me re-read: R1 is a -> b - c with no a-c => b -> c
        # For the undirected edge i--j: find k such that k->i and no edge k-j
        for (k in seq_len(p)) {
          if (k == i || k == j) next
          if (is_directed(k, i) && !has_edge(k, j)) {
            orient(i, j)
            changed <- TRUE
            if (verbose) cat(sprintf("  Meek R1: %d -> %d -> %d\n", k, i, j))
            break
          }
        }
        if (is_directed(i, j)) next
        
        # R2: i -> k -> j and i -- j => orient i -> j
        for (k in seq_len(p)) {
          if (k == i || k == j) next
          if (is_directed(i, k) && is_directed(k, j)) {
            orient(i, j)
            changed <- TRUE
            if (verbose) cat(sprintf("  Meek R2: %d -> %d -> %d\n", i, k, j))
            break
          }
        }
        if (is_directed(i, j)) next
        
        # R3: i -- k1 -> j, i -- k2 -> j, k1 and k2 not adjacent => i -> j
        adj_i_undir <- which(vapply(seq_len(p), function(x) is_undirected(i, x), logical(1)))
        adj_i_undir <- setdiff(adj_i_undir, j)
        if (length(adj_i_undir) > 0) {
          dir_to_j <- adj_i_undir[vapply(adj_i_undir, function(x) is_directed(x, j), logical(1))]
        } else {
          dir_to_j <- integer(0)
        }
        if (length(dir_to_j) >= 2) {
          pairs <- combn(dir_to_j, 2)
          for (col in seq_len(ncol(pairs))) {
            if (!has_edge(pairs[1, col], pairs[2, col])) {
              orient(i, j)
              changed <- TRUE
              if (verbose) cat(sprintf("  Meek R3: orient %d -> %d\n", i, j))
              break
            }
          }
        }
      }
    }
  }
  
  pdag
}

# ------------------------------------------------------------------------------
# 9. Utility: print and plot
# ------------------------------------------------------------------------------
print_cpdag <- function(result) {
  cat("CPDAG adjacency matrix:\n")
  amat <- result$amat
  p <- nrow(amat)
  vn <- result$varnames
  
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i == j) next
      if (amat[i, j] == 1 && amat[j, i] == 1 && i < j) {
        cat(sprintf("  %s -- %s\n", vn[i], vn[j]))
      } else if (amat[i, j] == 1 && amat[j, i] == 0) {
        cat(sprintf("  %s -> %s\n", vn[i], vn[j]))
      }
    }
  }
  cat(sprintf("\nEffective n (Kish): %.1f\n", result$n_eff))
}

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================
if (FALSE) {
  set.seed(42)
  n <- 500
  
  # Simulate ordinal data (5 variables, 5-point Likert)
  Z1 <- rnorm(n)
  Z2 <- 0.6 * Z1 + rnorm(n, sd = 0.8)
  Z3 <- 0.5 * Z2 + rnorm(n, sd = 0.87)
  Z4 <- 0.4 * Z1 + 0.3 * Z3 + rnorm(n, sd = 0.75)
  Z5 <- rnorm(n)
  
  ordinalize <- function(z, k = 5) {
    as.integer(cut(z, breaks = quantile(z, probs = seq(0, 1, length.out = k + 1)),
                   include.lowest = TRUE))
  }
  
  data <- cbind(
    X1 = ordinalize(Z1),
    X2 = ordinalize(Z2),
    X3 = ordinalize(Z3),
    X4 = ordinalize(Z4),
    X5 = ordinalize(Z5)
  )
  
  # Survey weights (e.g., inverse probability weights)
  weights <- runif(n, 0.5, 2.0)
  
  # Run
  result <- pc_stable_weighted(data, weights, alpha = 0.05, verbose = TRUE)
  print_cpdag(result)
  
  # Access the adjacency matrix
  result$amat
  
  # Access the weighted polychoric correlation matrix
  result$cormat
}

if(FALSE){
  # ============================================================================
  # Example 2: Larger case (10 variables, 15 true edges, n=6000)
  #
  # True DAG (10 nodes, 15 edges):
  #   X1  -> X2, X3, X4
  #   X2  -> X5, X6
  #   X3  -> X5, X7
  #   X4  -> X7, X8
  #   X5  -> X9
  #   X6  -> X9, X10
  #   X7  -> X8, X10
  #   X8  -> X10
  #   X9  -> X10
  #
  # Edges (15 total):
  #   1->2, 1->3, 1->4, 2->5, 2->6,
  #   3->5, 3->7, 4->7, 4->8, 5->9,
  #   6->9, 6->10, 7->8, 7->10, 8->10
  # ============================================================================
  set.seed(123)
  n <- 6000
  
  # Coefficients chosen to keep signals moderate
  X1  <- rnorm(n)
  X2  <- 0.50 * X1 + rnorm(n, sd = 0.87)
  X3  <- 0.45 * X1 + rnorm(n, sd = 0.89)
  X4  <- 0.40 * X1 + rnorm(n, sd = 0.92)
  X5  <- 0.35 * X2 + 0.35 * X3 + rnorm(n, sd = 0.87)
  X6  <- 0.50 * X2 + rnorm(n, sd = 0.87)
  X7  <- 0.40 * X3 + 0.35 * X4 + rnorm(n, sd = 0.82)
  X8  <- 0.30 * X4 + 0.35 * X7 + rnorm(n, sd = 0.85)
  X9  <- 0.40 * X5 + 0.35 * X6 + rnorm(n, sd = 0.83)
  X10 <- 0.25 * X6 + 0.30 * X7 + 0.25 * X8 + rnorm(n, sd = 0.80)
  
  ordinalize <- function(z, k = 5) {
    as.integer(cut(z, breaks = quantile(z, probs = seq(0, 1, length.out = k + 1)),
                   include.lowest = TRUE))
  }
  
  data2 <- cbind(
    X1 = ordinalize(X1), X2 = ordinalize(X2),  X3 = ordinalize(X3),
    X4 = ordinalize(X4), X5 = ordinalize(X5),  X6 = ordinalize(X6),
    X7 = ordinalize(X7), X8 = ordinalize(X8),  X9 = ordinalize(X9),
    X10 = ordinalize(X10)
  )
  
  # Survey weights: simulate stratified design with unequal probabilities
  strata <- sample(1:4, n, replace = TRUE, prob = c(0.1, 0.2, 0.3, 0.4))
  base_w <- c(4, 2, 1.33, 1)[strata]  # inverse of selection prob
  weights2 <- base_w * runif(n, 0.8, 1.2)  # add some noise
  
  cat("=== Example 2: 10 variables, 15 true edges, n=6000 ===\n\n")
  
  result2 <- pc_stable_weighted(data2, weights2, alpha = 0.01, verbose = TRUE)
  print_cpdag(result2)
  
  # Compare to true edges
  true_edges <- rbind(
    c(1,2), c(1,3), c(1,4), c(2,5), c(2,6),
    c(3,5), c(3,7), c(4,7), c(4,8), c(5,9),
    c(6,9), c(6,10), c(7,8), c(7,10), c(8,10)
  )
  
  cat("\n--- Recovery comparison ---\n")
  # Check skeleton (ignoring direction)
  skeleton <- (result2$amat | t(result2$amat)) * 1L
  n_true <- nrow(true_edges)
  recovered <- 0
  for (r in seq_len(n_true)) {
    i <- true_edges[r, 1]; j <- true_edges[r, 2]
    status <- if (skeleton[i, j] == 1) "FOUND" else "MISSED"
    if (skeleton[i, j] == 1) recovered <- recovered + 1
    cat(sprintf("  %s -> %s : %s\n", result2$varnames[i], result2$varnames[j], status))
  }
  
  # Count false positives
  fp <- 0
  p <- ncol(data2)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      if (skeleton[i, j] == 1) {
        is_true <- any(apply(true_edges, 1, function(e) {
          (e[1] == i & e[2] == j) | (e[1] == j & e[2] == i)
        }))
        if (!is_true) {
          fp <- fp + 1
          cat(sprintf("  FALSE POSITIVE: %s -- %s\n", result2$varnames[i], result2$varnames[j]))
        }
      }
    }
  }
  
  cat(sprintf("\nSkeleton: %d/%d true edges recovered, %d false positives\n",
              recovered, n_true, fp))
}
