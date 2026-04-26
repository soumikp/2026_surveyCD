# ==============================================================================
# swa-oBN: Survey-Weighted, Covariate-Adjusted Ordinal Bayesian Network
# with Expert Constraints and Two-Step Scoring
#
# Extensions over Ni & Mallick (2022) OCD:
#   (a) Covariate adjustment (Z) in ordinal regression
#   (b) Survey weights with effective sample size in BIC
#   (c) Expert constraints via blacklist / whitelist
#   (d) Two-step scoring: penalized pseudo-BIC for skeleton discovery
#       (add/delete), unpenalized pseudo-likelihood ratio for edge
#       orientation (reversal). This separation ensures DAG selection
#       consistency for heterogeneous ordinal levels (Theorem 1).
#
# Weight convention:
#   User-supplied weights w_i = 1/pi_i are normalized internally to
#   sum to n (the sample size). The pseudo-log-likelihood computed by
#   polr/glm under these normalized weights targets the sample-scaled
#   quantity (n/N) * sum_{i in U_N} log p(x_i | G, theta, z_i).
#   The design-effect correction (n_eff / n) rescales this to the
#   effective-information scale before applying the log(n_eff) penalty,
#   following Lumley & Scott (2015).
#
# Note on polr vs svyolr:
#   We use MASS::polr with normalized weights rather than survey::svyolr.
#   For pseudo-log-likelihood computation (which is all we need for BIC
#   scoring), the two are equivalent: both maximize sum_i w_i log p(y_i|.).
#   The difference is in standard error computation — svyolr uses the
#   sandwich estimator for design-based inference, while polr does not.
#   Since we only extract logLik (not standard errors) from the fitted
#   model, polr is sufficient and substantially faster. The same applies
#   to stats::glm vs survey::svyglm for binary nodes.
# ==============================================================================

library(MASS)
library(stats)
library(igraph)
library(gRbase)

# ==============================================================================
# SECTION 1: GRAPH UTILITIES
# ==============================================================================

#' Check if adding edge j -> i preserves DAG property
#' @noRd
is_dag_after_add <- function(i, j, gam) {
  if (gam[i, j]) return(TRUE)
  gam[i, j] <- 1
  return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam)))
}

#' Check if reversing the edge between i and j preserves DAG property
#' @noRd
is_dag_after_rev <- function(i, j, gam) {
  if (gam[i, j] == gam[j, i]) return(FALSE)
  tmp <- gam[i, j]
  gam[j, i] <- tmp
  gam[i, j] <- !tmp
  return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam)))
}

#' Check if an edge operation is allowed by expert constraints
#' @param from Source node index
#' @param to Target node index
#' @param operation One of "add", "del", "rev"
#' @param blacklist Matrix (k x 2) of forbidden directed edges
#' @param whitelist Matrix (k x 2) of required directed edges
#' @noRd
is_allowed <- function(from, to, operation, blacklist = NULL, whitelist = NULL) {
  
  if (operation == "add") {
    if (!is.null(blacklist)) {
      if (any(blacklist[, 1] == from & blacklist[, 2] == to)) return(FALSE)
    }
    return(TRUE)
  }
  
  if (operation == "del") {
    if (!is.null(whitelist)) {
      if (any(whitelist[, 1] == from & whitelist[, 2] == to)) return(FALSE)
    }
    return(TRUE)
  }
  
  if (operation == "rev") {
    if (!is.null(whitelist)) {
      if (any(whitelist[, 1] == from & whitelist[, 2] == to)) return(FALSE)
    }
    if (!is.null(blacklist)) {
      if (any(blacklist[, 1] == to & blacklist[, 2] == from)) return(FALSE)
    }
    return(TRUE)
  }
  
  return(TRUE)
}

# ==============================================================================
# SECTION 2: EFFECTIVE SAMPLE SIZE
# ==============================================================================

#' Compute effective sample size (Kish, 1965)
#' n_eff = (sum(w))^2 / sum(w^2)
#' Falls back to n when weights are uniform
#' @noRd
compute_n_eff <- function(weights) {
  sw <- sum(weights)
  sw2 <- sum(weights^2)
  if (sw2 == 0) return(length(weights))
  return(sw^2 / sw2)
}

# ==============================================================================
# SECTION 3: SCORE FUNCTION (WEIGHTED ORDINAL REGRESSION + COVARIATES)
# ==============================================================================

#' Fit a single node's ordinal regression and return its score
#'
#' Model: P(Y <= l | X_pa, Z) = F(gamma_l - beta' X_pa - delta' Z - alpha)
#'
#' The function returns either the full pseudo-BIC (penalized) or just the
#' rescaled pseudo-log-likelihood (unpenalized), controlled by the `penalize`
#' argument. This separation is the core of the two-step scoring approach:
#'
#'   - penalize = TRUE  (default): Used for skeleton discovery (add/delete).
#'     Returns  -2 * ll * (n_eff/n) + k * log(n_eff)
#'
#'   - penalize = FALSE: Used for edge orientation (reversal within MEC).
#'     Returns  -2 * ll * (n_eff/n)
#'     The penalty is omitted because within-MEC discrimination is driven
#'     entirely by the ordinal likelihood (Theorem 1), and including the
#'     penalty introduces finite-sample bias when ordinal levels are
#'     heterogeneous across variables.
#'
#' @param y_vec Factor response for this node
#' @param x_df  Data frame of parent node columns (NULL if no parents)
#' @param z_df  Data frame of covariate columns (NULL if no covariates)
#' @param weights Numeric vector of observation weights (normalized to sum to n)
#' @param ic    "bic" or "aic" (only affects penalized scores)
#' @param method Link function: "probit" or "logistic"
#' @param n_eff Effective sample size (used in BIC penalty and deff correction)
#' @param nq_y  Number of levels in y_vec
#' @param penalize Logical. If TRUE, return full pseudo-BIC. If FALSE, return
#'        unpenalized pseudo-log-likelihood (rescaled). Default TRUE.
#' @noRd
node_score <- function(y_vec, x_df, z_df, weights, ic, method, n_eff, nq_y,
                       penalize = TRUE) {
  
  # --- Build data frame and formula ---
  data_list <- list(Y = y_vec)
  pred_names <- character(0)
  
  if (!is.null(x_df) && ncol(x_df) > 0) {
    data_list <- c(data_list, as.list(x_df))
    pred_names <- c(pred_names, colnames(x_df))
  }
  if (!is.null(z_df) && ncol(z_df) > 0) {
    data_list <- c(data_list, as.list(z_df))
    pred_names <- c(pred_names, colnames(z_df))
  }
  
  df <- as.data.frame(data_list)
  df$.w <- weights
  
  if (length(pred_names) > 0) {
    form <- stats::as.formula(paste("Y ~", paste(pred_names, collapse = " + ")))
  } else {
    form <- stats::as.formula("Y ~ 1")
  }
  
  # --- Score computation ---
  # The pseudo-log-likelihood from polr/glm under weights normalized to sum
  # to n is on the "sample scale". The deff correction (n_eff / n) rescales
  # it to the effective-information scale so that likelihood and penalty are
  # commensurate (Lumley & Scott, 2015).
  calc_score <- function(model) {
    ll <- as.numeric(stats::logLik(model))
    n_raw <- nrow(df)
    deff_correction <- n_eff / n_raw  # = n_eff / n, typically < 1
    
    # Rescaled negative pseudo-log-likelihood (always computed)
    neg2ll_rescaled <- -2 * ll * deff_correction
    
    if (!penalize) {
      return(neg2ll_rescaled)
    }
    
    # Full pseudo-BIC: add complexity penalty
    k <- length(stats::coef(model))
    if (inherits(model, "polr")) k <- k + length(model$zeta)
    
    if (ic == "bic") {
      return(neg2ll_rescaled + k * log(n_eff))
    } else {
      # AIC: no deff correction on likelihood for AIC (standard convention)
      return(-2 * ll + 2 * k)
    }
  }
  
  # --- Fit model ---
  score <- Inf
  fitted <- FALSE
  
  if (nq_y > 2) {
    # Ordinal regression via MASS::polr
    tryCatch({
      mod <- MASS::polr(form, data = df, weights = .w,
                        method = method, Hess = FALSE)
      score <- calc_score(mod)
      fitted <- TRUE
    }, error = function(e) {})
    
    # Retry with random start if first attempt failed
    if (!fitted) {
      tryCatch({
        X_mat <- stats::model.matrix(form, df)
        n_coef <- ncol(X_mat) - 1
        n_zeta <- nq_y - 1
        mod <- MASS::polr(form, data = df, weights = .w,
                          method = method, Hess = FALSE,
                          start = sort(stats::rnorm(n_coef + n_zeta)))
        score <- calc_score(mod)
        fitted <- TRUE
      }, error = function(e) {})
    }
  } else {
    # Binary case: use glm with probit or logit link
    tryCatch({
      mod <- stats::glm(form, data = df, weights = .w,
                        family = stats::binomial(link = method))
      score <- calc_score(mod)
      fitted <- TRUE
    }, error = function(e) {})
    
    if (!fitted) {
      tryCatch({
        X_mat <- stats::model.matrix(form, df)
        n_coef <- ncol(X_mat)
        mod <- stats::glm(form, data = df, weights = .w,
                          family = stats::binomial(link = method),
                          start = stats::rnorm(n_coef))
        score <- calc_score(mod)
        fitted <- TRUE
      }, error = function(e) {})
    }
  }
  
  if (is.na(score)) score <- Inf
  return(score)
}

# ==============================================================================
# SECTION 4: CACHED SCORE WRAPPER
# ==============================================================================

#' Wrapper around node_score with hash-map caching
#'
#' Cache key: "node_{i}_parents_{sorted parent indices}_{pen|unpen}"
#' This avoids redundant ordinal regressions across search iterations.
#' The penalize flag is part of the cache key because the same parent-set
#' configuration yields different scores under penalized vs unpenalized scoring.
#' @noRd
cached_node_score <- function(i, parent_vec, y, Z, weights, ic, method,
                              n_eff, nq, cache = NULL, penalize = TRUE) {
  
  parent_vec <- (parent_vec > 0)
  
  # --- Cache lookup ---
  parent_str <- paste(which(parent_vec), collapse = "_")
  pen_str <- if (penalize) "pen" else "unpen"
  cache_key <- paste0("node_", i, "_pa_", parent_str, "_", pen_str)
  if (!is.null(cache) && exists(cache_key, envir = cache)) {
    return(get(cache_key, envir = cache))
  }
  
  # --- Build inputs ---
  y_vec <- y[, i]
  nq_y <- nq[i]
  
  x_df <- NULL
  if (any(parent_vec)) {
    x_df <- y[, parent_vec, drop = FALSE]
    colnames(x_df) <- paste0("X", which(parent_vec))
  }
  
  z_df <- NULL
  if (!is.null(Z)) {
    z_df <- as.data.frame(Z)
    if (is.null(colnames(z_df))) colnames(z_df) <- paste0("Z", seq_len(ncol(z_df)))
  }
  
  score <- node_score(y_vec, x_df, z_df, weights, ic, method, n_eff, nq_y,
                      penalize = penalize)
  
  # --- Cache store ---
  if (!is.null(cache)) assign(cache_key, score, envir = cache)
  
  return(score)
}

# ==============================================================================
# SECTION 5: GREEDY HILL-CLIMBING SEARCH (TWO-STEP SCORING)
# ==============================================================================

#' Greedy hill-climbing DAG search with two-step scoring
#'
#' The scoring strategy implements the theoretical separation between
#' skeleton discovery and edge orientation:
#'
#'   ADD / DELETE moves change the skeleton (adjacency structure).
#'     -> Scored by full pseudo-BIC (penalized = TRUE).
#'     -> Penalty prevents overfitting (Theorem 2b).
#'
#'   REVERSAL moves change edge direction within the current skeleton.
#'     -> Scored by unpenalized pseudo-likelihood ratio (penalized = FALSE).
#'     -> Ordinal identifiability drives orientation (Theorem 1).
#'     -> No penalty needed; avoids finite-sample bias from heterogeneous levels.
#'
#' @param y Data frame of ordinal factor columns (the causal variables)
#' @param Z Data frame of covariates (adjusted for, not in the DAG)
#' @param weights Numeric vector of survey weights (normalized to sum to n)
#' @param gam Initial adjacency matrix (q x q logical)
#' @param ic "bic" or "aic"
#' @param method Link function
#' @param blacklist Matrix (k x 2) of forbidden edges (col1 -> col2)
#' @param whitelist Matrix (k x 2) of required edges (col1 -> col2)
#' @param verbose Print progress
#' @param maxit Maximum iterations
#' @noRd
greedy_search <- function(y, Z = NULL, weights = NULL, gam = NULL,
                          ic = "bic", method = "probit",
                          blacklist = NULL, whitelist = NULL,
                          verbose = FALSE, maxit = 100) {
  n <- nrow(y)
  q <- ncol(y)
  nq <- sapply(seq_len(q), function(i) nlevels(y[, i]))
  n_eff <- compute_n_eff(weights)
  
  cache <- new.env(hash = TRUE)
  
  if (is.null(gam)) {
    gam <- matrix(FALSE, q, q)
  } else {
    gam <- (gam != 0)
  }
  
  # Enforce whitelist edges in the initial graph
  if (!is.null(whitelist)) {
    for (r in seq_len(nrow(whitelist))) {
      from <- whitelist[r, 1]
      to <- whitelist[r, 2]
      gam[to, from] <- TRUE  # gam[i, j] = TRUE means j -> i
    }
  }
  
  # Build index of candidate neighbors for each node
  ind_q <- matrix(0L, q, q - 1)
  for (i in seq_len(q)) {
    ind_q[i, ] <- setdiff(seq_len(q), i)
  }
  
  # Score initial graph with PENALIZED scores (for skeleton tracking)
  ic_best_pen <- numeric(q)
  for (i in seq_len(q)) {
    ic_best_pen[i] <- cached_node_score(i, gam[i, ], y, Z, weights, ic, method,
                                        n_eff, nq, cache, penalize = TRUE)
  }
  
  # Also maintain UNPENALIZED scores (for reversal comparisons)
  ic_best_unpen <- numeric(q)
  for (i in seq_len(q)) {
    ic_best_unpen[i] <- cached_node_score(i, gam[i, ], y, Z, weights, ic, method,
                                          n_eff, nq, cache, penalize = FALSE)
  }
  
  iter <- 0
  improving <- TRUE
  
  while (improving && iter < maxit) {
    iter <- iter + 1
    improving <- FALSE
    best_gain <- 0
    best_action <- NULL
    
    gam_tmp <- gam
    
    # --- Scan ADD and DELETE moves (scored by penalized pseudo-BIC) ---
    for (i in seq_len(q)) {
      for (jj in seq_len(q - 1)) {
        j <- ind_q[i, jj]  # j is the potential parent, i is the child
        is_present <- gam[i, j]
        
        if (is_present) {
          if (!is_allowed(from = j, to = i, "del", blacklist, whitelist)) next
        } else {
          if (!is_allowed(from = j, to = i, "add", blacklist, whitelist)) next
          if (!is_dag_after_add(i, j, gam)) next
        }
        
        gam_tmp[i, j] <- !is_present
        ic_new <- cached_node_score(i, gam_tmp[i, ], y, Z, weights, ic,
                                    method, n_eff, nq, cache, penalize = TRUE)
        gain <- ic_best_pen[i] - ic_new
        if (is.na(gain)) gain <- -Inf
        
        if (gain > best_gain) {
          best_gain <- gain
          best_action <- list(
            type = ifelse(is_present, "del", "add"),
            node_i = i, node_j = j,
            ic_new_pen_i = ic_new
          )
        }
        gam_tmp[i, j] <- is_present
      }
    }
    
    # --- Scan REVERSAL moves (scored by unpenalized pseudo-likelihood) ---
    for (i in seq_len(q)) {
      for (jj in seq_len(q - 1)) {
        j <- ind_q[i, jj]
        if (!is_dag_after_rev(i, j, gam)) next
        
        # Determine which direction currently exists
        if (gam[i, j] && !gam[j, i]) {
          from_cur <- j; to_cur <- i
        } else if (gam[j, i] && !gam[i, j]) {
          from_cur <- i; to_cur <- j
        } else {
          next
        }
        
        if (!is_allowed(from_cur, to_cur, "rev", blacklist, whitelist)) next
        
        # Apply reversal in temporary matrix
        tmp_val <- gam_tmp[i, j]
        gam_tmp[j, i] <- tmp_val
        gam_tmp[i, j] <- !tmp_val
        
        # Score reversals with UNPENALIZED pseudo-likelihood
        ic_new_unpen_i <- cached_node_score(i, gam_tmp[i, ], y, Z, weights, ic,
                                            method, n_eff, nq, cache,
                                            penalize = FALSE)
        ic_new_unpen_j <- cached_node_score(j, gam_tmp[j, ], y, Z, weights, ic,
                                            method, n_eff, nq, cache,
                                            penalize = FALSE)
        
        gain <- (ic_best_unpen[i] - ic_new_unpen_i) +
          (ic_best_unpen[j] - ic_new_unpen_j)
        if (is.na(gain)) gain <- -Inf
        
        if (gain > best_gain) {
          best_gain <- gain
          best_action <- list(
            type = "rev",
            node_i = i, node_j = j,
            ic_new_unpen_i = ic_new_unpen_i,
            ic_new_unpen_j = ic_new_unpen_j
          )
        }
        
        # Undo reversal
        gam_tmp[i, j] <- tmp_val
        gam_tmp[j, i] <- !tmp_val
      }
    }
    
    # --- Apply best move ---
    if (best_gain > 0 && !is.null(best_action)) {
      improving <- TRUE
      a <- best_action
      
      if (a$type == "add") {
        gam[a$node_i, a$node_j] <- TRUE
        ic_best_pen[a$node_i] <- a$ic_new_pen_i
        # Update unpenalized score for the affected node
        ic_best_unpen[a$node_i] <- cached_node_score(
          a$node_i, gam[a$node_i, ], y, Z, weights, ic, method,
          n_eff, nq, cache, penalize = FALSE)
        
      } else if (a$type == "del") {
        gam[a$node_i, a$node_j] <- FALSE
        ic_best_pen[a$node_i] <- a$ic_new_pen_i
        ic_best_unpen[a$node_i] <- cached_node_score(
          a$node_i, gam[a$node_i, ], y, Z, weights, ic, method,
          n_eff, nq, cache, penalize = FALSE)
        
      } else if (a$type == "rev") {
        tmp <- gam[a$node_i, a$node_j]
        gam[a$node_j, a$node_i] <- tmp
        gam[a$node_i, a$node_j] <- !tmp
        # Update unpenalized scores (these were computed during the scan)
        ic_best_unpen[a$node_i] <- a$ic_new_unpen_i
        ic_best_unpen[a$node_j] <- a$ic_new_unpen_j
        # Update penalized scores (needed for subsequent add/delete decisions)
        ic_best_pen[a$node_i] <- cached_node_score(
          a$node_i, gam[a$node_i, ], y, Z, weights, ic, method,
          n_eff, nq, cache, penalize = TRUE)
        ic_best_pen[a$node_j] <- cached_node_score(
          a$node_j, gam[a$node_j, ], y, Z, weights, ic, method,
          n_eff, nq, cache, penalize = TRUE)
      }
      
      if (verbose) {
        cat(sprintf("Iter %d | %s | BIC = %.2f | Cache: %d entries\n",
                    iter, a$type, sum(ic_best_pen), length(ls(cache))))
      }
    }
  }
  
  if (iter == maxit) {
    warning("Maximum iterations reached. Algorithm may not have converged.")
  }
  
  return(list(gam = gam + 0, ic_best = sum(ic_best_pen)))
}

# ==============================================================================
# SECTION 6: EXHAUSTIVE SEARCH (q = 2 or 3 only)
# ==============================================================================

#' Exhaustive search over all DAGs for q = 2 or 3
#'
#' For bivariate exhaustive search (q = 2), this function uses unpenalized
#' pseudo-likelihood to compare the three possible DAGs (null, X1->X2, X2->X1),
#' since all non-null DAGs have the same skeleton and differ only in orientation.
#' For q = 3, the full pseudo-BIC is used for skeleton determination and
#' unpenalized comparison is applied within each equivalence class.
#'
#' Note: For simplicity, this implementation uses full pseudo-BIC for all
#' q = 3 DAGs. The two-step logic is primarily relevant for greedy search
#' where add/delete and reversal moves are distinct operations.
#' @noRd
exhaustive_search <- function(y, Z = NULL, weights = NULL, gam_list = NULL,
                              ic = "bic", method = "probit",
                              blacklist = NULL, whitelist = NULL) {
  n <- nrow(y)
  q <- ncol(y)
  nq <- sapply(seq_len(q), function(i) nlevels(y[, i]))
  n_eff <- compute_n_eff(weights)
  
  cache <- new.env(hash = TRUE)
  
  # --- Enumerate all DAGs for q = 2 or 3 ---
  if (is.null(gam_list)) {
    if (q == 2) {
      gam_list <- array(0, c(2, 2, 3))
      # DAG 1: empty
      # DAG 2: 2 -> 1
      gam_list[1, 2, 2] <- 1
      # DAG 3: 1 -> 2
      gam_list[2, 1, 3] <- 1
    } else if (q == 3) {
      gam_list <- array(0, c(3, 3, 25))
      # All 25 three-node DAGs (same enumeration as original OCD)
      gam_list[1, 2, 2] <- 1
      gam_list[1, 3, 3] <- 1
      gam_list[2, 3, 4] <- 1
      gam_list[2, 1, 5] <- 1
      gam_list[3, 1, 6] <- 1
      gam_list[3, 2, 7] <- 1
      gam_list[1, 2, 8] <- gam_list[1, 3, 8] <- 1
      gam_list[2, 3, 9] <- gam_list[2, 1, 9] <- 1
      gam_list[3, 1, 10] <- gam_list[3, 2, 10] <- 1
      gam_list[2, 1, 11] <- gam_list[3, 1, 11] <- 1
      gam_list[3, 2, 12] <- gam_list[1, 2, 12] <- 1
      gam_list[1, 3, 13] <- gam_list[2, 3, 13] <- 1
      gam_list[1, 2, 14] <- gam_list[2, 3, 14] <- 1
      gam_list[1, 3, 15] <- gam_list[3, 2, 15] <- 1
      gam_list[2, 3, 16] <- gam_list[3, 1, 16] <- 1
      gam_list[2, 1, 17] <- gam_list[1, 3, 17] <- 1
      gam_list[3, 2, 18] <- gam_list[2, 1, 18] <- 1
      gam_list[3, 1, 19] <- gam_list[1, 2, 19] <- 1
      gam_list[1, 2, 20] <- gam_list[1, 3, 20] <- gam_list[2, 3, 20] <- 1
      gam_list[1, 2, 21] <- gam_list[1, 3, 21] <- gam_list[3, 2, 21] <- 1
      gam_list[2, 1, 22] <- gam_list[2, 3, 22] <- gam_list[1, 3, 22] <- 1
      gam_list[2, 1, 23] <- gam_list[2, 3, 23] <- gam_list[3, 1, 23] <- 1
      gam_list[3, 1, 24] <- gam_list[3, 2, 24] <- gam_list[1, 2, 24] <- 1
      gam_list[3, 1, 25] <- gam_list[3, 2, 25] <- gam_list[2, 1, 25] <- 1
    } else {
      stop("Exhaustive search only supports q = 2 or 3. Use greedy search.")
    }
  }
  
  n_dags <- dim(gam_list)[3]
  
  # --- Filter DAGs by constraints ---
  valid <- rep(TRUE, n_dags)
  for (m in seq_len(n_dags)) {
    gam_m <- gam_list[, , m]
    # Check blacklist: no forbidden edge should be present
    if (!is.null(blacklist) && valid[m]) {
      for (r in seq_len(nrow(blacklist))) {
        from <- blacklist[r, 1]; to <- blacklist[r, 2]
        if (gam_m[to, from] == 1) { valid[m] <- FALSE; break }
      }
    }
    if (!valid[m]) next
    # Check whitelist: all required edges must be present
    if (!is.null(whitelist)) {
      for (r in seq_len(nrow(whitelist))) {
        from <- whitelist[r, 1]; to <- whitelist[r, 2]
        if (gam_m[to, from] != 1) { valid[m] <- FALSE; break }
      }
    }
  }
  
  if (!any(valid)) stop("No DAG satisfies the given constraints.")
  
  # --- Score all valid DAGs ---
  # For q = 2, use unpenalized scoring (all non-null DAGs share a skeleton,
  # so discrimination is purely about orientation — Theorem 1 applies).
  # For q = 3, use penalized BIC (DAGs may differ in skeleton).
  use_penalized <- (q > 2)
  
  IC <- rep(Inf, n_dags)
  for (m in which(valid)) {
    IC[m] <- 0
    gam_m <- gam_list[, , m]
    for (i in seq_len(q)) {
      IC[m] <- IC[m] + cached_node_score(i, gam_m[i, ], y, Z, weights, ic,
                                         method, n_eff, nq, cache,
                                         penalize = use_penalized)
    }
  }
  
  mi <- which.min(IC)
  return(list(gam = gam_list[, , mi] + 0, ic_best = IC[mi]))
}

# ==============================================================================
# SECTION 7: MULTI-START WRAPPER
# ==============================================================================

#' Run greedy search with multiple random restarts
#'
#' @param nstart Number of random initializations (first is always empty graph)
#' @noRd
multi_start_search <- function(y, Z = NULL, weights = NULL,
                               ic = "bic", method = "probit",
                               blacklist = NULL, whitelist = NULL,
                               nstart = 5, verbose = FALSE, maxit = 100) {
  q <- ncol(y)
  results <- vector("list", nstart)
  
  for (r in seq_len(nstart)) {
    if (r == 1) {
      # Start from empty graph
      gam_init <- matrix(FALSE, q, q)
    } else {
      # Random sparse initialization
      net <- bnlearn::random.graph(
        nodes = as.character(seq_len(q)),
        method = "ordered", num = 1, prob = 1 / q
      )
      gam_init <- matrix(FALSE, q, q)
      if (nrow(net$arcs) > 0) {
        gam_init[cbind(
          as.numeric(net$arcs[, 2]),
          as.numeric(net$arcs[, 1])
        )] <- TRUE
      }
    }
    
    results[[r]] <- greedy_search(
      y, Z, weights, gam = gam_init, ic = ic, method = method,
      blacklist = blacklist, whitelist = whitelist,
      verbose = verbose, maxit = maxit
    )
  }
  
  best_idx <- which.min(sapply(results, function(x) x$ic_best))
  return(results[[best_idx]])
}

# ==============================================================================
# SECTION 8: MAIN INTERFACE
# ==============================================================================

#' swa-oBN: Main function for survey-weighted, covariate-adjusted ordinal
#' causal discovery with expert constraints and two-step scoring
#'
#' @param y Data frame where each column is an ordered factor.
#'          These are the causal variables (nodes in the DAG).
#' @param Z Optional data frame of covariates to adjust for.
#'          Not treated as DAG nodes. Can be continuous or categorical.
#' @param weights Optional numeric vector of survey/sampling weights.
#'          When NULL, uniform weights are used. Weights are normalized
#'          internally to sum to n (the sample size).
#' @param search "greedy" for hill-climbing or "exhaust" for brute force (q<=3).
#' @param ic "bic" (default, uses n_eff in penalty) or "aic".
#' @param link Link function: "probit" (default) or "logistic".
#' @param blacklist Optional matrix (k x 2) of forbidden edges.
#'          Each row (from, to) means from -> to is forbidden.
#' @param whitelist Optional matrix (k x 2) of required edges.
#'          Each row (from, to) means from -> to must be present.
#' @param nstart Number of random restarts for greedy search (default 5).
#' @param boot Number of bootstrap samples for uncertainty quantification.
#'          When NULL, returns a single point estimate.
#' @param G Optional initial DAG adjacency matrix.
#' @param verbose Print search progress.
#' @param maxit Max iterations per greedy run.
#'
#' @return If boot is NULL: list with $gam (adjacency matrix) and $ic_best.
#'         If boot > 0: list with $gam (point estimate), $ic_best,
#'         $edge_probs (bootstrap edge inclusion probabilities),
#'         and $boot_gams (list of bootstrap adjacency matrices).
#'
#' @examples
#' # Basic usage
#' fit <- swa_obn(y = my_ordinal_data)
#'
#' # With all extensions
#' fit <- swa_obn(
#'   y = my_ordinal_data,
#'   Z = data.frame(age = age_vec, sex = sex_vec),
#'   weights = survey_wts,
#'   blacklist = matrix(c(3, 1), ncol = 2),  # forbid 3 -> 1
#'   whitelist = matrix(c(1, 2), ncol = 2),  # require 1 -> 2
#'   boot = 500
#' )
#' @export
swa_obn <- function(y,
                    Z = NULL,
                    weights = NULL,
                    search = "greedy",
                    ic = "bic",
                    link = "probit",
                    blacklist = NULL,
                    whitelist = NULL,
                    nstart = 5,
                    boot = NULL,
                    G = NULL,
                    verbose = FALSE,
                    maxit = 100) {
  
  # --- Input validation ---
  y <- as.data.frame(y)
  for (j in seq_len(ncol(y))) {
    if (!is.factor(y[, j])) y[, j] <- as.ordered(y[, j])
    if (!is.ordered(y[, j])) y[, j] <- as.ordered(y[, j])
  }
  
  n <- nrow(y)
  q <- ncol(y)
  
  if (is.null(weights)) {
    weights <- rep(1, n)
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != n) stop("weights must have same length as nrow(y)")
    if (any(weights < 0)) stop("weights must be non-negative")
    # Normalize so weights sum to n (stabilizes likelihood scale).
    # The deff correction in node_score rescales to the effective-information
    # scale before applying the BIC penalty.
    weights <- weights * (n / sum(weights))
  }
  
  if (!is.null(Z)) {
    Z <- as.data.frame(Z)
    if (nrow(Z) != n) stop("Z must have same number of rows as y")
  }
  
  # --- Validate constraints ---
  if (!is.null(blacklist)) {
    blacklist <- as.matrix(blacklist)
    if (ncol(blacklist) != 2) stop("blacklist must have 2 columns (from, to)")
  }
  if (!is.null(whitelist)) {
    whitelist <- as.matrix(whitelist)
    if (ncol(whitelist) != 2) stop("whitelist must have 2 columns (from, to)")
  }
  
  # --- Single fit function ---
  fit_once <- function(y_in, Z_in, w_in) {
    if (search == "exhaust") {
      return(exhaustive_search(y_in, Z_in, w_in, gam_list = G, ic = ic,
                               method = link, blacklist = blacklist,
                               whitelist = whitelist))
    } else {
      return(multi_start_search(y_in, Z_in, w_in, ic = ic, method = link,
                                blacklist = blacklist, whitelist = whitelist,
                                nstart = nstart, verbose = verbose,
                                maxit = maxit))
    }
  }
  
  # --- Point estimate ---
  fit_point <- fit_once(y, Z, weights)
  
  if (is.null(boot)) {
    return(fit_point)
  }
  
  # --- Bootstrap for uncertainty quantification ---
  # Naive pairs bootstrap: resample (y_i, w_i, z_i) tuples with replacement.
  # Survey weights are preserved for each resampled observation. The n_eff
  # is recomputed within each bootstrap resample.
  #
  # Note: This is not a design-appropriate bootstrap (e.g., Rao-Wu rescaled).
  # For DAG selection consistency, each sufficiently large resample recovers
  # G* via Theorems 1-2 regardless of the resampling scheme, so edge inclusion
  # probabilities converge to 0 or 1 for all edges. In finite samples,
  # design-appropriate resampling may improve calibration of inclusion
  # probabilities at intermediate values.
  boot_gams <- vector("list", boot)
  adj_sum <- matrix(0, q, q)
  
  for (b in seq_len(boot)) {
    idx <- sample(n, replace = TRUE)
    y_b <- y[idx, , drop = FALSE]
    w_b <- weights[idx]
    Z_b <- if (!is.null(Z)) Z[idx, , drop = FALSE] else NULL
    
    boot_fit <- fit_once(y_b, Z_b, w_b)
    boot_gams[[b]] <- boot_fit$gam
    adj_sum <- adj_sum + boot_fit$gam
    
    if (verbose && b %% 50 == 0) {
      cat(sprintf("Bootstrap %d / %d complete\n", b, boot))
    }
  }
  
  edge_probs <- adj_sum / boot
  colnames(edge_probs) <- colnames(y)
  rownames(edge_probs) <- colnames(y)
  
  return(list(
    gam = fit_point$gam,
    ic_best = fit_point$ic_best,
    edge_probs = edge_probs,
    boot_gams = boot_gams
  ))
}