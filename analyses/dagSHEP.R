# ==============================================================================
# FINAL SHEP PIPELINE: Copula PC + Global Tiered OrdCD (ALL ORDINAL DATA)
# ==============================================================================
library(pcalg)
library(survey)
library(MASS)
library(gRbase)
library(igraph)

# ------------------------------------------------------------------------------
# 1. Copula Transform for Phase 1 Skeleton
# ------------------------------------------------------------------------------
svy_gaussian_copula_transform <- function(x, weights) {
  n <- length(x)
  ord <- order(x)
  w_cdf <- cumsum(weights[ord]) / sum(weights)
  w_cdf <- w_cdf * (n / (n + 1))
  cdf_orig <- numeric(n)
  cdf_orig[ord] <- w_cdf
  return(qnorm(cdf_orig))
}

# ------------------------------------------------------------------------------
# 2. Tier-Constrained OrdCD Functions
# ------------------------------------------------------------------------------
admissible_tiered <- function(i, j, gam, tiers) {
  # If the edge exists, we must ALWAYS be allowed to evaluate deleting it!
  if (gam[i, j]) return(TRUE) 
  
  # Forbid ADDING backward edges
  if (!is.null(tiers) && tiers[[i]] > tiers[[j]]) return(FALSE) 
  
  gam[i, j] <- 1
  return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam)))
}

admissible_rev_tiered <- function(i, j, gam, tiers) {
  # 1. ONLY evaluate a reversal if the edge i -> j actually exists!
  if (!gam[i, j]) return(FALSE) 
  
  # 2. We are reversing i -> j into j -> i. 
  # So j becomes the NEW cause. If the new cause (j) is a higher tier than 
  # the new effect (i), this is a backward arrow. Reject it!
  if (!is.null(tiers) && tiers[[j]] > tiers[[i]]) return(FALSE)
  
  # 3. Apply the reversal and check for cycles
  tmp <- gam[i, j]
  gam[j, i] <- tmp
  gam[i, j] <- !tmp
  return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam)))
}

svy_mypolr <- function(formula, data, weights, ic = "bic", method = "probit") {
  des <- svydesign(ids = ~1, weights = ~weights, data = data)
  IC <- 1e9 # Massive finite penalty prevents NaN propagation
  n_eff <- sum(weights)
  
  tryCatch({
    fit <- svyolr(formula, design = des, method = method)
    k <- length(fit$coefficients) + length(fit$zeta)
    IC <- fit$deviance + k * log(n_eff)
  }, error = function(e) {
    # If the model fails due to sparse data, it silently defaults to the 1e9 penalty
  })
  
  return(IC)
}

get_node_score_svy <- function(i, gam, data, weights, ic, method) {
  # Look at the column for node i to find its parents (gam[, i]), NOT the row
  parents <- which(gam[, i] > 0) 
  target <- colnames(data)[i]
  
  if (length(parents) > 0) {
    preds <- colnames(data)[parents]
    form <- as.formula(paste(target, "~", paste(preds, collapse = " + ")))
  } else {
    form <- as.formula(paste(target, "~ 1"))
  }
  return(svy_mypolr(form, data, weights, ic, method))
}

OrdCD_svy_greedy <- function(data, weights, tiers = NULL, gam = NULL, ic = "bic", method = "probit", maxit = 50) {
  q <- ncol(data)
  if (is.null(gam)) gam <- matrix(FALSE, q, q) else gam <- (gam != 0)
  tier_vec <- sapply(colnames(data), function(x) tiers[[x]])
  ic_best <- sapply(1:q, function(i) get_node_score_svy(i, gam, data, weights, ic, method))
  
  iter <- 0
  ic_improv <- 1
  
  while (ic_improv > 0 && iter < maxit) {
    iter <- iter + 1
    ic_improv <- -Inf
    best_move <- NULL
    
    # 1. Evaluate Add/Delete
    for (i in 1:q) {
      for (j in setdiff(1:q, i)) {
        if (admissible_tiered(i, j, gam, tier_vec)) {
          gam_temp <- gam
          gam_temp[i, j] <- !gam[i, j]
          new_score <- get_node_score_svy(i, gam_temp, data, weights, ic, method)
          improv <- ic_best[i] - new_score
          
          # SAFETY CHECK: Ignore if Inf - Inf creates NaN
          if (!is.na(improv) && improv > ic_improv) {
            ic_improv <- improv
            best_move <- list(type = if(gam[i,j]) "del" else "add", i=i, j=j, score_i=new_score)
          }
        }
      }
    }
    
    # 2. Evaluate Reversals
    for (i in 1:q) {
      for (j in setdiff(1:q, i)) {
        if (admissible_rev_tiered(i, j, gam, tier_vec)) {
          gam_temp <- gam
          gam_temp[j, i] <- gam[i, j]
          gam_temp[i, j] <- !gam[i, j]
          score_i <- get_node_score_svy(i, gam_temp, data, weights, ic, method)
          score_j <- get_node_score_svy(j, gam_temp, data, weights, ic, method)
          improv <- (ic_best[i] - score_i) + (ic_best[j] - score_j)
          
          # SAFETY CHECK: Ignore if Inf - Inf creates NaN
          if (!is.na(improv) && improv > ic_improv) {
            ic_improv <- improv
            best_move <- list(type = "rev", i=i, j=j, score_i=score_i, score_j=score_j)
          }
        }
      }
    }
    
    # 3. Apply Best Move
    if (ic_improv > 0) {
      i <- best_move$i; j <- best_move$j
      if (best_move$type == "add") { gam[i, j] <- TRUE; ic_best[i] <- best_move$score_i }
      if (best_move$type == "del") { gam[i, j] <- FALSE; ic_best[i] <- best_move$score_i }
      if (best_move$type == "rev") {
        tmp <- gam[i, j]; gam[i, j] <- !tmp; gam[j, i] <- tmp
        ic_best[i] <- best_move$score_i; ic_best[j] <- best_move$score_j
      }
    }
  }
  return(list(gam = gam + 0, ic_best = sum(ic_best)))
}

OrdCD_svy_with_restarts <- function(data, weights, tiers = NULL, gam = NULL, ic = "bic", method = "probit", maxit = 50, nstart = 10) {
  
  best_global_gam <- gam 
  if (is.null(best_global_gam)) best_global_gam <- matrix(0, ncol(data), ncol(data))
  best_global_ic <- Inf
  var_names <- colnames(data)
  
  for (start in 1:nstart) {
    current_gam <- gam
    
    if (!is.null(current_gam)) {
      # A. SKELETON CLEANSER: Force any existing directed edges to obey Tiers
      if (!is.null(tiers)) {
        for (u in 1:ncol(data)) {
          for (v in 1:ncol(data)) {
            if (current_gam[u, v] == 1 && current_gam[v, u] == 0) {
              if (tiers[[var_names[u]]] > tiers[[var_names[v]]]) {
                current_gam[u, v] <- 0
                current_gam[v, u] <- 1 # Flip it forward
              }
            }
          }
        }
      }
      
      # B. Resolve UNDIRECTED edges randomly but strictly respecting Tiers
      undirected_edges <- which(current_gam == 1 & t(current_gam) == 1, arr.ind = TRUE)
      undirected_edges <- undirected_edges[undirected_edges[,1] < undirected_edges[,2], , drop = FALSE]
      
      if (nrow(undirected_edges) > 0) {
        for (k in 1:nrow(undirected_edges)) {
          u <- undirected_edges[k, 1]
          v <- undirected_edges[k, 2]
          
          if (!is.null(tiers) && tiers[[var_names[u]]] < tiers[[var_names[v]]]) {
            current_gam[u, v] <- 1; current_gam[v, u] <- 0
          } else if (!is.null(tiers) && tiers[[var_names[v]]] < tiers[[var_names[u]]]) {
            current_gam[u, v] <- 0; current_gam[v, u] <- 1
          } else {
            if (sample(c(TRUE, FALSE), 1)) {
              current_gam[u, v] <- 1; current_gam[v, u] <- 0
            } else {
              current_gam[u, v] <- 0; current_gam[v, u] <- 1
            }
          }
        }
      }
    }
    
    if (!gRbase::is.DAG(igraph::graph_from_adjacency_matrix(current_gam))) next
    
    res <- OrdCD_svy_greedy(data, weights, tiers, current_gam, ic, method, maxit)
    
    if (!is.na(res$ic_best) && res$ic_best < best_global_ic) {
      best_global_ic <- res$ic_best
      best_global_gam <- res$gam
    }
  }
  
  return(list(gam = best_global_gam, ic_best = best_global_ic))
}

# ------------------------------------------------------------------------------
# 3. MASTER PIPELINE WRAPPER (ALL ORDINAL)
# ------------------------------------------------------------------------------
run_shep_pipeline_ordinal <- function(data, weight_col, tiers, alpha = 0.05, maxit = 50) {
  weights <- data[[weight_col]]
  survey_vars <- data[, setdiff(colnames(data), weight_col)]
  p <- ncol(survey_vars)
  
  # Ensure all variables are factors (ordinal)
  for (var in colnames(survey_vars)) {
    survey_vars[[var]] <- as.factor(survey_vars[[var]])
  }
  
  cat("PHASE 1: Building Copula PC Skeleton...\n")
  Z_matrix <- matrix(NA, nrow = nrow(survey_vars), ncol = p)
  for (j in 1:p) Z_matrix[, j] <- svy_gaussian_copula_transform(as.numeric(survey_vars[[j]]), weights)
  Z_cor <- cov.wt(Z_matrix, wt = weights, cor = TRUE)$cor
  
  pc_res <- pc(suffStat = list(C = Z_cor, n = nrow(survey_vars)), 
               indepTest = gaussCItest, alpha = alpha, 
               labels = colnames(survey_vars), verbose = FALSE)
  
  cpdag_adj <- as(pc_res@graph, "matrix")
  
  cat("PHASE 2: Running Tiered Global OrdCD (with 10 random restarts)...\n")
  obn_res <- OrdCD_svy_with_restarts(data = survey_vars, 
                                     weights = weights, 
                                     tiers = tiers, 
                                     gam = cpdag_adj, 
                                     ic = "bic", 
                                     method = "probit", 
                                     maxit = maxit,
                                     nstart = 10)
  cat("Done.\n")
  return(list(initial_pc = cpdag_adj, final_adj = obn_res$gam, final_bic = obn_res$ic_best))
}