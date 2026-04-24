# ==============================================================================
# File: svy_obn.R
# Description: Survey-Weighted Ordinal Bayesian Network (oBN)
# ==============================================================================
library(survey)
library(gRbase)
library(igraph)

# Helper: Check if adding an edge is admissible (keeps graph acyclic)
admissible <- function(i, j, gam) {
  if (gam[i, j]) return(TRUE)
  gam[i, j] <- 1
  return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam)))
}

# Helper: Check if reversing an edge is admissible
admissible_rev <- function(i, j, gam) {
  if (gam[i, j] == gam[j, i]) return(FALSE)
  tmp <- gam[i, j]
  gam[j, i] <- tmp
  gam[i, j] <- !tmp
  return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam)))
}

# Core: Survey-Weighted Scorer computing Pseudo-BIC / Pseudo-AIC
svy_mypolr_robust <- function(y_target, x_predictors, weights, ic, method, nq_y) {
  if (is.null(x_predictors) || ncol(x_predictors) == 0) {
    dat <- data.frame(Y = y_target)
    form <- Y ~ 1
  } else {
    colnames(x_predictors) <- paste0("X", 1:ncol(x_predictors))
    dat <- data.frame(Y = y_target, x_predictors)
    form <- as.formula(paste("Y ~", paste(colnames(x_predictors), collapse = " + ")))
  }
  
  des <- svydesign(ids = ~1, weights = ~weights, data = dat)
  IC <- Inf
  n_eff <- sum(weights)
  
  if (nq_y > 2) {
    tryCatch({
      fit <- svyolr(form, design = des, method = method)
      k <- length(fit$coefficients) + length(fit$zeta)
      if (ic == "bic") IC <- -2 * fit$deviance + k * log(n_eff)
      if (ic == "aic") IC <- -2 * fit$deviance + 2 * k
    }, error = function(e) {})
  } else {
    tryCatch({
      fit <- svyglm(form, design = des, family = quasibinomial(link = method))
      k <- length(coef(fit))
      if (ic == "bic") IC <- fit$deviance + k * log(n_eff)
      if (ic == "aic") IC <- fit$deviance + 2 * k
    }, error = function(e) {})
  }
  return(IC)
}

# Helper: Get score for a specific node
get_node_score <- function(i, gam, y, weights, ic, method, nq) {
  parents <- which(gam[i, ] > 0)
  if (length(parents) > 0) {
    return(svy_mypolr_robust(y[, i], y[, parents, drop = FALSE], weights, ic, method, nq[i]))
  } else {
    return(svy_mypolr_robust(y[, i], NULL, weights, ic, method, nq[i]))
  }
}

# Main Hill-Climbing Algorithm
svy_obn_greedy <- function(y, weights, gam = NULL, ic = "bic", method = "probit", maxit = 50) {
  q <- ncol(y)
  nq <- sapply(1:q, function(i) nlevels(as.factor(y[, i])))
  if (is.null(gam)) gam <- matrix(FALSE, q, q) else gam <- (gam != 0)
  
  ic_best <- sapply(1:q, function(i) get_node_score(i, gam, y, weights, ic, method, nq))
  
  iter <- 0
  ic_improv <- 1
  
  while (ic_improv > 0 && iter < maxit) {
    iter <- iter + 1
    ic_improv <- -Inf
    best_move <- NULL
    
    # 1. Evaluate Add/Delete
    for (i in 1:q) {
      for (j in setdiff(1:q, i)) {
        if (admissible(i, j, gam)) {
          gam_temp <- gam
          gam_temp[i, j] <- !gam[i, j] # toggle edge
          new_score <- get_node_score(i, gam_temp, y, weights, ic, method, nq)
          improv <- ic_best[i] - new_score
          
          if (improv > ic_improv) {
            ic_improv <- improv
            best_move <- list(type = if(gam[i,j]) "del" else "add", i=i, j=j, score_i=new_score)
          }
        }
      }
    }
    
    # 2. Evaluate Reversals
    for (i in 1:q) {
      for (j in setdiff(1:q, i)) {
        if (admissible_rev(i, j, gam)) {
          gam_temp <- gam
          gam_temp[j, i] <- gam[i, j]
          gam_temp[i, j] <- !gam[i, j]
          
          score_i <- get_node_score(i, gam_temp, y, weights, ic, method, nq)
          score_j <- get_node_score(j, gam_temp, y, weights, ic, method, nq)
          
          improv <- (ic_best[i] - score_i) + (ic_best[j] - score_j)
          
          if (improv > ic_improv) {
            ic_improv <- improv
            best_move <- list(type = "rev", i=i, j=j, score_i=score_i, score_j=score_j)
          }
        }
      }
    }
    
    # 3. Apply Best Move
    if (ic_improv > 0) {
      i <- best_move$i
      j <- best_move$j
      if (best_move$type == "add") { gam[i, j] <- TRUE; ic_best[i] <- best_move$score_i }
      if (best_move$type == "del") { gam[i, j] <- FALSE; ic_best[i] <- best_move$score_i }
      if (best_move$type == "rev") {
        tmp <- gam[i, j]; gam[i, j] <- !tmp; gam[j, i] <- tmp
        ic_best[i] <- best_move$score_i
        ic_best[j] <- best_move$score_j
      }
    }
  }
  return(list(gam = gam + 0, ic_best = sum(ic_best)))
}

#' Main Survey-Weighted oBN Wrapper
svy_obn <- function(y, weights, ic = "bic", link = "probit", nstart = 1, maxit = 50) {
  best_gam <- NULL
  global_best_ic <- Inf
  q <- ncol(y)
  
  for (r in 1:nstart) {
    gam_ini <- if(r == 1) matrix(FALSE, q, q) else bnlearn::random.graph(as.character(1:q), num=1, prob=1/q)
    if(r > 1) gam_ini <- as.matrix(gam_ini$arcs)
    
    fit <- svy_obn_greedy(y, weights, gam = gam_ini, ic = ic, method = link, maxit = maxit)
    if (fit$ic_best < global_best_ic) {
      global_best_ic <- fit$ic_best
      best_gam <- fit$gam
    }
  }
  return(list(gam = best_gam, ic_best = global_best_ic))
}