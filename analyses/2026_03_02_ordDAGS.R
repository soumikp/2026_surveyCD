# ==============================================================================
# Survey-Weighted Covariate-Adjusted Ordinal Causal Discovery (swa-oBN)
# SPEED-OPTIMIZED & ROBUST VERSION
# ==============================================================================

admissible = function(i, j, gam_short_old) {
  if (gam_short_old[i, j]) {
    return(TRUE)
  } else {
    gam_short_old[i, j] = 1
    return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam_short_old)))
  }
}

admissible_rev = function(i, j, gam_short_old) {
  if (gam_short_old[i, j] == gam_short_old[j, i]) {
    return(FALSE)
  } else {
    tmp = gam_short_old[i, j]
    gam_short_old[j, i] = tmp
    gam_short_old[i, j] = !tmp
    return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam_short_old)))
  }
}

# ------------------------------------------------------------------------------
# Core Regression Wrapper (SPEED OPTIMIZED & ROBUST)
# ------------------------------------------------------------------------------
mypolr = function(y_vec, x_df, z_df, weights, ic, method, nq_y) {
  
  # 1. Fast, robust data binding (Avoids cbind() dispatch bugs with NULLs)
  data_list = list(Y = y_vec)
  if (!is.null(x_df)) data_list = c(data_list, as.list(x_df))
  if (!is.null(z_df)) data_list = c(data_list, as.list(z_df))
  df = as.data.frame(data_list)
  
  # 2. Fast formula generation (Y ~ . is much faster than paste parsing)
  form = stats::as.formula("Y ~ .")
  
  calc_ic <- function(model) {
    ll <- as.numeric(stats::logLik(model))
    k <- length(stats::coef(model))
    if (inherits(model, "polr")) {
      k <- k + length(model$zeta)
    }
    n <- nrow(df) 
    
    if (ic == "bic") {
      return(-2 * ll + k * log(n))
    } else {
      return(-2 * ll + 2 * k)
    }
  }
  
  boolFalse = FALSE
  IC = Inf
  
  # 3. SPEED FIX: Reduce max_tries to 1. 
  # If a model fails to fit cleanly on the first or second try, 
  # it's a statistically bad edge. Skip it and save time!
  tries = 0
  max_tries = 1 
  
  if (nq_y > 2) {
    tryCatch({
      # Hess = FALSE saves computation time at the end of the BFGS optimization
      mod = MASS::polr(form, data = df, weights = weights, method = method, Hess = FALSE)
      IC = calc_ic(mod)
      boolFalse <- TRUE
    }, error = function(e) {})
    
    while (!boolFalse && tries < max_tries) {
      tries = tries + 1
      tryCatch({
        X_mat = stats::model.matrix(form, df)
        n_coef = ncol(X_mat) - 1 
        n_zeta = nq_y - 1
        mod = MASS::polr(form, data = df, weights = weights, method = method, Hess = FALSE,
                         start = sort(stats::rnorm(n_coef + n_zeta)))
        IC = calc_ic(mod)
        boolFalse <- TRUE
      }, error = function(e) {})
    }
  } else {
    tryCatch({
      mod = stats::glm(form, data = df, weights = weights, family = stats::binomial(link = method))
      IC = calc_ic(mod)
      boolFalse <- TRUE
    }, error = function(e) {})
    
    while (!boolFalse && tries < max_tries) {
      tries = tries + 1
      tryCatch({
        X_mat = stats::model.matrix(form, df)
        n_coef = ncol(X_mat)
        mod = stats::glm(form, data = df, weights = weights, family = stats::binomial(link = method),
                         start = stats::rnorm(n_coef))
        IC = calc_ic(mod)
        boolFalse <- TRUE
      }, error = function(e) {})
    }
  }
  
  # Prevent NA propagation if the model completely fails
  if (is.na(IC)) IC = Inf 
  return(IC)
}

# ------------------------------------------------------------------------------
# Data Preparation Wrapper 
# ------------------------------------------------------------------------------
mypolr_wrap = function(i, gam_i, y, Z, weights, ic, method, nq) {
  y_vec = y[, i]
  nq_y = nq[i]
  
  x_df = NULL
  if (sum(gam_i) > 0) {
    x_df = y[, gam_i, drop = FALSE]
    colnames(x_df) = paste0("X", which(gam_i))
  }
  
  z_df = NULL
  if (!is.null(Z)) {
    z_df = as.data.frame(Z) # Ensure Z is a dataframe
    if (is.null(colnames(z_df))) colnames(z_df) = paste0("Z", 1:ncol(z_df))
  }
  
  return(mypolr(y_vec, x_df, z_df, weights, ic, method, nq_y))
}

# ------------------------------------------------------------------------------
# Greedy Search Algorithm
# ------------------------------------------------------------------------------
oBN_greedy = function(y, Z = NULL, weights = NULL, gam = NULL, ic = "bic", method = "probit", verbose = FALSE, maxit = 50) {
  n = nrow(y)
  q = ncol(y)
  nq = sapply(1:q, function(i) nlevels(y[, i]))
  
  if (is.null(gam)) {
    gam = matrix(FALSE, q, q)
  } else {
    gam = (gam != 0)
  }
  
  ind_q = matrix(0, q, q - 1)
  for (i in 1:q) {
    if (i == 1) ind_noi = 2:q
    else if (i == q) ind_noi = 1:(q - 1)
    else ind_noi = c(1:(i - 1), (i + 1):q)
    ind_q[i, ] = ind_noi
  }
  
  iter = 0
  ic_improv = 1
  act_ind = c(NA, NA)
  state = "add"
  
  ic_best = rep(0, q)
  for (i in 1:q) {
    ic_best[i] = mypolr_wrap(i, gam[i, ], y, Z, weights, ic, method, nq)
  }
  
  while (ic_improv > 0 && iter < maxit) {
    iter = iter + 1
    ic_improv = -Inf
    ic_improv_rev = rep(-Inf, 2)
    gam_new = gam
    ic_improv_new = -Inf
    ic_improv_rev_new = rep(-Inf, 2)
    ic_best_new = -Inf
    ic_rev_best_new = rep(-Inf, 2)
    
    for (i in 1:q) {
      for (j in 1:(q - 1)) {
        if (admissible(i, ind_q[i, j], gam)) {
          is_delete = gam[i, ind_q[i, j]]
          
          gam_new[i, ind_q[i, j]] = !is_delete
          ic_best_new = mypolr_wrap(i, gam_new[i, ], y, Z, weights, ic, method, nq)
          
          ic_improv_new = ic_best[i] - ic_best_new
          if (is.na(ic_improv_new)) ic_improv_new = -Inf # Prevent NA failure
          
          if (ic_improv_new > ic_improv) {
            ic_improv = ic_improv_new
            act_ind = c(i, ind_q[i, j])
            state = ifelse(is_delete, "del", "add")
          }
          gam_new[i, ind_q[i, j]] = is_delete
        }
      }
    }
    
    for (i in 1:q) {
      for (j in 1:(q - 1)) {
        if (admissible_rev(i, ind_q[i, j], gam)) {
          tmp = gam_new[i, ind_q[i, j]]
          gam_new[ind_q[i, j], i] = tmp
          gam_new[i, ind_q[i, j]] = !tmp
          
          ic_rev_best_new[1] = mypolr_wrap(i, gam_new[i, ], y, Z, weights, ic, method, nq)
          ic_rev_best_new[2] = mypolr_wrap(ind_q[i, j], gam_new[ind_q[i, j], ], y, Z, weights, ic, method, nq)
          
          ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
          ic_improv_rev_new[2] = ic_best[ind_q[i, j]] - ic_rev_best_new[2]
          
          ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
          if (is.na(ic_improv_new)) ic_improv_new = -Inf # Prevent NA failure
          
          if (ic_improv_new > ic_improv) {
            ic_improv = ic_improv_new
            ic_improv_rev = ic_improv_rev_new
            act_ind = c(i, ind_q[i, j])
            state = "rev"
          }
          gam_new[i, ind_q[i, j]] = tmp
          gam_new[ind_q[i, j], i] = !tmp
        }
      }
    }
    
    if (ic_improv > 0) {
      if (state == "add") {
        gam[act_ind[1], act_ind[2]] = TRUE
        ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv
      } else if (state == "del") {
        gam[act_ind[1], act_ind[2]] = FALSE
        ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv
      } else if (state == "rev") {
        tmp = gam[act_ind[1], act_ind[2]]
        gam[act_ind[2], act_ind[1]] = tmp
        gam[act_ind[1], act_ind[2]] = !tmp
        ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv_rev[1]
        ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv_rev[2]
      }
    }
    
    if (verbose && iter %% 1 == 0) {
      print(paste(iter, " iterations have completed", sep = ""))
      print("The current DAG adjacency matrix is")
      print(gam + 0)
      print(paste("with ", ic, " = ", sum(ic_best), sep = ""))
    }
  }
  
  if (iter == maxit) {
    warning("The maximum number of iterations was reached. The algorithm has not converged.")
  }
  return(list(gam = gam + 0, ic_best = sum(ic_best)))
}

# ------------------------------------------------------------------------------
# Greedy Search CPDAG Wrapper 
# ------------------------------------------------------------------------------
oBN_greedy_CPDAG = function(y, Z = NULL, weights = NULL, gam = NULL, ic = "bic", edge_list = NULL, method = "probit", verbose = FALSE, maxit = 50) {
  n = nrow(y)
  q = ncol(y)
  nq = sapply(1:q, function(i) nlevels(y[, i]))
  
  if (is.null(gam)) {
    gam = matrix(FALSE, q, q)
  } else {
    gam = (gam != 0)
  }
  
  ind_q = vector("list", q)
  nind_q = rep(0, q)
  for (e in 1:nrow(edge_list)) {
    i = edge_list[e, 1]
    ind_q[[i]] = c(ind_q[[i]], edge_list[e, 2])
    nind_q[i] = nind_q[i] + 1
  }
  
  iter = 0
  ic_improv = 1
  act_ind = c(NA, NA)
  
  ic_best = rep(0, q)
  for (i in 1:q) {
    ic_best[i] = mypolr_wrap(i, gam[i, ], y, Z, weights, ic, method, nq)
  }
  
  while (ic_improv > 0 && iter < maxit) {
    iter = iter + 1
    ic_improv = -Inf
    ic_improv_rev = rep(-Inf, 2)
    gam_new = gam
    ic_improv_new = -Inf
    ic_improv_rev_new = rep(-Inf, 2)
    ic_rev_best_new = rep(-Inf, 2)
    
    for (i in 1:q) {
      if (nind_q[i] > 0) {
        for (j in 1:(nind_q[i])) {
          if (admissible_rev(i, ind_q[[i]][j], gam)) {
            tmp = gam_new[i, ind_q[[i]][j]]
            gam_new[ind_q[[i]][j], i] = tmp
            gam_new[i, ind_q[[i]][j]] = !tmp
            
            ic_rev_best_new[1] = mypolr_wrap(i, gam_new[i, ], y, Z, weights, ic, method, nq)
            ic_rev_best_new[2] = mypolr_wrap(ind_q[[i]][j], gam_new[ind_q[[i]][j], ], y, Z, weights, ic, method, nq)
            
            ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
            ic_improv_rev_new[2] = ic_best[ind_q[[i]][j]] - ic_rev_best_new[2]
            ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
            
            if (is.na(ic_improv_new)) ic_improv_new = -Inf # Prevent NA failure
            
            if (ic_improv_new > ic_improv) {
              ic_improv = ic_improv_new
              ic_improv_rev = ic_improv_rev_new
              act_ind = c(i, ind_q[[i]][j])
            }
            gam_new[i, ind_q[[i]][j]] = tmp
            gam_new[ind_q[[i]][j], i] = !tmp
          }
        }
      }
    }
    
    if (ic_improv > 0) {
      tmp = gam[act_ind[1], act_ind[2]]
      gam[act_ind[2], act_ind[1]] = tmp
      gam[act_ind[1], act_ind[2]] = !tmp
      ic_best[act_ind[1]] = ic_best[act_ind[1]] - ic_improv_rev[1]
      ic_best[act_ind[2]] = ic_best[act_ind[2]] - ic_improv_rev[2]
    }
    
    if (verbose && iter %% 1 == 0) {
      print(paste(iter, " iterations have completed", sep = ""))
      print("The current DAG adjacency matrix is")
      print(gam + 0)
      print(paste("with ", ic, " = ", sum(ic_best), sep = ""))
    }
  }
  
  if (iter == maxit) {
    warning("The maximum number of iterations was reached. The algorithm has not converged.")
  }
  return(list(gam = gam + 0, ic_best = sum(ic_best)))
}

# ------------------------------------------------------------------------------
# Initialization & Random Restarts Wrapper
# ------------------------------------------------------------------------------
oBN_greedy_wrap = function(y, Z = NULL, weights = NULL, ic = "bic", edge_list = NULL, method = "probit", nstart = 1, verbose = FALSE, maxit = 50) {
  q = ncol(y)
  gam_list  = vector("list", nstart)
  ic_best_list = rep(NA, nstart)
  
  if (nstart == 1) {
    if (is.null(edge_list)) {
      fit = oBN_greedy(y, Z, weights, gam = NULL, ic = ic, method = method, verbose = verbose, maxit = maxit)
    } else {
      gam = matrix(0, q, q)
      gam[edge_list] = 1
      und_edge = which(as.matrix(Matrix::tril(gam * t(gam))) == 1, arr.ind = TRUE)
      for (i in 1:nrow(und_edge)) {
        if (stats::rbinom(1, 1, .5) == 1) {
          gam[und_edge[i, 1], und_edge[i, 2]] = 0
        } else {
          gam[und_edge[i, 2], und_edge[i, 1]] = 0
        }
      }
      fit = oBN_greedy_CPDAG(y, Z, weights, gam = gam, ic = ic, edge_list = und_edge, method = method, verbose = verbose, maxit = maxit)
    }
    gam_list[[1]] = fit$gam
    ic_best_list[1] = fit$ic_best
  } else {
    if (is.null(edge_list)) {
      netlist = bnlearn::random.graph(nodes = as.character(1:q), method = "ordered", num = nstart - 1, prob = 1 / q)
      if (nstart == 2) netlist = list(netlist)
      
      for (i in 1:nstart) {
        gam = matrix(FALSE, q, q)
        if (i != 1) {
          gam[apply(netlist[[i - 1]]$arcs, 2, as.numeric)] = TRUE
        }
        fit = oBN_greedy(y, Z, weights, gam = gam, ic = ic, method = method, verbose = verbose, maxit = maxit)
        gam_list[[i]] = fit$gam
        ic_best_list[i] = fit$ic_best
      }
    } else {
      gam = matrix(0, q, q)
      gam[edge_list] = 1
      und_edge = which(as.matrix(Matrix::tril(gam * t(gam))) == 1, arr.ind = TRUE)
      for (i in 1:nstart) {
        gam_ini = gam
        for (ii in 1:nrow(und_edge)) {
          if (stats::rbinom(1, 1, .5) == 1) {
            gam_ini[und_edge[ii, 1], und_edge[ii, 2]] = 0
          } else {
            gam_ini[und_edge[ii, 2], und_edge[ii, 1]] = 0
          }
        }
        fit = oBN_greedy_CPDAG(y, Z, weights, gam = gam_ini, ic = ic, edge_list = und_edge, method = method, verbose = verbose, maxit = maxit)
        gam_list[[i]] = fit$gam
        ic_best_list[i] = fit$ic_best
      }
    }
  }
  
  i = which.min(ic_best_list)
  return(list(gam = gam_list[[i]], ic_best = ic_best_list[i]))
}

# ------------------------------------------------------------------------------
# Exhaustive Search (Small N Only)
# ------------------------------------------------------------------------------
oBN_exhaust = function(y, Z = NULL, weights = NULL, gam_list = NULL, ic = "bic", method = "probit") {
  n = nrow(y)
  q = ncol(y)
  nl = sapply(1:q, function(i) nlevels(y[, i]))
  
  if (is.null(gam_list)) {
    if (q == 2) {
      gam_list = array(0, c(q, q, 3))
      gam_list[1, 2, 2] = 1
      gam_list[2, 1, 3] = 1
    } else if (q == 3) {
      gam_list = array(0, c(q, q, 25))
      gam_list[1, 2, 2] = 1
      gam_list[1, 3, 3] = 1
      gam_list[2, 3, 4] = 1
      gam_list[2, 1, 5] = 1
      gam_list[3, 1, 6] = 1
      gam_list[3, 2, 7] = 1
      gam_list[1, 2, 8] = gam_list[1, 3, 8] = 1
      gam_list[2, 3, 9] = gam_list[2, 1, 9] = 1
      gam_list[3, 1, 10] = gam_list[3, 2, 10] = 1
      gam_list[2, 1, 11] = gam_list[3, 1, 11] = 1
      gam_list[3, 2, 12] = gam_list[1, 2, 12] = 1
      gam_list[1, 3, 13] = gam_list[2, 3, 13] = 1
      gam_list[1, 2, 14] = gam_list[2, 3, 14] = 1
      gam_list[1, 3, 15] = gam_list[3, 2, 15] = 1
      gam_list[2, 3, 16] = gam_list[3, 1, 16] = 1
      gam_list[2, 1, 17] = gam_list[1, 3, 17] = 1
      gam_list[3, 2, 18] = gam_list[2, 1, 18] = 1
      gam_list[3, 1, 19] = gam_list[1, 2, 19] = 1
      gam_list[1, 2, 20] = gam_list[1, 3, 20] = gam_list[2, 3, 20] = 1
      gam_list[1, 2, 21] = gam_list[1, 3, 21] = gam_list[3, 2, 21] = 1
      gam_list[2, 1, 22] = gam_list[2, 3, 22] = gam_list[1, 3, 22] = 1
      gam_list[2, 1, 23] = gam_list[2, 3, 23] = gam_list[3, 1, 23] = 1
      gam_list[3, 1, 24] = gam_list[3, 2, 24] = gam_list[1, 2, 24] = 1
      gam_list[3, 1, 25] = gam_list[3, 2, 25] = gam_list[2, 1, 25] = 1
    } else {
      stop("The number of nodes must be 2 or 3")
    }
  }
  
  IC = rep(0, dim(gam_list)[3])
  for (m in 1:length(IC)) {
    gam = gam_list[, , m]
    for (i in 1:q) {
      gam_tmp = gam[i, ]
      IC[m] = IC[m] + mypolr_wrap(i, gam_tmp, y, Z, weights, ic, method, nl)
    }
  }
  
  mi = which.min(IC)
  return(list(gam = gam_list[, , mi] + 0, ic_best = IC[mi]))
}

# ------------------------------------------------------------------------------
# Main Dispatch Functions
# ------------------------------------------------------------------------------
OCD = function(y, Z = NULL, weights = NULL, search = "greedy", ic = "bic", edge_list = NULL, link = "probit", G = NULL, nstart = 1, verbose = FALSE, maxit = 50) {
  if (search == "exhaust") {
    G = oBN_exhaust(y, Z, weights, G, ic, link)
  } else {
    G = oBN_greedy_wrap(y, Z, weights, ic, edge_list, link, nstart, verbose, maxit)
  }
  return(G)
}

#' @title Survey-Weighted Covariate-Adjusted Ordinal Causal Discovery
#'
#' @description Estimate a causal directed acyclic graph (DAG) for ordinal categorical data with greedy or exhaustive search, incorporating sample weights and fixed covariates.
#'
#' @param y a data frame with each column being an ordinal categorical variable, which must be a factor.
#' @param Z an optional data frame of covariates (continuous or categorical). These nodes are locked in the background and are automatically adjusted for in every evaluation.
#' @param weights an optional numeric vector of survey weights. 
#' @param search the search method used to find the best-scored DAG. "greedy" or "exhaust".
#' @param ic the information criterion (AIC or BIC) used to score DAGs. 
#' @param edge_list an edge list of a CPDAG.
#' @param link the link function for ordinal regression. Default is "probit".
#' @param G a list of DAG adjacency matrices for "exhaust" search.
#' @param nstart number of random graph initializations for the "greedy" search.
#' @param verbose if TRUE, messages are printed during the run of the greedy search algorithm.
#' @param maxit the maximum number of iterations for the greedy search algorithm.
#' @param boot the number of bootstrap samples. 
#' @return A list with "boot" elements (or 1 if no boot). Each element is a list with 'gam' and 'ic_best'.
#'
#' @export
OrdCD = function(y,
                 Z = NULL,
                 weights = NULL,
                 search = "greedy",
                 ic = "bic",
                 edge_list = NULL,
                 link = "probit",
                 G = NULL,
                 nstart = 1,
                 verbose = FALSE,
                 maxit = 50,
                 boot = NULL) {
  
  # AUTO-COERCE DATAFRAMES TO PROPER FORMATS
  if (!is.null(weights)) {
    if (is.data.frame(weights) || is.matrix(weights)) {
      weights = as.numeric(weights[, 1]) # Extract as vector
    }
  } else {
    weights = rep(1, nrow(y))
  }
  
  if (!is.null(Z)) {
    Z = as.data.frame(Z)
  }
  
  if (is.null(boot)) {
    G_out = OCD(y, Z, weights, search, ic, edge_list, link, G, nstart, verbose, maxit)
  } else {
    G_out = vector("list", boot)
    for (b in 1:boot) {
      idx = sample(nrow(y), replace = TRUE)
      
      y_boot = y[idx, , drop = FALSE]
      w_boot = weights[idx]
      Z_boot = if (!is.null(Z)) Z[idx, , drop = FALSE] else NULL
      
      G_out[[b]] = OCD(y_boot, Z_boot, w_boot, search, ic, edge_list, link, G, nstart, verbose, maxit)
    }
  }
  return(G_out)
}

library(igraph)
library(ggraph)

plot_directed_graph <- function(conn_matrix, node_names = NULL, 
                                layout_type = 'fr', 
                                node_color = "#4C72B0", 
                                plot_title = "Directed Connectivity Graph") {
  
  # 1. Handle node names
  if (!is.null(node_names)) {
    # Check if the length of provided names matches the matrix dimensions
    if (length(node_names) != ncol(conn_matrix)) {
      stop("Error: The length of 'node_names' must match the number of columns in the matrix.")
    }
    colnames(conn_matrix) <- node_names
    rownames(conn_matrix) <- node_names
  } else if (is.null(colnames(conn_matrix))) {
    # Fallback: assign generic names if none exist
    generic_names <- paste0("Node_", 1:ncol(conn_matrix))
    colnames(conn_matrix) <- generic_names
    rownames(conn_matrix) <- generic_names
  }
  
  # 2. Convert to igraph object
  graph_obj <- graph_from_adjacency_matrix(
    conn_matrix, 
    mode = "directed", 
    diag = FALSE # Set to TRUE if you want self-loops
  )
  
  # 3. Generate the plot
  p <- ggraph(graph_obj, layout = layout_type) + 
    geom_edge_link(
      arrow = arrow(length = unit(4, 'mm'), type = "closed"),
      end_cap = circle(8, 'mm'), 
      color = "darkgray",
      width = 0.8
    ) +
    geom_node_point(
      size = 16, 
      color = node_color, 
      alpha = 0.9
    ) +
    geom_node_text(
      aes(label = name), 
      color = "black", 
      fontface = "bold",
      size = 4
    ) +
    theme_graph(base_family = "sans") +
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Return the plot object so it can be viewed or saved
  return(p)
}

