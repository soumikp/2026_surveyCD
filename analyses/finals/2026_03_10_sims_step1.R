# ==============================================================================
# Survey-Weighted Covariate-Adjusted Ordinal Causal Discovery (swa-oBN)
# SPEED-OPTIMIZED, ROBUST, WITH HASH-MAP SCORE CACHING
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
# Core Regression Wrapper 
# ------------------------------------------------------------------------------
mypolr = function(y_vec, x_df, z_df, weights, ic, method, nq_y) {
  
  data_list = list(Y = y_vec)
  pred_names = character(0)
  
  if (!is.null(x_df) && ncol(x_df) > 0) {
    data_list = c(data_list, as.list(x_df))
    pred_names = c(pred_names, colnames(x_df))
  }
  if (!is.null(z_df) && ncol(z_df) > 0) {
    data_list = c(data_list, as.list(z_df))
    pred_names = c(pred_names, colnames(z_df))
  }
  
  df = as.data.frame(data_list)
  
  # SAFE WEIGHT INJECTION: Put weights directly inside df to avoid R scoping bugs
  df$.weights = weights
  
  # Explicitly build formula so .weights is NOT treated as a predictor
  if (length(pred_names) > 0) {
    form = stats::as.formula(paste("Y ~", paste(pred_names, collapse = " + ")))
  } else {
    form = stats::as.formula("Y ~ 1")
  }
  
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
  tries = 0
  max_tries = 1 
  
  if (nq_y > 2) {
    tryCatch({
      mod = MASS::polr(form, data = df, weights = .weights, method = method, Hess = FALSE)
      IC = calc_ic(mod)
      boolFalse <- TRUE
    }, error = function(e) {})
    
    while (!boolFalse && tries < max_tries) {
      tries = tries + 1
      tryCatch({
        X_mat = stats::model.matrix(form, df)
        n_coef = ncol(X_mat) - 1 
        n_zeta = nq_y - 1
        mod = MASS::polr(form, data = df, weights = .weights, method = method, Hess = FALSE,
                         start = sort(stats::rnorm(n_coef + n_zeta)))
        IC = calc_ic(mod)
        boolFalse <- TRUE
      }, error = function(e) {})
    }
  } else {
    tryCatch({
      mod = stats::glm(form, data = df, weights = .weights, family = stats::binomial(link = method))
      IC = calc_ic(mod)
      boolFalse <- TRUE
    }, error = function(e) {})
    
    while (!boolFalse && tries < max_tries) {
      tries = tries + 1
      tryCatch({
        X_mat = stats::model.matrix(form, df)
        n_coef = ncol(X_mat)
        mod = stats::glm(form, data = df, weights = .weights, family = stats::binomial(link = method),
                         start = stats::rnorm(n_coef))
        IC = calc_ic(mod)
        boolFalse <- TRUE
      }, error = function(e) {})
    }
  }
  
  if (is.na(IC)) IC = Inf 
  return(IC)
}
# ------------------------------------------------------------------------------
# Data Preparation & CACHE Wrapper 
# ------------------------------------------------------------------------------
mypolr_wrap = function(i, gam_i, y, Z, weights, ic, method, nq, score_cache = NULL) {
  
  # Ensure gam_i is a logical vector (Exhaustive search passes 0/1 numerics)
  gam_i = (gam_i > 0)
  
  # --- CACHE CHECK ---
  parent_str = paste(which(gam_i), collapse = "_")
  cache_key = paste0("node_", i, "_parents_", parent_str)
  
  if (!is.null(score_cache) && exists(cache_key, envir = score_cache)) {
    return(get(cache_key, envir = score_cache))
  }
  # -------------------
  
  y_vec = y[, i]
  nq_y = nq[i]
  
  x_df = NULL
  if (any(gam_i)) {
    x_df = y[, gam_i, drop = FALSE]
    colnames(x_df) = paste0("X", which(gam_i))
  }
  
  z_df = NULL
  if (!is.null(Z)) {
    z_df = as.data.frame(Z) 
    if (is.null(colnames(z_df))) colnames(z_df) = paste0("Z", 1:ncol(z_df))
  }
  
  IC = mypolr(y_vec, x_df, z_df, weights, ic, method, nq_y)
  
  # --- SAVE TO CACHE ---
  if (!is.null(score_cache)) {
    assign(cache_key, IC, envir = score_cache)
  }
  # ---------------------
  
  return(IC)
}

# ------------------------------------------------------------------------------
# Greedy Search Algorithm
# ------------------------------------------------------------------------------
oBN_greedy = function(y, Z = NULL, weights = NULL, gam = NULL, ic = "bic", method = "probit", verbose = FALSE, maxit = 50) {
  n = nrow(y)
  q = ncol(y)
  nq = sapply(1:q, function(i) nlevels(y[, i]))
  
  # Initialize the Hash Map for caching scores
  score_cache <- new.env(hash = TRUE)
  
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
    ic_best[i] = mypolr_wrap(i, gam[i, ], y, Z, weights, ic, method, nq, score_cache)
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
          ic_best_new = mypolr_wrap(i, gam_new[i, ], y, Z, weights, ic, method, nq, score_cache)
          
          ic_improv_new = ic_best[i] - ic_best_new
          if (is.na(ic_improv_new)) ic_improv_new = -Inf 
          
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
          
          ic_rev_best_new[1] = mypolr_wrap(i, gam_new[i, ], y, Z, weights, ic, method, nq, score_cache)
          ic_rev_best_new[2] = mypolr_wrap(ind_q[i, j], gam_new[ind_q[i, j], ], y, Z, weights, ic, method, nq, score_cache)
          
          ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
          ic_improv_rev_new[2] = ic_best[ind_q[i, j]] - ic_rev_best_new[2]
          
          ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
          if (is.na(ic_improv_new)) ic_improv_new = -Inf 
          
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
      print(paste("Cache size:", length(ls(score_cache)), "unique models evaluated."))
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
  
  score_cache <- new.env(hash = TRUE)
  
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
    ic_best[i] = mypolr_wrap(i, gam[i, ], y, Z, weights, ic, method, nq, score_cache)
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
            
            ic_rev_best_new[1] = mypolr_wrap(i, gam_new[i, ], y, Z, weights, ic, method, nq, score_cache)
            ic_rev_best_new[2] = mypolr_wrap(ind_q[[i]][j], gam_new[ind_q[[i]][j], ], y, Z, weights, ic, method, nq, score_cache)
            
            ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
            ic_improv_rev_new[2] = ic_best[ind_q[[i]][j]] - ic_rev_best_new[2]
            ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
            
            if (is.na(ic_improv_new)) ic_improv_new = -Inf 
            
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
      print(paste("Cache size:", length(ls(score_cache)), "unique models evaluated."))
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
      IC[m] = IC[m] + mypolr_wrap(i, gam_tmp, y, Z, weights, ic, method, nl, score_cache = NULL)
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
  
  if (!is.null(weights)) {
    if (is.data.frame(weights) || is.matrix(weights)) {
      weights = as.numeric(weights[, 1]) 
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