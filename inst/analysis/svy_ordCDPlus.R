# =============================================================================
# ordCD_updated.R
# Extended Ordinal Causal Discovery
# New features:
#   1. Allow 2-level ordinal variables (binary)
#   2. Survey weights with effective sample size in BIC
#   3. Confounders Z (mixed numeric/factor, not in DAG)
#   4. Blacklist edges (directed)
# =============================================================================

# --- Blacklist helper --------------------------------------------------------
# Convert blacklist (data.frame or matrix with "from","to") to a logical matrix
.make_blacklist_mat = function(blacklist, q, colnames_y = NULL) {
  bl = matrix(FALSE, q, q)
  if (is.null(blacklist)) return(bl)
  blacklist = as.data.frame(blacklist, stringsAsFactors = FALSE)
  stopifnot(all(c("from", "to") %in% names(blacklist)))
  fr = blacklist$from
  to = blacklist$to
  # Support both numeric indices and column names
  if (is.character(fr) && !is.null(colnames_y)) {
    fr = match(fr, colnames_y)
    to = match(to, colnames_y)
    if (any(is.na(fr)) || any(is.na(to)))
      stop("Blacklist contains names not found in column names of y.")
  }
  fr = as.integer(fr)
  to = as.integer(to)
  for (k in seq_along(fr)) {
    # bl[i,j]=TRUE means edge j->i is blacklisted
    # "from" = j (parent), "to" = i (child), so bl[to, from]
    bl[to[k], fr[k]] = TRUE
  }
  return(bl)
}

# --- Effective sample size (Kish) -------------------------------------------
.eff_n = function(w) {
  if (is.null(w)) return(NULL)
  (sum(w))^2 / sum(w^2)
}

# --- Custom weighted BIC/AIC ------------------------------------------------
# Uses effective sample size for penalty when weights are present
.weighted_ic = function(fit, ic, n_eff) {
  ll = as.numeric(stats::logLik(fit))
  k = attr(stats::logLik(fit), "df")
  if (ic == "bic") {
    return(-2 * ll + k * log(n_eff))
  } else {
    return(-2 * ll + 2 * k)
  }
}

# --- DAG admissibility checks -----------------------------------------------
admissible = function(i, j, gam_short_old, bl_mat) {
  # Check if adding/keeping edge j->i is admissible
  if (gam_short_old[i, j]) {
    # Deleting an edge is always admissible
    return(TRUE)
  } else {
    # Adding j->i: check blacklist first
    if (bl_mat[i, j]) return(FALSE)
    gam_short_old[i, j] = 1
    return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(gam_short_old)))
  }
}

admissible_rev = function(i, j, gam_short_old, bl_mat) {
  if (gam_short_old[i, j] == gam_short_old[j, i]) {
    return(FALSE)
  } else {
    # Determine the reversal: if i->j exists, we reverse to j->i (and vice versa)
    # After reversal, new edges must not be blacklisted
    tmp_ij = gam_short_old[i, j]
    tmp_ji = gam_short_old[j, i]
    # After swap: gam[j,i] = tmp_ij, gam[i,j] = tmp_ji
    # The new edge being created: if tmp_ij was 1 (i.e., j->i existed),
    # after reversal gam[j,i]=1 means i->j is created
    # We need to check the NEW directed edge is not blacklisted
    new_gam = gam_short_old
    new_gam[j, i] = tmp_ij
    new_gam[i, j] = tmp_ji
    # Identify the newly created edge direction
    if (tmp_ij && !tmp_ji) {
      # Was j->i (gam[i,j]=1), now becomes i->j (gam[j,i]=1)
      if (bl_mat[j, i]) return(FALSE)
    } else if (tmp_ji && !tmp_ij) {
      # Was i->j (gam[j,i]=1), now becomes j->i (gam[i,j]=1)
      if (bl_mat[i, j]) return(FALSE)
    }
    return(gRbase::is.DAG(igraph::graph_from_adjacency_matrix(new_gam)))
  }
}

# --- Ordinal / binary regression wrapper ------------------------------------
# Now supports: weights, confounders Z, effective sample size for IC
mypolr = function(formula, data, ic, method, nq_y, nq_x,
                  w = NULL, n_eff = NULL) {
  boolFalse = FALSE
  # Use raw n if no weights
  if (is.null(n_eff)) n_eff = nrow(data)
  
  if (nq_y > 2) {
    # --- Ordinal response (>=3 levels) ---
    tryCatch({
      if (is.null(w)) {
        fit = MASS::polr(formula, data = data, method = method)
      } else {
        fit = MASS::polr(formula, data = data, method = method, weights = w)
      }
      IC = .weighted_ic(fit, ic, n_eff)
      boolFalse <- TRUE
    }, error = function(e) {})
    while (!boolFalse) {
      tryCatch({
        npar = nq_y - 1 + sum(nq_x - 1)
        # Account for Z columns already in data: polr will estimate their params too
        # but npar here is just for the start vector length; polr figures out the rest
        # We need total number of predictors in formula
        # Safer: let polr estimate, just provide random starts for threshold+beta
        start_len = nq_y - 1 + (ncol(data) - 1)  # rough: thresholds + all predictors
        # But we only need to match what polr expects
        # polr start = c(coefficients, zeta); length = n_predictors + (nq_y - 1)
        n_pred = ncol(data) - 1  # all columns except response
        start_vec = c(sort(stats::rnorm(n_pred)), sort(stats::rnorm(nq_y - 1)))
        if (is.null(w)) {
          fit = MASS::polr(formula, data = data, start = start_vec, method = method)
        } else {
          fit = MASS::polr(formula, data = data, start = start_vec,
                           method = method, weights = w)
        }
        IC = .weighted_ic(fit, ic, n_eff)
        boolFalse <- TRUE
      }, error = function(e) {})
    }
  } else {
    # --- Binary response (2 levels) ---
    tryCatch({
      if (is.null(w)) {
        fit = stats::glm(formula, data = data,
                         family = stats::binomial(link = method))
      } else {
        fit = stats::glm(formula, data = data,
                         family = stats::binomial(link = method), weights = w)
      }
      IC = .weighted_ic(fit, ic, n_eff)
      boolFalse <- TRUE
    }, error = function(e) {})
    while (!boolFalse) {
      tryCatch({
        n_pred = ncol(data) - 1
        start_vec = sort(stats::rnorm(1 + n_pred))
        if (is.null(w)) {
          fit = stats::glm(formula, data = data,
                           start = start_vec,
                           family = stats::binomial(link = method))
        } else {
          fit = stats::glm(formula, data = data,
                           start = start_vec,
                           family = stats::binomial(link = method), weights = w)
        }
        IC = .weighted_ic(fit, ic, n_eff)
        boolFalse <- TRUE
      }, error = function(e) {})
    }
  }
  return(IC)
}

# --- Null model IC (intercept only, possibly with Z) -------------------------
.null_ic = function(y_col, nq_i, ic, method, w, n_eff, Z) {
  if (is.null(n_eff)) n_eff = length(y_col)
  if (is.null(Z)) {
    dat_null = data.frame(.Y = y_col)
    frm = .Y ~ 1
  } else {
    dat_null = data.frame(.Y = y_col, Z)
    frm = stats::as.formula(paste(".Y ~", paste(names(Z), collapse = " + ")))
  }
  if (nq_i > 2) {
    if (is.null(w)) {
      fit = MASS::polr(frm, data = dat_null, method = method)
    } else {
      fit = MASS::polr(frm, data = dat_null, method = method, weights = w)
    }
  } else {
    if (is.null(w)) {
      fit = stats::glm(frm, data = dat_null,
                       family = stats::binomial(link = method))
    } else {
      fit = stats::glm(frm, data = dat_null,
                       family = stats::binomial(link = method), weights = w)
    }
  }
  return(.weighted_ic(fit, ic, n_eff))
}

# --- Build regression data for node i given parent mask ----------------------
# Combines parent columns from y with confounder columns from Z
.build_reg_data = function(y, i, parent_mask, Z) {
  # parent_mask is a logical vector of length q
  if (is.null(Z)) {
    if (sum(parent_mask) > 0) {
      dat = data.frame(.Y = y[, i], y[, parent_mask, drop = FALSE])
    } else {
      dat = data.frame(.Y = y[, i])
    }
  } else {
    if (sum(parent_mask) > 0) {
      dat = data.frame(.Y = y[, i], y[, parent_mask, drop = FALSE], Z)
    } else {
      dat = data.frame(.Y = y[, i], Z)
    }
  }
  return(dat)
}

# --- Compute IC for node i given parent mask ---------------------------------
.node_ic = function(y, i, parent_mask, nq, ic, method, w, n_eff, Z) {
  has_parents = sum(parent_mask) > 0
  has_Z = !is.null(Z)
  if (!has_parents && !has_Z) {
    # Pure intercept-only model
    return(.null_ic(y[, i], nq[i], ic, method, w, n_eff, Z = NULL))
  }
  dat = .build_reg_data(y, i, parent_mask, Z)
  frm = .Y ~ .
  nq_x_parents = if (has_parents) nq[parent_mask] else integer(0)
  return(mypolr(frm, data = dat, ic = ic, method = method,
                nq_y = nq[i], nq_x = nq_x_parents, w = w, n_eff = n_eff))
}


# =============================================================================
# Core hill-climbing algorithm (full search)
# =============================================================================
oBN_greedy = function(y, gam = NULL, ic = "bic", method = "probit",
                      verbose = FALSE, maxit = 50,
                      w = NULL, n_eff = NULL, Z = NULL, bl_mat = NULL) {
  n = nrow(y)
  q = ncol(y)
  nq = rep(0, q)
  for (i in 1:q) nq[i] = nlevels(y[, i])
  
  if (is.null(gam)) {
    gam = matrix(FALSE, q, q)
  } else {
    gam = (gam != 0)
  }
  if (is.null(bl_mat)) bl_mat = matrix(FALSE, q, q)
  if (is.null(n_eff)) n_eff_use = n else n_eff_use = n_eff
  
  # Enforce blacklist on initial graph
  gam = gam & (!bl_mat)
  
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
  
  # Initialize IC for each node
  ic_best = rep(0, q)
  for (i in 1:q) {
    ic_best[i] = .node_ic(y, i, gam[i, ], nq, ic, method, w, n_eff_use, Z)
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
    
    # --- Add / Delete moves ---
    for (i in 1:q) {
      for (j in 1:(q - 1)) {
        if (admissible(i, ind_q[i, j], gam, bl_mat)) {
          if (gam[i, ind_q[i, j]]) {
            # Delete
            gam_new[i, ind_q[i, j]] = FALSE
            ic_best_new = .node_ic(y, i, gam_new[i, ], nq, ic, method,
                                   w, n_eff_use, Z)
            ic_improv_new = ic_best[i] - ic_best_new
            if (ic_improv_new > ic_improv) {
              ic_improv = ic_improv_new
              act_ind = c(i, ind_q[i, j])
              state = "del"
            }
            gam_new[i, ind_q[i, j]] = TRUE
          } else {
            # Add
            gam_new[i, ind_q[i, j]] = TRUE
            ic_best_new = .node_ic(y, i, gam_new[i, ], nq, ic, method,
                                   w, n_eff_use, Z)
            ic_improv_new = ic_best[i] - ic_best_new
            if (ic_improv_new > ic_improv) {
              ic_improv = ic_improv_new
              act_ind = c(i, ind_q[i, j])
              state = "add"
            }
            gam_new[i, ind_q[i, j]] = FALSE
          }
        }
      }
    }
    
    # --- Reverse moves ---
    for (i in 1:q) {
      for (j in 1:(q - 1)) {
        if (admissible_rev(i, ind_q[i, j], gam, bl_mat)) {
          tmp = gam_new[i, ind_q[i, j]]
          gam_new[ind_q[i, j], i] = tmp
          gam_new[i, ind_q[i, j]] = !tmp
          
          ic_rev_best_new[1] = .node_ic(y, i, gam_new[i, ], nq, ic, method,
                                        w, n_eff_use, Z)
          ic_rev_best_new[2] = .node_ic(y, ind_q[i, j], gam_new[ind_q[i, j], ],
                                        nq, ic, method, w, n_eff_use, Z)
          
          ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
          ic_improv_rev_new[2] = ic_best[ind_q[i, j]] - ic_rev_best_new[2]
          ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
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
    
    # --- Apply best move ---
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
      message(paste0(iter, " iterations completed | ", ic, " = ", sum(ic_best)))
    }
  }
  
  if (iter == maxit) {
    warning("Maximum iterations reached. Algorithm may not have converged.")
  }
  return(list(gam = gam + 0, ic_best = sum(ic_best)))
}


# =============================================================================
# Hill-climbing restricted to CPDAG edge reversals
# =============================================================================
oBN_greedy_CPDAG = function(y, gam = NULL, ic = "bic", edge_list = NULL,
                            method = "probit", verbose = FALSE, maxit = 50,
                            w = NULL, n_eff = NULL, Z = NULL, bl_mat = NULL) {
  n = nrow(y)
  q = ncol(y)
  nq = rep(0, q)
  for (i in 1:q) nq[i] = nlevels(y[, i])
  
  if (is.null(gam)) gam = matrix(FALSE, q, q) else gam = (gam != 0)
  if (is.null(bl_mat)) bl_mat = matrix(FALSE, q, q)
  if (is.null(n_eff)) n_eff_use = n else n_eff_use = n_eff
  
  gam = gam & (!bl_mat)
  
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
    ic_best[i] = .node_ic(y, i, gam[i, ], nq, ic, method, w, n_eff_use, Z)
  }
  
  while (ic_improv > 0 && iter < maxit) {
    iter = iter + 1
    ic_improv = -Inf
    ic_improv_rev = rep(-Inf, 2)
    gam_new = gam
    ic_improv_rev_new = rep(-Inf, 2)
    ic_rev_best_new = rep(-Inf, 2)
    
    for (i in 1:q) {
      if (nind_q[i] > 0) {
        for (j in 1:(nind_q[i])) {
          if (admissible_rev(i, ind_q[[i]][j], gam, bl_mat)) {
            tmp = gam_new[i, ind_q[[i]][j]]
            gam_new[ind_q[[i]][j], i] = tmp
            gam_new[i, ind_q[[i]][j]] = !tmp
            
            ic_rev_best_new[1] = .node_ic(y, i, gam_new[i, ], nq, ic, method,
                                          w, n_eff_use, Z)
            ic_rev_best_new[2] = .node_ic(y, ind_q[[i]][j],
                                          gam_new[ind_q[[i]][j], ], nq, ic,
                                          method, w, n_eff_use, Z)
            
            ic_improv_rev_new[1] = ic_best[i] - ic_rev_best_new[1]
            ic_improv_rev_new[2] = ic_best[ind_q[[i]][j]] - ic_rev_best_new[2]
            ic_improv_new = ic_improv_rev_new[1] + ic_improv_rev_new[2]
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
      message(paste0(iter, " iterations completed | ", ic, " = ", sum(ic_best)))
    }
  }
  
  if (iter == maxit) {
    warning("Maximum iterations reached. Algorithm may not have converged.")
  }
  return(list(gam = gam + 0, ic_best = sum(ic_best)))
}


# =============================================================================
# Multiple random restarts wrapper
# =============================================================================
oBN_greedy_wrap = function(y, ic = "bic", edge_list = NULL, method = "probit",
                           nstart = 1, verbose = FALSE, maxit = 50,
                           w = NULL, n_eff = NULL, Z = NULL, bl_mat = NULL) {
  q = ncol(y)
  gam_list = vector("list", nstart)
  ic_best_list = rep(NA, nstart)
  
  if (nstart == 1) {
    if (is.null(edge_list)) {
      fit = oBN_greedy(y, gam = NULL, ic = ic, method = method,
                       verbose = verbose, maxit = maxit,
                       w = w, n_eff = n_eff, Z = Z, bl_mat = bl_mat)
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
      # Enforce blacklist on initial CPDAG orientation
      if (!is.null(bl_mat)) gam = gam * (!bl_mat)
      fit = oBN_greedy_CPDAG(y, gam = gam, ic = ic, edge_list = und_edge,
                             method = method, verbose = verbose, maxit = maxit,
                             w = w, n_eff = n_eff, Z = Z, bl_mat = bl_mat)
    }
    gam_list[[1]] = fit$gam
    ic_best_list[1] = fit$ic_best
  } else {
    if (is.null(edge_list)) {
      netlist = bnlearn::random.graph(
        nodes = as.character(1:q), method = "ordered",
        num = nstart - 1, prob = 1 / q
      )
      if (nstart == 2) netlist = list(netlist)
      for (i in 1:nstart) {
        gam = matrix(FALSE, q, q)
        if (i != 1) {
          gam[apply(netlist[[i - 1]]$arcs, 2, as.numeric)] = TRUE
        }
        # Enforce blacklist on random initial graph
        if (!is.null(bl_mat)) gam = gam & (!bl_mat)
        fit = oBN_greedy(y, gam = gam, ic = ic, method = method,
                         verbose = verbose, maxit = maxit,
                         w = w, n_eff = n_eff, Z = Z, bl_mat = bl_mat)
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
        if (!is.null(bl_mat)) gam_ini = gam_ini * (!bl_mat)
        fit = oBN_greedy_CPDAG(y, gam = gam_ini, ic = ic, edge_list = und_edge,
                               method = method, verbose = verbose, maxit = maxit,
                               w = w, n_eff = n_eff, Z = Z, bl_mat = bl_mat)
        gam_list[[i]] = fit$gam
        ic_best_list[i] = fit$ic_best
      }
    }
  }
  
  i = which.min(ic_best_list)
  return(list(gam = gam_list[[i]], ic_best = ic_best_list[i]))
}


# =============================================================================
# Exhaustive search (q = 2 or 3)
# =============================================================================
oBN_exhaust = function(y, gam_list = NULL, ic = "bic", method = "probit",
                       w = NULL, n_eff = NULL, Z = NULL, bl_mat = NULL) {
  n = nrow(y)
  q = ncol(y)
  if (is.null(bl_mat)) bl_mat = matrix(FALSE, q, q)
  if (is.null(n_eff)) n_eff_use = n else n_eff_use = n_eff
  
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
      stop("Exhaustive search requires q = 2 or 3.")
    }
  }
  
  # Filter out DAGs that violate the blacklist
  nq = rep(0, q)
  for (i in 1:q) nq[i] = nlevels(y[, i])
  
  n_graphs = dim(gam_list)[3]
  IC = rep(Inf, n_graphs)
  
  for (m in 1:n_graphs) {
    gam = gam_list[, , m]
    # Skip if any blacklisted edge is present
    if (any(gam & bl_mat)) next
    IC[m] = 0
    for (i in 1:q) {
      parent_mask = as.logical(gam[i, ])
      IC[m] = IC[m] + .node_ic(y, i, parent_mask, nq, ic, method,
                               w, n_eff_use, Z)
    }
  }
  
  mi = which.min(IC)
  return(list(gam = gam_list[, , mi] + 0, ic_best = IC[mi]))
}


# =============================================================================
# Internal single-run interface
# =============================================================================
OCD = function(y, search = "greedy", ic = "bic", edge_list = NULL,
               link = "probit", G = NULL, nstart = 1,
               verbose = FALSE, maxit = 50,
               w = NULL, n_eff = NULL, Z = NULL, bl_mat = NULL) {
  if (search == "exhaust") {
    G = oBN_exhaust(y, G, ic, link, w = w, n_eff = n_eff, Z = Z,
                    bl_mat = bl_mat)
  } else {
    G = oBN_greedy_wrap(y, ic, edge_list, link, nstart, verbose, maxit,
                        w = w, n_eff = n_eff, Z = Z, bl_mat = bl_mat)
  }
  return(G)
}


# =============================================================================
# Main exported function
# =============================================================================
#' @title Causal Discovery for Ordinal Categorical Data (Extended)
#'
#' @description Estimate a causal DAG for ordinal categorical data (including
#'   binary variables) with greedy or exhaustive search. Supports survey
#'   weights, confounders, blacklisted edges, and bootstrapping.
#'
#' @param y a data.frame with each column being a factor (ordinal, >=2 levels).
#' @param search "greedy" (default) or "exhaust" (for q <= 3).
#' @param ic information criterion: "bic" (default) or "aic".
#' @param edge_list edge list of a CPDAG to restrict search space.
#' @param link link function for ordinal regression: "probit" (default),
#'   "logistic", "loglog", "cloglog", "cauchit".
#' @param G list of DAG adjacency matrices for exhaustive search.
#' @param nstart number of random restarts for greedy search.
#' @param verbose if TRUE, print progress messages.
#' @param maxit maximum iterations for greedy search (default 50).
#' @param boot number of bootstrap samples (default NULL = no bootstrap).
#' @param weights optional numeric vector of survey/sampling weights (length n).
#'   Affects both the weighted likelihood and the BIC penalty via Kish's
#'   effective sample size: n_eff = (sum(w))^2 / sum(w^2).
#' @param Z optional data.frame of confounders (numeric and/or factor columns).
#'   These are included as covariates in every regression but are NOT part of
#'   the DAG — no edges to/from Z are learned.
#' @param blacklist optional two-column data.frame or matrix with columns
#'   "from" and "to" specifying directed edges to forbid. Values can be
#'   column indices (integer) or column names of y. Only the specified
#'   direction is blocked.
#'
#' @return If boot is NULL: a list with gam (adjacency matrix) and ic_best.
#'   If boot > 0: a list of length boot, each element a list with gam and
#'   ic_best. The average adjacency matrix gives edge inclusion probabilities.
#'
#' @export
#'
#' @examples
#' set.seed(2020)
#' n = 1000; q = 3
#' y = u = matrix(0, n, q)
#' u[, 1] = 4 * rnorm(n)
#' y[, 1] = (u[, 1] > 1) + (u[, 1] > 2)
#' for (j in 2:q) {
#'   u[, j] = 2 * y[, j-1] + rnorm(n)
#'   y[, j] = (u[, j] > 1) + (u[, j] > 2)
#' }
#' A = matrix(0, q, q); A[2,1] = A[3,2] = 1
#' y = as.data.frame(y)
#' for (j in 1:q) y[, j] = as.factor(y[, j])
#'
#' # Basic usage
#' G = OrdCD(y)
#'
#' # With survey weights
#' w = runif(n, 0.5, 2)
#' G_w = OrdCD(y, weights = w)
#'
#' # With confounders
#' Z = data.frame(age = rnorm(n), gender = factor(sample(0:1, n, TRUE)))
#' G_z = OrdCD(y, Z = Z)
#'
#' # With blacklist (forbid edge 1->3)
#' bl = data.frame(from = 1, to = 3)
#' G_bl = OrdCD(y, blacklist = bl)
#'
#' # With bootstrapping
#' \dontrun{
#' G_boot = OrdCD(y, boot = 100, weights = w, Z = Z)
#' G_avg = Reduce("+", lapply(G_boot, function(x) x$gam)) / length(G_boot)
#' }
OrdCD = function(y,
                 search = "greedy",
                 ic = "bic",
                 edge_list = NULL,
                 link = "probit",
                 G = NULL,
                 nstart = 1,
                 verbose = FALSE,
                 maxit = 50,
                 boot = NULL,
                 weights = NULL,
                 Z = NULL,
                 blacklist = NULL) {
  
  # --- Input validation ------------------------------------------------------
  q = ncol(y)
  n = nrow(y)
  
  # Validate y: all columns must be factors with >= 2 levels
  for (j in 1:q) {
    if (!is.factor(y[, j]))
      stop(paste0("Column ", j, " of y must be a factor."))
    if (nlevels(y[, j]) < 2)
      stop(paste0("Column ", j, " of y must have at least 2 levels."))
  }
  
  # Validate weights
  w = NULL
  n_eff = NULL
  if (!is.null(weights)) {
    stopifnot(is.numeric(weights), length(weights) == n, all(weights > 0))
    w = weights
    n_eff = .eff_n(w)
  }
  
  # Validate Z
  if (!is.null(Z)) {
    Z = as.data.frame(Z)
    stopifnot(nrow(Z) == n)
    # Ensure column names don't clash with y
    if (any(names(Z) %in% names(y))) {
      names(Z) = paste0(".Z_", seq_len(ncol(Z)))
    }
    # If Z columns have no names, assign them
    if (is.null(names(Z))) names(Z) = paste0(".Z_", seq_len(ncol(Z)))
  }
  
  # Build blacklist matrix
  bl_mat = .make_blacklist_mat(blacklist, q, colnames(y))
  
  # --- Run -------------------------------------------------------------------
  if (is.null(boot)) {
    result = OCD(y, search, ic, edge_list, link, G, nstart, verbose, maxit,
                 w = w, n_eff = n_eff, Z = Z, bl_mat = bl_mat)
  } else {
    result = vector("list", boot)
    for (b in 1:boot) {
      idx = sample(n, replace = TRUE)
      y_b = y[idx, , drop = FALSE]
      # Re-level factors to drop empty levels in bootstrap sample
      for (j in 1:q) y_b[, j] = droplevels(y_b[, j])
      # Check all variables still have >= 2 levels after droplevels
      skip = FALSE
      for (j in 1:q) {
        if (nlevels(y_b[, j]) < 2) { skip = TRUE; break }
      }
      if (skip) { b = b - 1; next }  # Resample if degenerate
      
      w_b = if (!is.null(w)) w[idx] else NULL
      n_eff_b = if (!is.null(w_b)) .eff_n(w_b) else NULL
      Z_b = if (!is.null(Z)) Z[idx, , drop = FALSE] else NULL
      
      result[[b]] = OCD(y_b, search, ic, edge_list, link, NULL, nstart,
                        verbose, maxit,
                        w = w_b, n_eff = n_eff_b, Z = Z_b, bl_mat = bl_mat)
      if (verbose) message(paste0("Bootstrap sample ", b, "/", boot, " done."))
    }
  }
  return(result)
}