# ==============================================================================
# sim_utils.R
# Simulation utilities for swa-oBN cluster simulations
#
# Contents:
#   1. Data generation from a known ordinal Bayesian network
#   2. Survey sampling with informative designs
#   3. Structural Hamming Distance and related metrics
#   4. Single replication runner (Study 1)
#   5. Confounded population generation + replication runner (Study 2)
#   6. Study 3 helpers (PC, FCI, LiNGAM)
# ==============================================================================


# ==============================================================================
# SECTION 1: DATA GENERATION FROM A KNOWN oBN
# ==============================================================================

#' Generate a random DAG adjacency matrix
generate_random_dag <- function(q, prob = 0.3) {
  order <- sample(q)
  gam <- matrix(0L, q, q)
  for (idx in 2:q) {
    node <- order[idx]
    for (pidx in 1:(idx - 1)) {
      parent <- order[pidx]
      if (runif(1) < prob) {
        gam[node, parent] <- 1L
      }
    }
  }
  return(gam)
}


#' Generate ordinal data from a known oBN model (cumulative probit)
generate_obn_data <- function(n, gam, sigma = 1.0, L = 3, Z = NULL,
                              delta_strength = 0, seed = NULL) {
  q <- nrow(gam)
  if (length(L) == 1) L <- rep(L, q)
  
  g <- igraph::graph_from_adjacency_matrix(t(gam), mode = "directed")
  topo <- as.integer(igraph::topo_sort(g))
  
  if (!is.null(seed)) set.seed(seed)
  
  params <- vector("list", q)
  for (j in seq_len(q)) {
    pa_j <- which(gam[j, ] == 1)
    n_levels <- L[j]
    
    if (n_levels > 2) {
      probs <- seq(0, 1, length.out = n_levels + 1)[2:n_levels]
      gamma_j <- qnorm(probs)
    } else {
      gamma_j <- 0
    }
    
    alpha_j <- 0
    
    beta_j <- list()
    if (length(pa_j) > 0) {
      for (k in pa_j) {
        beta_j[[as.character(k)]] <- rnorm(1, 0, sigma)
      }
    }
    
    delta_j <- NULL
    if (!is.null(Z) && delta_strength > 0) {
      delta_j <- rnorm(ncol(Z), 0, delta_strength)
    }
    
    params[[j]] <- list(
      gamma = gamma_j, alpha = alpha_j,
      beta  = beta_j,  delta = delta_j,
      parents = pa_j
    )
  }
  
  X <- matrix(0L, n, q)
  for (j in topo) {
    par <- params[[j]]
    pa_j <- par$parents
    n_levels <- L[j]
    
    eta <- rep(par$alpha, n)
    if (length(pa_j) > 0) {
      for (k in pa_j) {
        eta <- eta + par$beta[[as.character(k)]] * X[, k]
      }
    }
    if (!is.null(Z) && !is.null(par$delta)) {
      eta <- eta + as.matrix(Z) %*% par$delta
    }
    
    u <- runif(n)
    for (i in seq_len(n)) {
      cum_probs <- pnorm(par$gamma - eta[i])
      cum_probs <- c(cum_probs, 1.0)
      X[i, j] <- findInterval(u[i], c(0, cum_probs))
      X[i, j] <- min(X[i, j], n_levels)
      X[i, j] <- max(X[i, j], 1L)
    }
  }
  
  df <- as.data.frame(X)
  colnames(df) <- paste0("X", seq_len(q))
  for (j in seq_len(q)) {
    df[, j] <- factor(df[, j], levels = seq_len(L[j]), ordered = TRUE)
  }
  
  return(list(
    data       = df,
    params     = params,
    gam        = gam,
    topo_order = topo,
    L          = L
  ))
}


# ==============================================================================
# SECTION 2: SURVEY SAMPLING WITH INFORMATIVE DESIGNS
# ==============================================================================

#' Apply informative survey sampling to a population
apply_survey_sampling <- function(pop_data, design = "moderate",
                                  oversample_vars = NULL,
                                  target_n = NULL) {
  N <- nrow(pop_data)
  
  if (design == "none") {
    if (is.null(target_n)) target_n <- round(N * 0.10)
    idx <- sample(N, target_n)
    return(list(
      sample_data   = pop_data[idx, , drop = FALSE],
      weights       = rep(N / target_n, target_n),
      pi            = rep(target_n / N, target_n),
      pop_indices   = idx,
      design_effect = 1.0,
      n_eff         = target_n
    ))
  }
  
  if (is.null(oversample_vars)) {
    oversample_vars <- c(1, ncol(pop_data))
  }
  
  target_score <- rowSums(sapply(oversample_vars, function(j) {
    as.numeric(pop_data[, j])
  }))
  
  threshold <- quantile(target_score, 0.75)
  is_target <- target_score >= threshold
  
  pi_base <- switch(design,
                    "mild"     = c(target = 0.30, other = 0.08),
                    "moderate" = c(target = 0.60, other = 0.05),
                    "extreme"  = c(target = 1.00, other = 0.03),
                    stop("Unknown design: ", design)
  )
  
  pi_i <- ifelse(is_target, pi_base["target"], pi_base["other"])
  
  if (!is.null(target_n)) {
    expected_n <- sum(pi_i)
    scale_factor <- target_n / expected_n
    pi_i <- pmin(pi_i * scale_factor, 1.0)
  }
  
  selected <- rbinom(N, 1, pi_i) == 1
  sample_data <- pop_data[selected, , drop = FALSE]
  pi_selected <- pi_i[selected]
  w_selected  <- 1 / pi_selected
  
  n_eff <- sum(w_selected)^2 / sum(w_selected^2)
  deff  <- sum(selected) / n_eff
  
  return(list(
    sample_data   = sample_data,
    weights       = w_selected,
    pi            = pi_selected,
    pop_indices   = which(selected),
    design_effect = deff,
    n_eff         = n_eff
  ))
}


# ==============================================================================
# SECTION 3: STRUCTURAL HAMMING DISTANCE AND METRICS
# ==============================================================================

compute_shd <- function(true_gam, est_gam) {
  true_gam <- (true_gam != 0) * 1L
  est_gam  <- (est_gam != 0) * 1L
  q <- nrow(true_gam)
  shd <- 0
  
  for (i in seq_len(q)) {
    for (j in (i + 1):min(q, q)) {
      if (j > q) next
      true_ij <- true_gam[i, j]; true_ji <- true_gam[j, i]
      est_ij  <- est_gam[i, j];  est_ji  <- est_gam[j, i]
      
      true_edge <- true_ij + true_ji
      est_edge  <- est_ij + est_ji
      
      if (true_edge == 0 && est_edge > 0) {
        shd <- shd + 1
      } else if (true_edge > 0 && est_edge == 0) {
        shd <- shd + 1
      } else if (true_edge > 0 && est_edge > 0) {
        if (true_ij != est_ij || true_ji != est_ji) {
          shd <- shd + 1
        }
      }
    }
  }
  return(shd)
}


compute_edge_metrics <- function(true_gam, est_gam) {
  true_gam <- (true_gam != 0) * 1L
  est_gam  <- (est_gam != 0) * 1L
  q <- nrow(true_gam)
  
  true_skel <- pmax(true_gam, t(true_gam))
  est_skel  <- pmax(est_gam, t(est_gam))
  ut <- upper.tri(true_gam)
  
  true_edges <- true_skel[ut]
  est_edges  <- est_skel[ut]
  
  tp_skel <- sum(true_edges == 1 & est_edges == 1)
  fp_skel <- sum(true_edges == 0 & est_edges == 1)
  fn_skel <- sum(true_edges == 1 & est_edges == 0)
  tn_skel <- sum(true_edges == 0 & est_edges == 0)
  
  tpr_skel <- ifelse(tp_skel + fn_skel > 0, tp_skel / (tp_skel + fn_skel), NA)
  fpr_skel <- ifelse(fp_skel + tn_skel > 0, fp_skel / (fp_skel + tn_skel), NA)
  precision_skel <- ifelse(tp_skel + fp_skel > 0, tp_skel / (tp_skel + fp_skel), NA)
  
  n_correct_skel <- 0
  n_correct_orient <- 0
  for (i in seq_len(q)) {
    for (j in (i + 1):min(q, q)) {
      if (j > q) next
      if (true_skel[i, j] == 1 && est_skel[i, j] == 1) {
        n_correct_skel <- n_correct_skel + 1
        if (true_gam[i, j] == est_gam[i, j] && true_gam[j, i] == est_gam[j, i]) {
          n_correct_orient <- n_correct_orient + 1
        }
      }
    }
  }
  
  orient_acc <- ifelse(n_correct_skel > 0, n_correct_orient / n_correct_skel, NA)
  n_total_true <- sum(true_skel[ut])
  correct_orient_rate <- ifelse(n_total_true > 0, n_correct_orient / n_total_true, NA)
  
  return(list(
    shd              = compute_shd(true_gam, est_gam),
    tpr_skeleton     = tpr_skel,
    fpr_skeleton     = fpr_skel,
    precision_skeleton = precision_skel,
    orientation_acc  = orient_acc,
    correct_orient_rate = correct_orient_rate,
    n_true_edges     = sum(true_skel[ut]),
    n_est_edges      = sum(est_skel[ut]),
    tp_skeleton      = tp_skel,
    fp_skeleton      = fp_skel,
    fn_skeleton      = fn_skel
  ))
}


# ==============================================================================
# SECTION 4: SINGLE REPLICATION RUNNER (Study 1)
# ==============================================================================

run_single_replication <- function(gam, N = 20000, sigma = 1.0, L = 3,
                                   design = "moderate", target_n = 1000,
                                   nstart = 10, maxit = 100,
                                   Z_pop = NULL, delta_strength = 0,
                                   oversample_vars = NULL,
                                   swa_obn_fn = NULL) {
  q <- nrow(gam)
  
  pop <- generate_obn_data(
    n = N, gam = gam, sigma = sigma, L = L,
    Z = Z_pop, delta_strength = delta_strength
  )
  
  if (is.null(oversample_vars)) oversample_vars <- c(1, q)
  
  svy <- apply_survey_sampling(
    pop_data = pop$data, design = design,
    oversample_vars = oversample_vars, target_n = target_n
  )
  
  y_sample <- svy$sample_data
  w_sample <- svy$weights
  n_actual <- nrow(y_sample)
  
  Z_sample <- NULL
  if (!is.null(Z_pop)) {
    Z_sample <- as.data.frame(Z_pop[svy$pop_indices, , drop = FALSE])
  }
  
  fit_naive <- tryCatch({
    swa_obn_fn(
      y = y_sample, Z = Z_sample, weights = NULL,
      search = "greedy", ic = "bic", link = "probit",
      nstart = nstart, boot = NULL, verbose = FALSE, maxit = maxit
    )
  }, error = function(e) list(gam = matrix(0, q, q)))
  
  fit_svy <- tryCatch({
    swa_obn_fn(
      y = y_sample, Z = Z_sample, weights = w_sample,
      search = "greedy", ic = "bic", link = "probit",
      nstart = nstart, boot = NULL, verbose = FALSE, maxit = maxit
    )
  }, error = function(e) list(gam = matrix(0, q, q)))
  
  metrics_naive <- compute_edge_metrics(gam, fit_naive$gam)
  metrics_svy   <- compute_edge_metrics(gam, fit_svy$gam)
  
  return(data.frame(
    n_actual     = n_actual,
    n_eff        = svy$n_eff,
    deff         = svy$design_effect,
    
    shd_naive    = metrics_naive$shd,
    tpr_naive    = metrics_naive$tpr_skeleton,
    fpr_naive    = metrics_naive$fpr_skeleton,
    prec_naive   = metrics_naive$precision_skeleton,
    orient_naive = metrics_naive$orientation_acc,
    cor_naive    = metrics_naive$correct_orient_rate,
    nedge_naive  = metrics_naive$n_est_edges,
    
    shd_svy      = metrics_svy$shd,
    tpr_svy      = metrics_svy$tpr_skeleton,
    fpr_svy      = metrics_svy$fpr_skeleton,
    prec_svy     = metrics_svy$precision_skeleton,
    orient_svy   = metrics_svy$orientation_acc,
    cor_svy      = metrics_svy$correct_orient_rate,
    nedge_svy    = metrics_svy$n_est_edges,
    
    n_true_edges = metrics_naive$n_true_edges,
    stringsAsFactors = FALSE
  ))
}


# ==============================================================================
# SECTION 5: CONFOUNDED POPULATION + REPLICATION (Study 2)
# ==============================================================================

generate_confounded_population <- function(N, gam, sigma, L,
                                           confounded_nodes, delta_strength) {
  q <- nrow(gam)
  Z_pop <- matrix(rnorm(N, 0, 1), ncol = 1)
  colnames(Z_pop) <- "Z"
  
  # Generate parameters (with Z effects for all nodes initially)
  pop <- generate_obn_data(
    n = N, gam = gam, sigma = sigma, L = L,
    Z = Z_pop, delta_strength = delta_strength
  )
  
  # Regenerate data: apply Z effect only to confounded nodes
  topo <- pop$topo_order
  X <- matrix(0L, N, q)
  
  for (j in topo) {
    par <- pop$params[[j]]
    pa_j <- par$parents
    
    eta <- rep(par$alpha, N)
    if (length(pa_j) > 0) {
      for (k in pa_j) {
        eta <- eta + par$beta[[as.character(k)]] * X[, k]
      }
    }
    if (j %in% confounded_nodes && !is.null(par$delta)) {
      eta <- eta + as.numeric(Z_pop %*% par$delta)
    }
    
    u <- runif(N)
    for (i in seq_len(N)) {
      cum_probs <- pnorm(par$gamma - eta[i])
      cum_probs <- c(cum_probs, 1.0)
      X[i, j] <- findInterval(u[i], c(0, cum_probs))
      X[i, j] <- min(X[i, j], L)
      X[i, j] <- max(X[i, j], 1L)
    }
  }
  
  df <- as.data.frame(X)
  colnames(df) <- paste0("X", seq_len(q))
  for (j in seq_len(q)) {
    df[, j] <- factor(df[, j], levels = seq_len(L), ordered = TRUE)
  }
  
  return(list(data = df, Z = Z_pop, params = pop$params, gam = gam))
}


run_study2_replication <- function(gam, N, sigma, L, design, target_n,
                                   confounded_nodes, delta_strength,
                                   nstart, maxit, swa_obn_fn = NULL) {
  q <- nrow(gam)
  
  pop <- generate_confounded_population(
    N = N, gam = gam, sigma = sigma, L = L,
    confounded_nodes = confounded_nodes,
    delta_strength = delta_strength
  )
  
  svy <- apply_survey_sampling(
    pop_data = pop$data, design = design,
    oversample_vars = c(1, q), target_n = target_n
  )
  
  y_sample <- svy$sample_data
  w_sample <- svy$weights
  Z_sample <- as.data.frame(pop$Z[svy$pop_indices, , drop = FALSE])
  colnames(Z_sample) <- "Z"
  
  # (1) Naive: no weights, no Z
  fit_naive <- tryCatch(
    swa_obn_fn(y = y_sample, Z = NULL, weights = NULL,
               search = "greedy", ic = "bic", link = "probit",
               nstart = nstart, verbose = FALSE, maxit = maxit),
    error = function(e) list(gam = matrix(0, q, q))
  )
  
  # (2) Z-only
  fit_z <- tryCatch(
    swa_obn_fn(y = y_sample, Z = Z_sample, weights = NULL,
               search = "greedy", ic = "bic", link = "probit",
               nstart = nstart, verbose = FALSE, maxit = maxit),
    error = function(e) list(gam = matrix(0, q, q))
  )
  
  # (3) Weights-only
  fit_w <- tryCatch(
    swa_obn_fn(y = y_sample, Z = NULL, weights = w_sample,
               search = "greedy", ic = "bic", link = "probit",
               nstart = nstart, verbose = FALSE, maxit = maxit),
    error = function(e) list(gam = matrix(0, q, q))
  )
  
  # (4) Full swa-oBN
  fit_full <- tryCatch(
    swa_obn_fn(y = y_sample, Z = Z_sample, weights = w_sample,
               search = "greedy", ic = "bic", link = "probit",
               nstart = nstart, verbose = FALSE, maxit = maxit),
    error = function(e) list(gam = matrix(0, q, q))
  )
  
  m_naive <- compute_edge_metrics(gam, fit_naive$gam)
  m_z     <- compute_edge_metrics(gam, fit_z$gam)
  m_w     <- compute_edge_metrics(gam, fit_w$gam)
  m_full  <- compute_edge_metrics(gam, fit_full$gam)
  
  return(data.frame(
    n_actual = nrow(y_sample),
    n_eff    = svy$n_eff,
    deff     = svy$design_effect,
    
    shd_naive = m_naive$shd,   shd_z = m_z$shd,
    shd_w = m_w$shd,           shd_full = m_full$shd,
    
    fpr_naive = m_naive$fpr_skeleton,  fpr_z = m_z$fpr_skeleton,
    fpr_w = m_w$fpr_skeleton,          fpr_full = m_full$fpr_skeleton,
    
    tpr_naive = m_naive$tpr_skeleton,  tpr_z = m_z$tpr_skeleton,
    tpr_w = m_w$tpr_skeleton,          tpr_full = m_full$tpr_skeleton,
    
    orient_naive = m_naive$orientation_acc,  orient_z = m_z$orientation_acc,
    orient_w = m_w$orientation_acc,          orient_full = m_full$orientation_acc,
    
    nedge_naive = m_naive$n_est_edges, nedge_z = m_z$n_est_edges,
    nedge_w = m_w$n_est_edges,         nedge_full = m_full$n_est_edges,
    
    n_true_edges = m_naive$n_true_edges,
    stringsAsFactors = FALSE
  ))
}


# ==============================================================================
# SECTION 6: STUDY 3 HELPERS (PC, FCI, LiNGAM)
# ==============================================================================

# Design parameters for Study 3 sampling
DESIGN_PARAMS_S3 <- list(
  none     = c(target = 0.10, other = 0.10),
  mild     = c(target = 0.30, other = 0.08),
  moderate = c(target = 0.60, other = 0.05),
  extreme  = c(target = 1.00, other = 0.03)
)


#' Survey sampling for 3-variable continuous data (Study 3)
sample_with_design <- function(pop_data, design_name) {
  N <- nrow(pop_data)
  target_score <- scale(as.numeric(pop_data[, 1])) +
    scale(as.numeric(pop_data[, 3]))
  is_target <- target_score > quantile(target_score, 0.75)
  
  params <- DESIGN_PARAMS_S3[[design_name]]
  pi_i <- ifelse(is_target, params["target"], params["other"])
  
  selected <- rbinom(N, 1, pi_i) == 1
  w <- 1 / pi_i[selected]
  n_eff <- sum(w)^2 / sum(w^2)
  
  list(
    data    = pop_data[selected, , drop = FALSE],
    weights = w,
    n       = sum(selected),
    n_eff   = n_eff,
    deff    = sum(selected) / n_eff
  )
}


#' Survey-weighted CI test for PC/FCI
svy_ci_test <- function(x, y, S, suffStat) {
  dat <- suffStat$data
  wts <- suffStat$weights
  
  dat$.wt <- wts
  des <- survey::svydesign(ids = ~1, weights = ~.wt, data = dat)
  
  var_x <- colnames(dat)[x]
  var_y <- colnames(dat)[y]
  
  if (length(S) == 0) {
    form <- as.formula(paste(var_x, "~", var_y))
  } else {
    var_S <- colnames(dat)[S]
    form <- as.formula(paste(var_x, "~", var_y, "+",
                             paste(var_S, collapse = " + ")))
  }
  
  fit <- try(survey::svyglm(form, design = des), silent = TRUE)
  if (inherits(fit, "try-error")) return(1)
  
  coefs <- summary(fit)$coefficients
  if (var_y %in% rownames(coefs)) {
    return(coefs[var_y, "Pr(>|t|)"])
  } else {
    return(1)
  }
}


#' Survey-weighted LiNGAM via pseudo-population bootstrap
svy_lingam <- function(X_mat, weights, B = 50) {
  n <- nrow(X_mat)
  p <- ncol(X_mat)
  edge_freq <- matrix(0, p, p)
  prob_wts <- weights / sum(weights)
  ok <- 0
  
  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE, prob = prob_wts)
    pseudo_X <- X_mat[idx, ]
    
    res <- try(pcalg::lingam(pseudo_X, verbose = FALSE), silent = TRUE)
    if (!inherits(res, "try-error")) {
      adj <- (t(res$Bpruned) != 0) * 1L
      edge_freq <- edge_freq + adj
      ok <- ok + 1
    }
  }
  
  if (ok == 0) return(matrix(0, p, p))
  
  consensus <- (edge_freq / ok > 0.5) * 1L
  rownames(consensus) <- colnames(consensus) <- colnames(X_mat)
  return(consensus)
}


#' Skeleton metrics for 3-node chain (true: X-Y, Y-Z, no X-Z)
skeleton_metrics_3node <- function(adj_mat) {
  skel <- pmax(adj_mat, t(adj_mat))
  skel <- (skel != 0) * 1L
  
  has_12 <- skel[1, 2] == 1
  has_23 <- skel[2, 3] == 1
  has_13 <- skel[1, 3] == 1
  
  shd <- 0
  if (!has_12) shd <- shd + 1
  if (!has_23) shd <- shd + 1
  if (has_13)  shd <- shd + 1
  
  list(
    shd         = shd,
    has_12      = has_12,
    has_23      = has_23,
    spurious_13 = has_13,
    correct     = (has_12 & has_23 & !has_13)
  )
}


#' Run one replication of Study 3
#'
#' @param algorithm "PC", "FCI", or "LiNGAM"
#' @param design    "none", "mild", "moderate", "extreme"
#' @param N_POP     Population size
#' @return data.frame with one row of metrics, or NULL on failure
run_study3_replication <- function(algorithm, design, N_POP = 20000) {
  
  # --- Generate population: chain X -> Y -> Z ---
  if (algorithm == "PC") {
    X <- rnorm(N_POP)
    Y <- 1.5 * X + rnorm(N_POP)
    Z <- 1.5 * Y + rnorm(N_POP)
  } else if (algorithm == "FCI") {
    X <- rnorm(N_POP)
    Y <- 0.8 * X + rnorm(N_POP)
    Z <- 0.8 * Y + rnorm(N_POP)
  } else {
    # LiNGAM: non-Gaussian noise for identifiability
    X <- runif(N_POP, -2, 2)
    Y <- 1.5 * X + runif(N_POP, -2, 2)
    Z <- 1.5 * Y + runif(N_POP, -2, 2)
  }
  pop <- data.frame(X = X, Y = Y, Z = Z)
  
  # --- Survey sampling ---
  svy <- sample_with_design(pop, design)
  sdat <- svy$data
  
  # --- Fit naive and survey-weighted versions ---
  if (algorithm == "PC") {
    suffStat_naive <- list(C = cor(sdat), n = nrow(sdat))
    pc_naive <- tryCatch(
      pcalg::pc(suffStat = suffStat_naive, indepTest = pcalg::gaussCItest,
                alpha = 0.05, labels = colnames(sdat), verbose = FALSE),
      error = function(e) NULL
    )
    suffStat_svy <- list(data = sdat, weights = svy$weights)
    pc_svy <- tryCatch(
      pcalg::pc(suffStat = suffStat_svy, indepTest = svy_ci_test,
                alpha = 0.05, labels = colnames(sdat), verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(pc_naive) || is.null(pc_svy)) return(NULL)
    adj_naive <- as(pc_naive@graph, "matrix")
    adj_svy   <- as(pc_svy@graph, "matrix")
    
  } else if (algorithm == "FCI") {
    suffStat_naive <- list(C = cor(sdat), n = nrow(sdat))
    fci_naive <- tryCatch(
      pcalg::fci(suffStat = suffStat_naive, indepTest = pcalg::gaussCItest,
                 alpha = 0.05, labels = colnames(sdat), p = 3, verbose = FALSE),
      error = function(e) NULL
    )
    suffStat_svy <- list(data = sdat, weights = svy$weights)
    fci_svy <- tryCatch(
      pcalg::fci(suffStat = suffStat_svy, indepTest = svy_ci_test,
                 alpha = 0.05, labels = colnames(sdat), p = 3, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(fci_naive) || is.null(fci_svy)) return(NULL)
    adj_naive <- (fci_naive@amat != 0) * 1L
    adj_svy   <- (fci_svy@amat != 0) * 1L
    
  } else {
    # LiNGAM
    lingam_naive <- tryCatch({
      res <- pcalg::lingam(as.matrix(sdat), verbose = FALSE)
      adj <- (t(res$Bpruned) != 0) * 1L
      rownames(adj) <- colnames(adj) <- colnames(sdat)
      adj
    }, error = function(e) NULL)
    lingam_svy <- tryCatch(
      svy_lingam(as.matrix(sdat), svy$weights, B = 50),
      error = function(e) NULL
    )
    if (is.null(lingam_naive) || is.null(lingam_svy)) return(NULL)
    adj_naive <- lingam_naive
    adj_svy   <- lingam_svy
  }
  
  m_naive <- skeleton_metrics_3node(adj_naive)
  m_svy   <- skeleton_metrics_3node(adj_svy)
  
  data.frame(
    n = svy$n, n_eff = svy$n_eff, deff = svy$deff,
    shd_naive = m_naive$shd, spurious_naive = m_naive$spurious_13,
    correct_naive = m_naive$correct,
    shd_svy = m_svy$shd, spurious_svy = m_svy$spurious_13,
    correct_svy = m_svy$correct,
    stringsAsFactors = FALSE
  )
}