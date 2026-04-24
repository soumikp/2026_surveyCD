# ==============================================================================
# sim_utils.R
# Simulation utilities for swa-oBN paper
#
# Contents:
#   1. Data generation from a known ordinal Bayesian network
#   2. Survey sampling with informative designs
#   3. Structural Hamming Distance and related metrics
#   4. Helper functions
# ==============================================================================

# ==============================================================================
# SECTION 1: DATA GENERATION FROM A KNOWN oBN
# ==============================================================================

#' Generate a random DAG adjacency matrix
#'
#' @param q Number of nodes
#' @param prob Edge probability (Erdos-Renyi on the upper triangle of a
#'        topological ordering)
#' @return A q x q adjacency matrix where gam[i, j] = 1 means j -> i
generate_random_dag <- function(q, prob = 0.3) {
  # Generate a random permutation as the topological order
  order <- sample(q)
  gam <- matrix(0L, q, q)
  
  for (idx in 2:q) {
    node <- order[idx]
    for (pidx in 1:(idx - 1)) {
      parent <- order[pidx]
      if (runif(1) < prob) {
        gam[node, parent] <- 1L  # parent -> node
      }
    }
  }
  return(gam)
}


#' Generate ordinal data from a known oBN model (cumulative probit)
#'
#' Data-generating process:
#'   For each node j in topological order:
#'     P(X_j <= l | X_pa(j), Z) = Phi(gamma_jl - sum_k beta_jk * X_k - delta_j' Z - alpha_j)
#'
#' @param n       Number of observations
#' @param gam     Adjacency matrix (q x q), gam[i,j]=1 means j -> i
#' @param sigma   Signal strength: beta_jk ~ N(0, sigma^2)
#' @param L       Number of ordinal levels per variable (scalar or vector of length q)
#' @param Z       Optional covariate matrix (n x p_z). If NULL, no covariates.
#' @param delta_strength  Strength of covariate effects: delta ~ N(0, delta_strength^2)
#' @param seed    Random seed for parameter generation (NULL = no seed for params)
#'
#' @return List with:
#'   $data      - data.frame of ordered factors (n x q)
#'   $params    - list of parameter vectors per node
#'   $gam       - the input DAG
#'   $topo_order - topological ordering used
generate_obn_data <- function(n, gam, sigma = 1.0, L = 3, Z = NULL,
                              delta_strength = 0, seed = NULL) {
  q <- nrow(gam)
  
  if (length(L) == 1) L <- rep(L, q)
  
  # --- Topological sort ---
  g <- igraph::graph_from_adjacency_matrix(t(gam), mode = "directed")
  topo <- as.integer(igraph::topo_sort(g))
  
  # --- Generate parameters ---
  if (!is.null(seed)) set.seed(seed)
  
  params <- vector("list", q)
  
  for (j in seq_len(q)) {
    pa_j <- which(gam[j, ] == 1)
    n_levels <- L[j]
    
    # Thresholds: equally spaced on the probit scale
    # gamma_1 = 0 (identification), gamma_{L} = Inf (implicit)
    if (n_levels > 2) {
      # Place thresholds at quantiles of standard normal
      probs <- seq(0, 1, length.out = n_levels + 1)[2:n_levels]
      gamma_j <- qnorm(probs)
    } else {
      gamma_j <- 0  # single threshold for binary
    }
    
    # Intercept
    alpha_j <- 0
    
    # Parent effects: for each parent k, generate (L_k - 1) effect parameters
    # beta_jk1 = 0 (identification), beta_jk2, ..., beta_jkL_k
    beta_j <- list()
    if (length(pa_j) > 0) {
      for (k in pa_j) {
        if (L[k] > 2) {
          # For ordinal parents: use a single slope (proportional odds)
          # This is the standard oBN parameterization where X_k enters linearly
          beta_j[[as.character(k)]] <- rnorm(1, 0, sigma)
        } else {
          beta_j[[as.character(k)]] <- rnorm(1, 0, sigma)
        }
      }
    }
    
    # Covariate effects
    delta_j <- NULL
    if (!is.null(Z) && delta_strength > 0) {
      delta_j <- rnorm(ncol(Z), 0, delta_strength)
    }
    
    params[[j]] <- list(
      gamma = gamma_j,
      alpha = alpha_j,
      beta  = beta_j,
      delta = delta_j,
      parents = pa_j
    )
  }
  
  # --- Generate data in topological order ---
  X <- matrix(0L, n, q)
  
  for (j in topo) {
    par <- params[[j]]
    pa_j <- par$parents
    n_levels <- L[j]
    
    # Compute linear predictor for each observation
    eta <- rep(par$alpha, n)
    
    if (length(pa_j) > 0) {
      for (k in pa_j) {
        # oBN uses the numeric coding of the parent directly
        eta <- eta + par$beta[[as.character(k)]] * X[, k]
      }
    }
    
    if (!is.null(Z) && !is.null(par$delta)) {
      eta <- eta + as.matrix(Z) %*% par$delta
    }
    
    # Generate from cumulative probit model
    # P(X_j <= l) = Phi(gamma_l - eta)
    u <- runif(n)
    for (i in seq_len(n)) {
      cum_probs <- pnorm(par$gamma - eta[i])
      cum_probs <- c(cum_probs, 1.0)  # gamma_L = Inf => P(X<=L) = 1
      
      # Find which level: smallest l such that u <= P(X <= l)
      X[i, j] <- findInterval(u[i], c(0, cum_probs)) 
      X[i, j] <- min(X[i, j], n_levels)
      X[i, j] <- max(X[i, j], 1L)
    }
  }
  
  # --- Convert to data frame of ordered factors ---
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
#'
#' Oversampling probability depends on observed variable values.
#' This creates the design-induced bias that motivates survey weighting.
#'
#' @param pop_data    Data frame of the full population
#' @param design      Character: "mild", "moderate", "extreme", or "none"
#' @param oversample_vars Column indices of variables used to define oversampling
#' @param target_n    Approximate target sample size (actual may vary due to
#'                    Bernoulli sampling)
#'
#' @return List with:
#'   $sample_data   - data frame of sampled observations
#'   $weights       - survey weights (1/pi_i)
#'   $pi            - inclusion probabilities
#'   $pop_indices   - which population rows were sampled
#'   $design_effect - actual n / n_eff
apply_survey_sampling <- function(pop_data, design = "moderate",
                                  oversample_vars = NULL,
                                  target_n = NULL) {
  N <- nrow(pop_data)
  
  if (design == "none") {
    # Simple random sampling
    if (is.null(target_n)) target_n <- round(N * 0.10)
    idx <- sample(N, target_n)
    return(list(
      sample_data   = pop_data[idx, , drop = FALSE],
      weights       = rep(N / target_n, target_n),
      pi            = rep(target_n / N, target_n),
      pop_indices   = idx,
      design_effect = 1.0,
      n_eff         = target_n  # all weights equal => n_eff = n
    ))
  }
  
  # --- Define oversampling based on variable values ---
  if (is.null(oversample_vars)) {
    # Default: use first and last variables
    oversample_vars <- c(1, ncol(pop_data))
  }
  
  # Compute a "target score": sum of numeric levels of oversample vars
  # High score = oversampled
  target_score <- rowSums(sapply(oversample_vars, function(j) {
    as.numeric(pop_data[, j])
  }))
  
  # Threshold for "target group": top quantile
  threshold <- quantile(target_score, 0.75)
  is_target <- target_score >= threshold
  
  # Set inclusion probabilities by design strength
  pi_base <- switch(design,
                    "mild"     = c(target = 0.30, other = 0.08),
                    "moderate" = c(target = 0.60, other = 0.05),
                    "extreme"  = c(target = 1.00, other = 0.03),
                    stop("Unknown design: ", design)
  )
  
  pi_i <- ifelse(is_target, pi_base["target"], pi_base["other"])
  
  # Scale to approximate target_n if provided
  if (!is.null(target_n)) {
    expected_n <- sum(pi_i)
    scale_factor <- target_n / expected_n
    pi_i <- pmin(pi_i * scale_factor, 1.0)  # cap at 1
  }
  
  # --- Bernoulli sampling ---
  selected <- rbinom(N, 1, pi_i) == 1
  
  sample_data <- pop_data[selected, , drop = FALSE]
  pi_selected <- pi_i[selected]
  w_selected  <- 1 / pi_selected
  
  # Design effect
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

#' Compute Structural Hamming Distance between two DAGs
#'
#' SHD = number of edge additions + deletions + reversals needed to
#' transform est_gam into true_gam.
#'
#' @param true_gam  True adjacency matrix (q x q)
#' @param est_gam   Estimated adjacency matrix (q x q)
#' @return Integer SHD
compute_shd <- function(true_gam, est_gam) {
  true_gam <- (true_gam != 0) * 1L
  est_gam  <- (est_gam != 0) * 1L
  q <- nrow(true_gam)
  
  shd <- 0
  
  for (i in seq_len(q)) {
    for (j in seq_len(q)) {
      if (i == j) next
      
      true_ij <- true_gam[i, j]
      true_ji <- true_gam[j, i]
      est_ij  <- est_gam[i, j]
      est_ji  <- est_gam[j, i]
      
      # Only count each pair once (i < j)
      if (i > j) next
      
      # True edge status: j->i if true_gam[i,j]=1, i->j if true_gam[j,i]=1
      true_edge  <- true_ij + true_ji  # 0 = no edge, 1 = one direction
      est_edge   <- est_ij + est_ji
      
      if (true_edge == 0 && est_edge == 0) {
        # Both absent: no error
      } else if (true_edge == 0 && est_edge > 0) {
        # Spurious edge: +1
        shd <- shd + 1
      } else if (true_edge > 0 && est_edge == 0) {
        # Missing edge: +1
        shd <- shd + 1
      } else {
        # Both have edge: check direction
        if (true_ij != est_ij || true_ji != est_ji) {
          # Wrong direction: +1 (reversal)
          shd <- shd + 1
        }
      }
    }
  }
  
  return(shd)
}


#' Compute detailed edge metrics
#'
#' @return List with SHD, TPR, FPR, precision, orientation_accuracy
compute_edge_metrics <- function(true_gam, est_gam) {
  true_gam <- (true_gam != 0) * 1L
  est_gam  <- (est_gam != 0) * 1L
  q <- nrow(true_gam)
  
  # --- Skeleton metrics (ignore direction) ---
  true_skel <- pmax(true_gam, t(true_gam))
  est_skel  <- pmax(est_gam, t(est_gam))
  
  # Use upper triangle only
  ut <- upper.tri(true_gam)
  
  true_edges  <- true_skel[ut]
  est_edges   <- est_skel[ut]
  
  tp_skel <- sum(true_edges == 1 & est_edges == 1)
  fp_skel <- sum(true_edges == 0 & est_edges == 1)
  fn_skel <- sum(true_edges == 1 & est_edges == 0)
  tn_skel <- sum(true_edges == 0 & est_edges == 0)
  
  tpr_skel <- ifelse(tp_skel + fn_skel > 0, tp_skel / (tp_skel + fn_skel), NA)
  fpr_skel <- ifelse(fp_skel + tn_skel > 0, fp_skel / (fp_skel + tn_skel), NA)
  precision_skel <- ifelse(tp_skel + fp_skel > 0, tp_skel / (tp_skel + fp_skel), NA)
  
  # --- Orientation accuracy (among correctly detected edges) ---
  n_correct_skel <- 0
  n_correct_orient <- 0
  
  for (i in seq_len(q)) {
    for (j in (i + 1):min(q, q)) {
      if (j > q) next
      # Both have this edge in skeleton?
      if (true_skel[i, j] == 1 && est_skel[i, j] == 1) {
        n_correct_skel <- n_correct_skel + 1
        # Same orientation?
        if (true_gam[i, j] == est_gam[i, j] && true_gam[j, i] == est_gam[j, i]) {
          n_correct_orient <- n_correct_orient + 1
        }
      }
    }
  }
  
  orient_acc <- ifelse(n_correct_skel > 0, n_correct_orient / n_correct_skel, NA)
  
  # --- Correct orientation rate (denominated by TOTAL true edges) ---
  # This penalizes both missing edges and wrong orientations,
  # avoiding the inflation artifact where biased samples make
  # detected edges easier to orient.
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
# SECTION 4: SINGLE REPLICATION RUNNER
# ==============================================================================

#' Run a single replication of the simulation
#'
#' Generates population, samples, fits naive and swa-oBN, returns metrics.
#'
#' @param gam        True DAG adjacency matrix
#' @param N          Population size
#' @param sigma      Signal strength
#' @param L          Ordinal levels
#' @param design     Survey design: "none", "mild", "moderate", "extreme"
#' @param target_n   Approximate sample size
#' @param nstart     Number of restarts for hill-climbing
#' @param maxit      Max iterations per restart
#' @param Z_pop      Population covariate matrix (NULL = no covariates)
#' @param delta_strength  Covariate effect strength in DGP
#' @param oversample_vars Which variables determine oversampling
#' @param swa_obn_fn The swa_obn function (pass it in to avoid sourcing issues)
#'
#' @return Data frame row with metrics for both methods
run_single_replication <- function(gam, N = 20000, sigma = 1.0, L = 3,
                                   design = "moderate", target_n = 1000,
                                   nstart = 10, maxit = 100,
                                   Z_pop = NULL, delta_strength = 0,
                                   oversample_vars = NULL,
                                   swa_obn_fn = swa_obn) {
  q <- nrow(gam)
  
  # --- Generate population ---
  pop <- generate_obn_data(
    n = N, gam = gam, sigma = sigma, L = L,
    Z = Z_pop, delta_strength = delta_strength
  )
  
  # --- Apply survey sampling ---
  if (is.null(oversample_vars)) {
    oversample_vars <- c(1, q)  # first and last node
  }
  
  svy <- apply_survey_sampling(
    pop_data = pop$data,
    design   = design,
    oversample_vars = oversample_vars,
    target_n = target_n
  )
  
  y_sample <- svy$sample_data
  w_sample <- svy$weights
  n_actual <- nrow(y_sample)
  
  # Covariates for sampled units
  Z_sample <- NULL
  if (!is.null(Z_pop)) {
    Z_sample <- as.data.frame(Z_pop[svy$pop_indices, , drop = FALSE])
  }
  
  # --- Fit naive oBN (no weights) ---
  fit_naive <- tryCatch({
    swa_obn_fn(
      y = y_sample, Z = Z_sample, weights = NULL,
      search = "greedy", ic = "bic", link = "probit",
      nstart = nstart, boot = NULL, verbose = FALSE, maxit = maxit
    )
  }, error = function(e) list(gam = matrix(0, q, q)))
  
  # --- Fit swa-oBN (with weights) ---
  fit_svy <- tryCatch({
    swa_obn_fn(
      y = y_sample, Z = Z_sample, weights = w_sample,
      search = "greedy", ic = "bic", link = "probit",
      nstart = nstart, boot = NULL, verbose = FALSE, maxit = maxit
    )
  }, error = function(e) list(gam = matrix(0, q, q)))
  
  # --- Compute metrics ---
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