# ==============================================================================
# sim_study2_covariates.R
# Study 2: Value of covariate adjustment in swa-oBN
#
# Design: A confounder Z affects both some need nodes and outcome nodes.
# If Z is not adjusted for, the learned DAG will have spurious edges
# between the confounded nodes. We compare four estimators:
#   (1) Naive oBN (no weights, no Z)
#   (2) Z-adjusted oBN (Z but no weights)
#   (3) Weighted oBN (weights but no Z)
#   (4) swa-oBN (weights + Z)  <-- full method
#
# The true DGP has Z -> {X1, X2, X5} and Z -> {X7, X8} (outcomes),
# creating confounding between needs and outcomes that Z adjustment removes.
# ==============================================================================

rm(list = ls())

source("sim_utils.R")
source("2026_03_11_step2.R")

library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(2026)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

N_POP    <- 20000
N_REPS   <- 50
Q_NODES  <- 8       # 6 "needs" + 2 "outcomes"
L_LEVELS <- 3
NSTART   <- 10
MAXIT    <- 100
SIGMA    <- 1.5

# --- Fixed true DAG ---
# Nodes 1-6: needs, Nodes 7-8: outcomes
# Structure: sparse chain among needs, needs -> outcomes
# gam[i,j] = 1 means j -> i
set.seed(42)
TRUE_DAG <- matrix(0L, Q_NODES, Q_NODES)
# Need edges: 1->2, 2->3, 3->4, 1->5, 5->6
TRUE_DAG[2, 1] <- 1  # 1 -> 2
TRUE_DAG[3, 2] <- 1  # 2 -> 3
TRUE_DAG[4, 3] <- 1  # 3 -> 4
TRUE_DAG[5, 1] <- 1  # 1 -> 5
TRUE_DAG[6, 5] <- 1  # 5 -> 6
# Need -> outcome edges: 3->7, 6->8, 4->8
TRUE_DAG[7, 3] <- 1  # 3 -> 7
TRUE_DAG[8, 6] <- 1  # 6 -> 8
TRUE_DAG[8, 4] <- 1  # 4 -> 8

cat("True DAG edges:", sum(TRUE_DAG), "\n")
cat("True DAG:\n")
print(TRUE_DAG)

# Confounder Z affects nodes 1, 2, 5 (upstream needs) and 7, 8 (outcomes)
CONFOUNDED_NODES <- c(1, 2, 5, 7, 8)
DELTA_STRENGTH   <- 1.2  # Strong confounding


# ==============================================================================
# DATA GENERATION WITH CONFOUNDING
# ==============================================================================

#' Generate population with explicit confounder Z
#'
#' Z is a continuous variable that affects CONFOUNDED_NODES.
#' If Z is not adjusted for, the associations between confounded need nodes
#' and confounded outcome nodes will be inflated, creating spurious edges.
generate_confounded_population <- function(N, gam, sigma, L, 
                                           confounded_nodes, delta_strength) {
  q <- nrow(gam)
  
  # Generate confounder
  Z_pop <- matrix(rnorm(N, 0, 1), ncol = 1)
  colnames(Z_pop) <- "Z"
  
  # Generate oBN data with Z affecting only confounded nodes
  # We do this by generating with delta_strength for all nodes,
  # then zeroing out the delta for non-confounded nodes
  pop <- generate_obn_data(
    n = N, gam = gam, sigma = sigma, L = L,
    Z = Z_pop, delta_strength = delta_strength
  )
  
  # Zero out covariate effects for non-confounded nodes
  for (j in seq_len(q)) {
    if (!(j %in% confounded_nodes)) {
      pop$params[[j]]$delta <- 0
    }
  }
  
  # Regenerate data with corrected parameters
  # (Need to regenerate because the initial generation used deltas for all nodes)
  pop2 <- generate_obn_data(
    n = N, gam = gam, sigma = sigma, L = L,
    Z = NULL, delta_strength = 0, seed = NULL
  )
  
  # Manually regenerate confounded nodes with Z effect
  # This is cleaner: generate without Z, then regenerate confounded nodes
  # Actually, let's use a different approach: generate all data from scratch
  # with selective Z effects
  
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
    
    # Add Z effect only for confounded nodes
    if (j %in% confounded_nodes && !is.null(par$delta)) {
      eta <- eta + as.numeric(Z_pop %*% par$delta)
    }
    
    # Generate from cumulative probit
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


# ==============================================================================
# SINGLE REPLICATION FOR STUDY 2
# ==============================================================================

run_study2_replication <- function(gam, N, sigma, L, design, target_n,
                                   confounded_nodes, delta_strength,
                                   nstart, maxit) {
  q <- nrow(gam)
  
  # --- Generate confounded population ---
  pop <- generate_confounded_population(
    N = N, gam = gam, sigma = sigma, L = L,
    confounded_nodes = confounded_nodes,
    delta_strength = delta_strength
  )
  
  # --- Survey sampling (oversample high X1 + high X8) ---
  svy <- apply_survey_sampling(
    pop_data = pop$data,
    design   = design,
    oversample_vars = c(1, q),
    target_n = target_n
  )
  
  y_sample <- svy$sample_data
  w_sample <- svy$weights
  Z_sample <- as.data.frame(pop$Z[svy$pop_indices, , drop = FALSE])
  colnames(Z_sample) <- "Z"
  
  # --- Fit four models ---
  
  # (1) Naive: no weights, no Z
  fit_naive <- tryCatch(
    swa_obn(y = y_sample, Z = NULL, weights = NULL,
            search = "greedy", ic = "bic", link = "probit",
            nstart = nstart, verbose = FALSE, maxit = maxit),
    error = function(e) list(gam = matrix(0, q, q))
  )
  
  # (2) Z-only: Z adjustment, no weights
  fit_z <- tryCatch(
    swa_obn(y = y_sample, Z = Z_sample, weights = NULL,
            search = "greedy", ic = "bic", link = "probit",
            nstart = nstart, verbose = FALSE, maxit = maxit),
    error = function(e) list(gam = matrix(0, q, q))
  )
  
  # (3) Weights-only: weights, no Z
  fit_w <- tryCatch(
    swa_obn(y = y_sample, Z = NULL, weights = w_sample,
            search = "greedy", ic = "bic", link = "probit",
            nstart = nstart, verbose = FALSE, maxit = maxit),
    error = function(e) list(gam = matrix(0, q, q))
  )
  
  # (4) Full swa-oBN: weights + Z
  fit_full <- tryCatch(
    swa_obn(y = y_sample, Z = Z_sample, weights = w_sample,
            search = "greedy", ic = "bic", link = "probit",
            nstart = nstart, verbose = FALSE, maxit = maxit),
    error = function(e) list(gam = matrix(0, q, q))
  )
  
  # --- Metrics ---
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
# RUN STUDY 2
# ==============================================================================

cat("\n========================================\n")
cat("STUDY 2: Covariate adjustment\n")
cat("========================================\n")

# 2a: Vary confounding strength
delta_strengths <- c(0, 0.5, 1.0, 1.5, 2.0)
results_2a <- data.frame()

for (ds in delta_strengths) {
  cat(sprintf("\n  Delta strength: %.1f\n", ds))
  
  for (rep in seq_len(N_REPS)) {
    if (rep %% 10 == 0) cat(sprintf("    rep %d / %d\n", rep, N_REPS))
    
    res <- tryCatch({
      run_study2_replication(
        gam = TRUE_DAG, N = N_POP, sigma = SIGMA, L = L_LEVELS,
        design = "moderate", target_n = 1000,
        confounded_nodes = CONFOUNDED_NODES,
        delta_strength = ds,
        nstart = NSTART, maxit = MAXIT
      )
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", e$message))
      NULL
    })
    
    if (!is.null(res)) {
      res$delta_strength <- ds
      res$rep <- rep
      results_2a <- rbind(results_2a, res)
    }
  }
}


# 2b: Vary design strength with fixed confounding
designs_2b <- c("none", "mild", "moderate", "extreme")
results_2b <- data.frame()

for (des in designs_2b) {
  cat(sprintf("\n  Design: %s (delta = %.1f)\n", des, DELTA_STRENGTH))
  
  for (rep in seq_len(N_REPS)) {
    if (rep %% 10 == 0) cat(sprintf("    rep %d / %d\n", rep, N_REPS))
    
    res <- tryCatch({
      run_study2_replication(
        gam = TRUE_DAG, N = N_POP, sigma = SIGMA, L = L_LEVELS,
        design = des, target_n = 1000,
        confounded_nodes = CONFOUNDED_NODES,
        delta_strength = DELTA_STRENGTH,
        nstart = NSTART, maxit = MAXIT
      )
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", e$message))
      NULL
    })
    
    if (!is.null(res)) {
      res$design <- des
      res$rep <- rep
      results_2b <- rbind(results_2b, res)
    }
  }
}


# ==============================================================================
# SUMMARIZE AND PLOT
# ==============================================================================

save(results_2a, results_2b, TRUE_DAG, file = "sim_study2_results.RData")

# --- Summary 2a: vary confounding ---
summary_2a <- results_2a %>%
  group_by(delta_strength) %>%
  summarize(
    n_reps = n(),
    across(starts_with("shd_"), list(mean = mean, se = ~ sd(.x)/sqrt(n())),
           .names = "{.col}_{.fn}"),
    across(starts_with("fpr_"), mean, .names = "{.col}_mean"),
    across(starts_with("nedge_"), mean, .names = "{.col}_mean"),
    .groups = "drop"
  )

cat("\n\n========================================\n")
cat("STUDY 2a: SHD by confounding strength\n")
cat("========================================\n")
print(as.data.frame(summary_2a %>% 
  select(delta_strength, ends_with("_mean")) %>%
  select(delta_strength, starts_with("shd_"))), digits = 3)

# --- Plot 2a: SHD by confounding strength, four methods ---
plot_data_2a <- summary_2a %>%
  select(delta_strength, 
         shd_naive_mean, shd_z_mean, shd_w_mean, shd_full_mean,
         shd_naive_se, shd_z_se, shd_w_se, shd_full_se) %>%
  pivot_longer(
    cols = -delta_strength,
    names_to = c("metric", "method", "stat"),
    names_pattern = "shd_(naive|z|w|full)_(mean|se)"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    method = factor(method,
      levels = c("naive", "z", "w", "full"),
      labels = c("Naive", "Z-adjusted", "Weighted", "swa-oBN (full)")
    )
  )

plot_2a <- plot_data_2a %>%
  ggplot(aes(x = delta_strength, y = mean, color = method, shape = method)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                width = 0.05) +
  labs(x = expression("Confounding Strength " * delta),
       y = "Structural Hamming Distance",
       color = "Method", shape = "Method",
       title = "Study 2a: Effect of Confounding Strength",
       subtitle = sprintf("q = %d, n ≈ 1000, σ = %.1f, moderate design, %d reps",
                          Q_NODES, SIGMA, N_REPS)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

ggsave("sim_study2a_confounding.pdf", plot_2a, width = 8, height = 6)


# --- Plot 2b: SHD by design, four methods, with strong confounding ---
summary_2b <- results_2b %>%
  group_by(design) %>%
  summarize(
    n_reps = n(),
    across(starts_with("shd_"), list(mean = mean, se = ~ sd(.x)/sqrt(n())),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  mutate(design = factor(design, levels = c("none", "mild", "moderate", "extreme")))

plot_data_2b <- summary_2b %>%
  select(design, 
         shd_naive_mean, shd_z_mean, shd_w_mean, shd_full_mean,
         shd_naive_se, shd_z_se, shd_w_se, shd_full_se) %>%
  pivot_longer(
    cols = -design,
    names_to = c("metric", "method", "stat"),
    names_pattern = "shd_(naive|z|w|full)_(mean|se)"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    method = factor(method,
      levels = c("naive", "z", "w", "full"),
      labels = c("Naive", "Z-adjusted", "Weighted", "swa-oBN (full)")
    )
  )

plot_2b <- plot_data_2b %>%
  ggplot(aes(x = design, y = mean, color = method, group = method, shape = method)) +
  geom_point(size = 3, position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                width = 0.15, position = position_dodge(0.3)) +
  labs(x = "Survey Design Strength",
       y = "Structural Hamming Distance",
       color = "Method", shape = "Method",
       title = "Study 2b: Joint Effect of Survey Design and Confounding",
       subtitle = sprintf("q = %d, n ≈ 1000, σ = %.1f, δ = %.1f, %d reps",
                          Q_NODES, SIGMA, DELTA_STRENGTH, N_REPS)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

ggsave("sim_study2b_design_confounding.pdf", plot_2b, width = 8, height = 6)


cat("\n\nStudy 2 complete. Results saved to sim_study2_results.RData\n")
