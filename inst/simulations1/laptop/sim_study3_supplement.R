# ==============================================================================
# sim_study3_supplement.R
# Study 3 (Supplementary): Survey weighting for other algorithm families
#
# Demonstrates that the survey-weighting problem is general, not oBN-specific.
# Three sub-studies:
#   3a. swa-PC vs naive PC (constraint-based, continuous)
#   3b. swa-FCI vs naive FCI (constraint-based with latent variables)
#   3c. swa-LiNGAM vs naive LiNGAM (asymmetry-based, continuous non-Gaussian)
#
# Each uses the same true chain X -> Y -> Z with informative survey sampling
# that creates a spurious X -- Z association. Replicated N_REPS times.
# Reports: skeleton SHD, spurious edge rate, and correct skeleton rate.
#
# NOTE: This is for the supplement. The main paper focuses on oBN (Studies 1-2).
# ==============================================================================

rm(list = ls())

library(pcalg)
library(survey)
library(MASS)

set.seed(2026)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

N_POP    <- 20000   # Population size
N_REPS   <- 100     # Replications (more reps since each is fast)
DESIGNS  <- c("none", "mild", "moderate", "extreme")

# Design strength parameters (probability of selection)
DESIGN_PARAMS <- list(
  none     = c(target = 0.10, other = 0.10),
  mild     = c(target = 0.30, other = 0.08),
  moderate = c(target = 0.60, other = 0.05),
  extreme  = c(target = 1.00, other = 0.03)
)


# ==============================================================================
# HELPER: Survey sampling for 3-variable continuous/ordinal data
# ==============================================================================

sample_with_design <- function(pop_data, design_name) {
  N <- nrow(pop_data)
  
  # Oversampling rule: top quartile of X + Z combined
  target_score <- scale(as.numeric(pop_data[, 1])) + scale(as.numeric(pop_data[, 3]))
  is_target <- target_score > quantile(target_score, 0.75)
  
  params <- DESIGN_PARAMS[[design_name]]
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


# ==============================================================================
# SURVEY-WEIGHTED CI TEST (for PC and FCI)
# ==============================================================================

svy_ci_test <- function(x, y, S, suffStat) {
  dat <- suffStat$data
  wts <- suffStat$weights
  
  dat$.wt <- wts
  des <- svydesign(ids = ~1, weights = ~.wt, data = dat)
  
  var_x <- colnames(dat)[x]
  var_y <- colnames(dat)[y]
  
  if (length(S) == 0) {
    form <- as.formula(paste(var_x, "~", var_y))
  } else {
    var_S <- colnames(dat)[S]
    form <- as.formula(paste(var_x, "~", var_y, "+",
                              paste(var_S, collapse = " + ")))
  }
  
  fit <- try(svyglm(form, design = des), silent = TRUE)
  if (inherits(fit, "try-error")) return(1)
  
  coefs <- summary(fit)$coefficients
  if (var_y %in% rownames(coefs)) {
    return(coefs[var_y, "Pr(>|t|)"])
  } else {
    return(1)
  }
}


# ==============================================================================
# svy-LiNGAM via pseudo-population bootstrap
# ==============================================================================

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


# ==============================================================================
# METRIC: Skeleton match for 3-node chain
# ==============================================================================

# True skeleton: X-Y, Y-Z, no X-Z
# Adjacency for true chain X->Y->Z:
#   true_adj[1,2]=1, true_adj[2,3]=1 (or symmetric for skeleton)

skeleton_metrics_3node <- function(adj_mat) {
  # Convert to skeleton (symmetric)
  skel <- pmax(adj_mat, t(adj_mat))
  skel <- (skel != 0) * 1L
  
  # True skeleton: 1-2 present, 2-3 present, 1-3 absent
  has_12 <- skel[1, 2] == 1
  has_23 <- skel[2, 3] == 1
  has_13 <- skel[1, 3] == 1  # this is the spurious edge
  
  # Skeleton SHD
  shd <- 0
  if (!has_12) shd <- shd + 1  # missing true edge
  if (!has_23) shd <- shd + 1  # missing true edge
  if (has_13)  shd <- shd + 1  # spurious edge
  
  return(list(
    shd         = shd,
    has_12      = has_12,
    has_23      = has_23,
    spurious_13 = has_13,
    correct     = (has_12 & has_23 & !has_13)
  ))
}


# ==============================================================================
# STUDY 3a: swa-PC vs naive PC (Gaussian data)
# ==============================================================================

cat("========================================\n")
cat("STUDY 3a: swa-PC vs naive PC\n")
cat("========================================\n")

results_3a <- data.frame()

for (des in DESIGNS) {
  cat(sprintf("\n  Design: %s\n", des))
  
  for (rep in seq_len(N_REPS)) {
    if (rep %% 25 == 0) cat(sprintf("    rep %d / %d\n", rep, N_REPS))
    
    # Generate Gaussian chain: X -> Y -> Z
    X <- rnorm(N_POP)
    Y <- 1.5 * X + rnorm(N_POP)
    Z <- 1.5 * Y + rnorm(N_POP)
    pop <- data.frame(X = X, Y = Y, Z = Z)
    
    # Sample
    svy <- sample_with_design(pop, des)
    sdat <- svy$data
    
    # Naive PC (unweighted Gaussian CI test)
    suffStat_naive <- list(C = cor(sdat), n = nrow(sdat))
    pc_naive <- tryCatch(
      pc(suffStat = suffStat_naive, indepTest = gaussCItest,
         alpha = 0.05, labels = colnames(sdat), verbose = FALSE),
      error = function(e) NULL
    )
    
    # swa-PC (survey-weighted CI test)
    suffStat_svy <- list(data = sdat, weights = svy$weights)
    pc_svy <- tryCatch(
      pc(suffStat = suffStat_svy, indepTest = svy_ci_test,
         alpha = 0.05, labels = colnames(sdat), verbose = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(pc_naive) && !is.null(pc_svy)) {
      adj_naive <- as(pc_naive@graph, "matrix")
      adj_svy   <- as(pc_svy@graph, "matrix")
      
      m_naive <- skeleton_metrics_3node(adj_naive)
      m_svy   <- skeleton_metrics_3node(adj_svy)
      
      results_3a <- rbind(results_3a, data.frame(
        design = des, rep = rep,
        n = svy$n, n_eff = svy$n_eff, deff = svy$deff,
        shd_naive = m_naive$shd, spurious_naive = m_naive$spurious_13,
        correct_naive = m_naive$correct,
        shd_svy = m_svy$shd, spurious_svy = m_svy$spurious_13,
        correct_svy = m_svy$correct
      ))
    }
  }
}


# ==============================================================================
# STUDY 3b: swa-FCI vs naive FCI (Gaussian data)
# ==============================================================================

cat("\n========================================\n")
cat("STUDY 3b: swa-FCI vs naive FCI\n")
cat("========================================\n")

results_3b <- data.frame()

for (des in DESIGNS) {
  cat(sprintf("\n  Design: %s\n", des))
  
  for (rep in seq_len(N_REPS)) {
    if (rep %% 25 == 0) cat(sprintf("    rep %d / %d\n", rep, N_REPS))
    
    X <- rnorm(N_POP)
    Y <- 0.8 * X + rnorm(N_POP)
    Z <- 0.8 * Y + rnorm(N_POP)
    pop <- data.frame(X = X, Y = Y, Z = Z)
    
    svy <- sample_with_design(pop, des)
    sdat <- svy$data
    
    # Naive FCI
    suffStat_naive <- list(C = cor(sdat), n = nrow(sdat))
    fci_naive <- tryCatch(
      fci(suffStat = suffStat_naive, indepTest = gaussCItest,
          alpha = 0.05, labels = colnames(sdat), p = 3, verbose = FALSE),
      error = function(e) NULL
    )
    
    # swa-FCI
    suffStat_svy <- list(data = sdat, weights = svy$weights)
    fci_svy <- tryCatch(
      fci(suffStat = suffStat_svy, indepTest = svy_ci_test,
          alpha = 0.05, labels = colnames(sdat), p = 3, verbose = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(fci_naive) && !is.null(fci_svy)) {
      adj_naive <- fci_naive@amat
      adj_svy   <- fci_svy@amat
      
      # For PAGs, any non-zero entry means edge present
      adj_naive_bin <- (adj_naive != 0) * 1L
      adj_svy_bin   <- (adj_svy != 0) * 1L
      
      m_naive <- skeleton_metrics_3node(adj_naive_bin)
      m_svy   <- skeleton_metrics_3node(adj_svy_bin)
      
      results_3b <- rbind(results_3b, data.frame(
        design = des, rep = rep,
        n = svy$n, n_eff = svy$n_eff, deff = svy$deff,
        shd_naive = m_naive$shd, spurious_naive = m_naive$spurious_13,
        correct_naive = m_naive$correct,
        shd_svy = m_svy$shd, spurious_svy = m_svy$spurious_13,
        correct_svy = m_svy$correct
      ))
    }
  }
}


# ==============================================================================
# STUDY 3c: swa-LiNGAM vs naive LiNGAM (non-Gaussian data)
# ==============================================================================

cat("\n========================================\n")
cat("STUDY 3c: swa-LiNGAM vs naive LiNGAM\n")
cat("========================================\n")

results_3c <- data.frame()

for (des in DESIGNS) {
  cat(sprintf("\n  Design: %s\n", des))
  
  for (rep in seq_len(N_REPS)) {
    if (rep %% 25 == 0) cat(sprintf("    rep %d / %d\n", rep, N_REPS))
    
    # Non-Gaussian noise (uniform) -- required for LiNGAM identifiability
    X <- runif(N_POP, -2, 2)
    Y <- 1.5 * X + runif(N_POP, -2, 2)
    Z <- 1.5 * Y + runif(N_POP, -2, 2)
    pop <- data.frame(X = X, Y = Y, Z = Z)
    
    svy <- sample_with_design(pop, des)
    sdat <- svy$data
    
    # Naive LiNGAM
    lingam_naive <- tryCatch({
      res <- lingam(as.matrix(sdat), verbose = FALSE)
      adj <- (t(res$Bpruned) != 0) * 1L
      rownames(adj) <- colnames(adj) <- colnames(sdat)
      adj
    }, error = function(e) NULL)
    
    # swa-LiNGAM
    lingam_svy <- tryCatch(
      svy_lingam(as.matrix(sdat), svy$weights, B = 50),
      error = function(e) NULL
    )
    
    if (!is.null(lingam_naive) && !is.null(lingam_svy)) {
      m_naive <- skeleton_metrics_3node(lingam_naive)
      m_svy   <- skeleton_metrics_3node(lingam_svy)
      
      results_3c <- rbind(results_3c, data.frame(
        design = des, rep = rep,
        n = svy$n, n_eff = svy$n_eff, deff = svy$deff,
        shd_naive = m_naive$shd, spurious_naive = m_naive$spurious_13,
        correct_naive = m_naive$correct,
        shd_svy = m_svy$shd, spurious_svy = m_svy$spurious_13,
        correct_svy = m_svy$correct
      ))
    }
  }
}


# ==============================================================================
# SUMMARIZE AND SAVE
# ==============================================================================

save(results_3a, results_3b, results_3c, file = "sim_study3_results.RData")

library(dplyr)
library(tidyr)
library(ggplot2)

summarize_study3 <- function(df, study_name) {
  df %>%
    group_by(design) %>%
    summarize(
      n_reps = n(),
      n_mean = mean(n),
      deff_mean = mean(deff),
      
      shd_naive_mean     = mean(shd_naive),
      shd_svy_mean       = mean(shd_svy),
      spurious_naive_rate = mean(spurious_naive),
      spurious_svy_rate   = mean(spurious_svy),
      correct_naive_rate  = mean(correct_naive),
      correct_svy_rate    = mean(correct_svy),
      .groups = "drop"
    ) %>%
    mutate(study = study_name)
}

summary_3a <- summarize_study3(results_3a, "PC")
summary_3b <- summarize_study3(results_3b, "FCI")
summary_3c <- summarize_study3(results_3c, "LiNGAM")

summary_all <- bind_rows(summary_3a, summary_3b, summary_3c) %>%
  mutate(design = factor(design, levels = c("none", "mild", "moderate", "extreme")))

cat("\n========================================\n")
cat("STUDY 3 SUMMARY: Spurious edge rate (X--Z)\n")
cat("========================================\n")
print(as.data.frame(summary_all %>%
  select(study, design, spurious_naive_rate, spurious_svy_rate,
         correct_naive_rate, correct_svy_rate)), digits = 3)


# --- Combined plot: spurious edge rate by design and algorithm family ---
plot_data <- summary_all %>%
  select(study, design, spurious_naive_rate, spurious_svy_rate) %>%
  pivot_longer(
    cols = starts_with("spurious_"),
    names_to = "method", values_to = "rate"
  ) %>%
  mutate(method = ifelse(method == "spurious_naive_rate",
                          "Naive (unweighted)", "Survey-weighted"))

plot_3 <- plot_data %>%
  ggplot(aes(x = design, y = rate, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ study) +
  labs(x = "Survey Design Strength",
       y = "Spurious Edge Rate (X -- Z)",
       fill = "Method",
       title = "Study 3: Survey Weighting Reduces Spurious Edges Across Algorithm Families",
       subtitle = sprintf("True graph: X → Y → Z. %d replications per scenario.", N_REPS)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 1))

ggsave("sim_study3_supplement.pdf", plot_3, width = 12, height = 5)


# --- Correct skeleton recovery rate ---
plot_data_correct <- summary_all %>%
  select(study, design, correct_naive_rate, correct_svy_rate) %>%
  pivot_longer(
    cols = starts_with("correct_"),
    names_to = "method", values_to = "rate"
  ) %>%
  mutate(method = ifelse(method == "correct_naive_rate",
                          "Naive (unweighted)", "Survey-weighted"))

plot_3_correct <- plot_data_correct %>%
  ggplot(aes(x = design, y = rate, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ study) +
  labs(x = "Survey Design Strength",
       y = "Correct Skeleton Recovery Rate",
       fill = "Method",
       title = "Study 3: Correct Skeleton Recovery Across Algorithm Families",
       subtitle = sprintf("True graph: X → Y → Z. %d replications per scenario.", N_REPS)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 1))

ggsave("sim_study3_supplement_correct.pdf", plot_3_correct, width = 12, height = 5)


cat("\n\nStudy 3 complete. Results saved to sim_study3_results.RData\n")
cat("Plots saved to sim_study3_supplement.pdf, sim_study3_supplement_correct.pdf\n")
