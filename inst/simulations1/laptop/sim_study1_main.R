# ==============================================================================
# sim_study1_main.R
# Study 1: Survey-weighted oBN vs naive oBN
#
# Replicates the structure of Ni et al. (2025) Figure 4, adding the
# survey-weighting dimension. Three sub-studies:
#   1a. Fixed n, fixed sigma, vary design strength
#   1b. Fixed design, fixed sigma, vary sample size
#   1c. Fixed design, fixed n, vary signal strength
#
# Output: results table + plots analogous to Ni et al. Figure 4
# ==============================================================================

rm(list = ls())

# --- Source dependencies ---
# Adjust paths as needed for your project structure
source("sim_utils.R")           # Data generation, sampling, metrics
source("2026_03_11_step2.R")    # swa_obn function

library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(2026)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

N_POP       <- 20000   # Population size per replication
N_REPS      <- 50      # Replications per scenario
Q_NODES     <- 8       # Number of DAG nodes
L_LEVELS    <- 3       # Ordinal levels per node (matching CAHPS 3-level items)
EDGE_PROB   <- 0.25    # DAG density (Erdos-Renyi probability)
NSTART      <- 2      # Hill-climbing restarts
MAXIT       <- 10     # Max iterations per restart

# Generate a FIXED true DAG used across all replications
# (so SHD is comparable across scenarios, as in Ni et al.)
set.seed(42)
TRUE_DAG <- generate_random_dag(Q_NODES, prob = EDGE_PROB)
cat("True DAG edges:", sum(TRUE_DAG), "\n")
cat("True DAG adjacency matrix:\n")
print(TRUE_DAG)

# Which variables determine oversampling:
# Use first and last node (simulates oversampling on a "demographic" 
# variable and an "outcome" variable)
OVERSAMPLE_VARS <- c(1, Q_NODES)


# ==============================================================================
# STUDY 1a: VARY DESIGN STRENGTH
# Fixed n ~ 1000, sigma = 1.5
# ==============================================================================

cat("\n========================================\n")
cat("STUDY 1a: Varying survey design strength\n")
cat("========================================\n")

designs   <- c("none", "mild", "moderate", "extreme")
sigma_1a  <- 1.5
target_n_1a <- 1000

results_1a <- data.frame()

for (des in designs) {
  cat(sprintf("\n  Design: %s\n", des))
  
  for (rep in seq_len(N_REPS)) {
    if (rep %% 10 == 0) cat(sprintf("    rep %d / %d\n", rep, N_REPS))
    
    res <- tryCatch({
      run_single_replication(
        gam = TRUE_DAG, N = N_POP, sigma = sigma_1a, L = L_LEVELS,
        design = des, target_n = target_n_1a,
        nstart = NSTART, maxit = MAXIT,
        oversample_vars = OVERSAMPLE_VARS,
        swa_obn_fn = swa_obn
      )
    }, error = function(e) {
      cat(sprintf("    ERROR in rep %d: %s\n", rep, e$message))
      NULL
    })
    
    if (!is.null(res)) {
      res$design <- des
      res$sigma  <- sigma_1a
      res$rep    <- rep
      res$study  <- "1a"
      results_1a <- rbind(results_1a, res)
    }
  }
}


# ==============================================================================
# STUDY 1b: VARY SAMPLE SIZE
# Fixed design = "moderate", sigma = 1.5
# ==============================================================================

cat("\n========================================\n")
cat("STUDY 1b: Varying sample size\n")
cat("========================================\n")

target_ns  <- c(300, 500, 1000, 2000, 4000)
sigma_1b   <- 1.5
design_1b  <- "moderate"

results_1b <- data.frame()

for (tn in target_ns) {
  cat(sprintf("\n  Target n: %d\n", tn))
  
  # Increase population size for larger samples
  N_pop_scaled <- max(N_POP, tn * 10)
  
  for (rep in seq_len(N_REPS)) {
    if (rep %% 10 == 0) cat(sprintf("    rep %d / %d\n", rep, N_REPS))
    
    res <- tryCatch({
      run_single_replication(
        gam = TRUE_DAG, N = N_pop_scaled, sigma = sigma_1b, L = L_LEVELS,
        design = design_1b, target_n = tn,
        nstart = NSTART, maxit = MAXIT,
        oversample_vars = OVERSAMPLE_VARS,
        swa_obn_fn = swa_obn
      )
    }, error = function(e) {
      cat(sprintf("    ERROR in rep %d: %s\n", rep, e$message))
      NULL
    })
    
    if (!is.null(res)) {
      res$design    <- design_1b
      res$sigma     <- sigma_1b
      res$target_n  <- tn
      res$rep       <- rep
      res$study     <- "1b"
      results_1b <- rbind(results_1b, res)
    }
  }
}


# ==============================================================================
# STUDY 1c: VARY SIGNAL STRENGTH
# Fixed design = "moderate", target n ~ 1000
# ==============================================================================

cat("\n========================================\n")
cat("STUDY 1c: Varying signal strength\n")
cat("========================================\n")

sigmas     <- c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0)
design_1c  <- "moderate"
target_n_1c <- 1000

results_1c <- data.frame()

for (sig in sigmas) {
  cat(sprintf("\n  Sigma: %.2f\n", sig))
  
  for (rep in seq_len(N_REPS)) {
    if (rep %% 10 == 0) cat(sprintf("    rep %d / %d\n", rep, N_REPS))
    
    res <- tryCatch({
      run_single_replication(
        gam = TRUE_DAG, N = N_POP, sigma = sig, L = L_LEVELS,
        design = design_1c, target_n = target_n_1c,
        nstart = NSTART, maxit = MAXIT,
        oversample_vars = OVERSAMPLE_VARS,
        swa_obn_fn = swa_obn
      )
    }, error = function(e) {
      cat(sprintf("    ERROR in rep %d: %s\n", rep, e$message))
      NULL
    })
    
    if (!is.null(res)) {
      res$design <- design_1c
      res$sigma  <- sig
      res$rep    <- rep
      res$study  <- "1c"
      results_1c <- rbind(results_1c, res)
    }
  }
}


# ==============================================================================
# COMBINE AND SUMMARIZE
# ==============================================================================

results_all <- bind_rows(results_1a, results_1b, results_1c)

# Save raw results
save(results_all, TRUE_DAG, file = "sim_study1_results.RData")

# --- Summary tables ---
summarize_results <- function(df, grouping_var) {
  df %>%
    group_by(across(all_of(grouping_var))) %>%
    summarize(
      n_reps         = n(),
      n_actual_mean  = mean(n_actual),
      n_eff_mean     = mean(n_eff, na.rm = TRUE),
      deff_mean      = mean(deff, na.rm = TRUE),
      
      shd_naive_mean = mean(shd_naive),
      shd_naive_se   = sd(shd_naive) / sqrt(n()),
      shd_svy_mean   = mean(shd_svy),
      shd_svy_se     = sd(shd_svy) / sqrt(n()),
      
      tpr_naive_mean = mean(tpr_naive, na.rm = TRUE),
      tpr_svy_mean   = mean(tpr_svy, na.rm = TRUE),
      
      fpr_naive_mean = mean(fpr_naive, na.rm = TRUE),
      fpr_svy_mean   = mean(fpr_svy, na.rm = TRUE),
      
      prec_naive_mean = mean(prec_naive, na.rm = TRUE),
      prec_svy_mean   = mean(prec_svy, na.rm = TRUE),
      
      orient_naive_mean = mean(orient_naive, na.rm = TRUE),
      orient_svy_mean   = mean(orient_svy, na.rm = TRUE),
      
      cor_naive_mean = mean(cor_naive, na.rm = TRUE),
      cor_svy_mean   = mean(cor_svy, na.rm = TRUE),
      
      .groups = "drop"
    )
}

cat("\n\n========================================\n")
cat("STUDY 1a SUMMARY: Varying design strength\n")
cat("========================================\n")
summary_1a <- summarize_results(results_1a, "design")
print(as.data.frame(summary_1a), digits = 3)

cat("\n\n========================================\n")
cat("STUDY 1b SUMMARY: Varying sample size\n")
cat("========================================\n")
summary_1b <- summarize_results(results_1b, "target_n")
print(as.data.frame(summary_1b), digits = 3)

cat("\n\n========================================\n")
cat("STUDY 1c SUMMARY: Varying signal strength\n")
cat("========================================\n")
summary_1c <- summarize_results(results_1c, "sigma")
print(as.data.frame(summary_1c), digits = 3)


# ==============================================================================
# PLOTS
# ==============================================================================

# --- Helper: reshape to long format for ggplot ---
reshape_for_plot <- function(summary_df, x_var, x_label) {
  summary_df %>%
    select(all_of(c(x_var, "shd_naive_mean", "shd_naive_se",
                    "shd_svy_mean", "shd_svy_se"))) %>%
    pivot_longer(
      cols = starts_with("shd_"),
      names_to = c("metric", "method", ".value"),
      names_pattern = "shd_(naive|svy)_(mean|se)"
    ) %>%
    mutate(
      method = ifelse(method == "naive", "Naive oBN", "swa-oBN"),
      x = .data[[x_var]]
    )
}

# --- Plot 1a: SHD vs design strength ---
plot_1a <- summary_1a %>%
  mutate(design = factor(design, levels = c("none", "mild", "moderate", "extreme"))) %>%
  pivot_longer(
    cols = c(shd_naive_mean, shd_svy_mean),
    names_to = "method", values_to = "shd"
  ) %>%
  mutate(
    se = ifelse(method == "shd_naive_mean",
                summary_1a$shd_naive_se[match(design, summary_1a$design)],
                summary_1a$shd_svy_se[match(design, summary_1a$design)]),
    method = ifelse(method == "shd_naive_mean", "Naive oBN", "swa-oBN")
  ) %>%
  ggplot(aes(x = design, y = shd, color = method, group = method)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  geom_line(position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = shd - 1.96 * se, ymax = shd + 1.96 * se),
                width = 0.15, position = position_dodge(0.2)) +
  labs(x = "Survey Design Strength",
       y = "Structural Hamming Distance",
       color = "Method",
       title = "Study 1a: Effect of Survey Design on DAG Recovery",
       subtitle = sprintf("q = %d nodes, n ≈ %d, σ = %.1f, %d replications",
                          Q_NODES, target_n_1a, sigma_1a, N_REPS)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

ggsave("sim_study1a_design.pdf", plot_1a, width = 8, height = 6)

# --- Plot 1b: SHD vs sample size ---
plot_data_1b <- summary_1b %>%
  pivot_longer(
    cols = c(shd_naive_mean, shd_svy_mean),
    names_to = "method", values_to = "shd"
  ) %>%
  mutate(
    se = ifelse(method == "shd_naive_mean",
                summary_1b$shd_naive_se[match(target_n, summary_1b$target_n)],
                summary_1b$shd_svy_se[match(target_n, summary_1b$target_n)]),
    method = ifelse(method == "shd_naive_mean", "Naive oBN", "swa-oBN")
  )

plot_1b <- plot_data_1b %>%
  ggplot(aes(x = target_n, y = shd, color = method, group = method)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = shd - 1.96 * se, ymax = shd + 1.96 * se),
                width = 0.05) +
  scale_x_log10(breaks = target_ns) +
  labs(x = "Sample Size (log scale)",
       y = "Structural Hamming Distance",
       color = "Method",
       title = "Study 1b: Effect of Sample Size on DAG Recovery",
       subtitle = sprintf("q = %d, σ = %.1f, design = %s, %d replications",
                          Q_NODES, sigma_1b, design_1b, N_REPS)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

ggsave("sim_study1b_samplesize.pdf", plot_1b, width = 8, height = 6)

# --- Plot 1c: SHD vs signal strength ---
plot_data_1c <- summary_1c %>%
  pivot_longer(
    cols = c(shd_naive_mean, shd_svy_mean),
    names_to = "method", values_to = "shd"
  ) %>%
  mutate(
    se = ifelse(method == "shd_naive_mean",
                summary_1c$shd_naive_se[match(sigma, summary_1c$sigma)],
                summary_1c$shd_svy_se[match(sigma, summary_1c$sigma)]),
    method = ifelse(method == "shd_naive_mean", "Naive oBN", "swa-oBN")
  )

plot_1c <- plot_data_1c %>%
  ggplot(aes(x = sigma, y = shd, color = method, group = method)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = shd - 1.96 * se, ymax = shd + 1.96 * se),
                width = 0.05) +
  labs(x = expression("Signal Strength " * sigma),
       y = "Structural Hamming Distance",
       color = "Method",
       title = "Study 1c: Effect of Signal Strength on DAG Recovery",
       subtitle = sprintf("q = %d, n ≈ %d, design = %s, %d replications",
                          Q_NODES, target_n_1c, design_1c, N_REPS)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

ggsave("sim_study1c_signal.pdf", plot_1c, width = 8, height = 6)


# --- Precision plot (primary: shows swa-oBN avoids spurious edges) ---
prec_data <- bind_rows(
  summary_1a %>%
    filter(design != "extreme") %>%
    mutate(x = design, study = "Design Strength"),
  summary_1c %>% mutate(x = as.character(sigma), study = "Signal Strength")
) %>%
  select(x, study, prec_naive_mean, prec_svy_mean) %>%
  pivot_longer(cols = starts_with("prec_"),
               names_to = "method", values_to = "precision") %>%
  mutate(method = ifelse(method == "prec_naive_mean", "Naive oBN", "swa-oBN"))

plot_precision <- prec_data %>%
  ggplot(aes(x = x, y = precision, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ study, scales = "free_x") +
  labs(x = "", y = "Precision (TP / (TP + FP))",
       fill = "Method",
       title = "Edge Precision: swa-oBN vs Naive oBN",
       subtitle = "Fraction of detected edges that are true edges") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 1))

ggsave("sim_study1_precision.pdf", plot_precision, width = 12, height = 6)

# --- Correct orientation rate (found AND correctly oriented / total true edges) ---
cor_data <- bind_rows(
  summary_1a %>%
    filter(design != "extreme") %>%
    mutate(x = design, study = "Design Strength"),
  summary_1c %>% mutate(x = as.character(sigma), study = "Signal Strength")
) %>%
  select(x, study, cor_naive_mean, cor_svy_mean) %>%
  pivot_longer(cols = starts_with("cor_"),
               names_to = "method", values_to = "rate") %>%
  mutate(method = ifelse(method == "cor_naive_mean", "Naive oBN", "swa-oBN"))

plot_cor <- cor_data %>%
  ggplot(aes(x = x, y = rate, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ study, scales = "free_x") +
  labs(x = "", y = "Correct Orientation Rate",
       fill = "Method",
       title = "Correct Orientation Rate: swa-oBN vs Naive oBN",
       subtitle = "Fraction of true edges found AND correctly oriented") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 1))

ggsave("sim_study1_orientation.pdf", plot_cor, width = 12, height = 6)


cat("\n\nStudy 1 complete. Results saved to sim_study1_results.RData\n")
cat("Plots saved to sim_study1a_design.pdf, sim_study1b_samplesize.pdf,\n")
cat("  sim_study1c_signal.pdf, sim_study1_precision.pdf, sim_study1_orientation.pdf\n")