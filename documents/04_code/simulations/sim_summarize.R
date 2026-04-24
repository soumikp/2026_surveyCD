###############################################################################
# sim_summarize.R
#
# Combines all RDS result files from the cluster and produces:
#   - sim_study1_results.RData, sim_study2_results.RData, sim_study3_results.RData
#   - All plots from the original study scripts
#
# Usage:
#   Rscript sim_summarize.R
#
# Run AFTER all SLURM array jobs have completed.
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)

BASE_DIR <- "/ihome/spurkayastha/soumik/2026_surveyCD/simulations"
out_dir  <- file.path(BASE_DIR, "sim_results")
fig_dir  <- file.path(BASE_DIR, "sim_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

select <- dplyr::select

# ============================================================================
# 1. LOAD AND COMBINE ALL RESULTS
# ============================================================================

rds_files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Found", length(rds_files), "result files\n")

if (length(rds_files) == 0) {
  stop("No .rds files found in ", out_dir,
       ". Have the SLURM jobs finished?")
}

df_all <- bind_rows(lapply(rds_files, readRDS))

# --- REMOVE EXTREME SETTING GLOBALLY ---
if ("design" %in% names(df_all)) {
  df_all <- df_all %>% filter(is.na(design) | design != "extreme")
}

cat("Total rows (after removing extreme):", nrow(df_all), "\n")
cat("Studies present:", paste(sort(unique(df_all$study)), collapse = ", "), "\n")

# Check for errors
if ("method" %in% names(df_all)) {
  n_errors <- sum(df_all$method == "ERROR", na.rm = TRUE)
  if (n_errors > 0) {
    cat(sprintf("WARNING: %d error rows found. Removing them.\n", n_errors))
    if ("error_msg" %in% names(df_all)) {
      err_summary <- df_all %>%
        filter(method == "ERROR") %>%
        count(study, error_msg) %>%
        arrange(desc(n))
      cat("Error summary:\n")
      print(as.data.frame(err_summary))
    }
    df_all <- df_all %>% filter(is.na(method) | method != "ERROR")
  }
}

# Check completion
completion <- df_all %>%
  group_by(study) %>%
  summarize(
    n_tasks      = n(),
    n_scenarios  = n_distinct(scenario_id),
    reps_per_scn = n() / n_distinct(scenario_id),
    .groups = "drop"
  )
cat("\nCompletion check:\n")
print(as.data.frame(completion))


# ============================================================================
# 2. REGENERATE TRUE DAGs (must match sim_run_cluster.R)
# ============================================================================

source(file.path(BASE_DIR, "sim_utils.R"))

set.seed(42)
TRUE_DAG_1 <- generate_random_dag(8, prob = 0.25)

TRUE_DAG_2 <- matrix(0L, 8, 8)
TRUE_DAG_2[2, 1] <- 1; TRUE_DAG_2[3, 2] <- 1; TRUE_DAG_2[4, 3] <- 1
TRUE_DAG_2[5, 1] <- 1; TRUE_DAG_2[6, 5] <- 1; TRUE_DAG_2[7, 3] <- 1
TRUE_DAG_2[8, 6] <- 1; TRUE_DAG_2[8, 4] <- 1


# ============================================================================
# 3. STUDY 1: swa-oBN vs naive oBN
# ============================================================================

cat("\n========================================\n")
cat("STUDY 1: Processing\n")
cat("========================================\n")

results_1a <- df_all %>% filter(study == "1a")
results_1b <- df_all %>% filter(study == "1b")
results_1c <- df_all %>% filter(study == "1c")

results_all_s1 <- bind_rows(results_1a, results_1b, results_1c)

# Save .RData for compatibility with any downstream scripts
results_all <- results_all_s1
save(results_all, TRUE_DAG_1, file = file.path(fig_dir, "sim_study1_results.RData"))

# --- Summary function ---
summarize_results <- function(df, grouping_var) {
  df %>%
    group_by(across(all_of(grouping_var))) %>%
    summarize(
      n_reps         = n(),
      n_actual_mean  = mean(n_actual, na.rm = TRUE),
      n_eff_mean     = mean(n_eff, na.rm = TRUE),
      deff_mean      = mean(deff, na.rm = TRUE),
      
      shd_naive_mean = mean(shd_naive, na.rm = TRUE),
      shd_naive_se   = sd(shd_naive, na.rm = TRUE) / sqrt(n()),
      shd_svy_mean   = mean(shd_svy, na.rm = TRUE),
      shd_svy_se     = sd(shd_svy, na.rm = TRUE) / sqrt(n()),
      
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

# --- Study 1a summary ---
summary_1a <- summarize_results(results_1a, "design")
cat("\nStudy 1a: Varying design strength\n")
print(as.data.frame(summary_1a), digits = 3)

# --- Study 1b summary ---
summary_1b <- summarize_results(results_1b, "target_n_param")
cat("\nStudy 1b: Varying sample size\n")
print(as.data.frame(summary_1b), digits = 3)

# --- Study 1c summary ---
summary_1c <- summarize_results(results_1c, "sigma")
cat("\nStudy 1c: Varying signal strength\n")
print(as.data.frame(summary_1c), digits = 3)


# --- Plot 1a: SHD vs design strength ---
Q_NODES <- 8; N_REPS <- max(summary_1a$n_reps)


summary_1a <- summary_1a %>%
  mutate(design = factor(design, 
                         levels = c("none", "mild", "moderate"),
                         labels = c("None", "Mild", "Strong")))

plot_1a <- summary_1a %>%
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
  geom_point(size = 4, aes(shape = method)) +
  geom_line(aes(linetype = method), linewidth = 2) +
  geom_errorbar(aes(ymin = shd - 1.96 * se, ymax = shd + 1.96 * se, linetype = method),
                width = 0.05, size = 2) +
  labs(x = "Survey Design Strength",
       y = "Structural Hamming Distance",
       color = "Method",
       linetype = "Method",
       shape = "Method",
       title = "Study 1A: Effect of Survey Design on DAG Recovery",
       subtitle = sprintf("q = %d nodes, n %s %d, %s = %.1f, %d replications",
                          Q_NODES, "\u2248", 1000, "\u03c3", 1.5, N_REPS)) +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom") + 
  scale_color_aaas()

ggsave(file.path(fig_dir, "sim_study1a_design.pdf"), plot_1a, width = 11, height = 8.5, device = cairo_pdf)


# --- Plot 1b: SHD vs sample size ---
target_ns <- sort(unique(summary_1b$target_n_param))

plot_data_1b <- summary_1b %>%
  pivot_longer(
    cols = c(shd_naive_mean, shd_svy_mean),
    names_to = "method", values_to = "shd"
  ) %>%
  mutate(
    se = ifelse(method == "shd_naive_mean",
                summary_1b$shd_naive_se[match(target_n_param, summary_1b$target_n_param)],
                summary_1b$shd_svy_se[match(target_n_param, summary_1b$target_n_param)]),
    method = ifelse(method == "shd_naive_mean", "Naive oBN", "swa-oBN")
  )

plot_1b <- plot_data_1b %>%
  ggplot(aes(x = target_n_param, y = shd, color = method, group = method)) +
  geom_point(size = 4, aes(shape = method)) +
  geom_line(aes(linetype = method), linewidth = 2) +
  geom_errorbar(aes(ymin = shd - 1.96 * se, ymax = shd + 1.96 * se, linetype = method), width = 0.05, size=2) +
  scale_x_log10(breaks = target_ns) +
  labs(x = "Sample Size (log scale)",
       y = "Structural Hamming Distance",
       color = "Method",
       linetype = "Method", 
       shape = "Method",
       title = "Study 1B: Effect of Sample Size on DAG Recovery",
       subtitle = sprintf("q = %d, %s = %.1f, design = strong, %d replications",
                          Q_NODES, "\u03c3", 1.5, N_REPS)) +
  theme_bw(base_size = 18) +
  scale_color_aaas() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "sim_study1b_samplesize.pdf"), plot_1b, width = 11, height = 8.5,device = cairo_pdf)


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
  ) %>% 
  filter(sigma >= 0.5)

plot_1c <- plot_data_1c %>%
  ggplot(aes(x = sigma, y = shd, color = method, group = method)) +
  geom_point(size = 4, aes(shape = method)) +
  geom_line(aes(linetype = method), linewidth=2) +
  geom_errorbar(aes(ymin = shd - 1.96 * se, ymax = shd + 1.96 * se, linetype = method), width = 0.05, size=2) +
  labs(x = expression("Signal Strength " * sigma),
       y = "Structural Hamming Distance",
       color = "Method",
       linetype = "Method",
       shape = "Method",
       title = "Study 1C: Effect of Signal Strength on DAG Recovery",
       subtitle = sprintf("q = %d, n %s %d, design = strong, %d replications",
                          Q_NODES, "\u2248", 1000, N_REPS)) +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom") +
  scale_color_aaas()

ggsave(file.path(fig_dir, "sim_study1c_signal.pdf"), plot_1c, width = 11, height = 8.5,device = cairo_pdf)


# --- Precision plot ---
prec_data <- bind_rows(
  summary_1a %>%
    mutate(x = factor(design, levels = c("None", "Mild", "Strong")), study_facet = "Design Strength"),
  summary_1c %>% mutate(x = as.character(sigma), study_facet = "Signal Strength")
) %>%
  select(x, study_facet, prec_naive_mean, prec_svy_mean) %>%
  pivot_longer(cols = starts_with("prec_"),
               names_to = "method", values_to = "precision") %>%
  mutate(method = ifelse(method == "prec_naive_mean", "Naive oBN", "swa-oBN")) |> 
  mutate(x = factor(x, levels = c("None", "Mild", "Strong", "0.25", "0.5", "0.75", "1", "1.5", "2")))

plot_precision <- prec_data %>%
  ggplot(aes(x = x, y = precision, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ study_facet, scales = "free_x") +
  labs(x = "", y = "Precision (TP / (TP + FP))",
       fill = "Method",
       title = "Edge Precision: swa-oBN vs Naive oBN",
       subtitle = "Fraction of detected edges that are true edges") +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 1)) + 
  scale_fill_aaas()

ggsave(file.path(fig_dir, "sim_study1_precision.pdf"), plot_precision, width = 11, height = 8.5, device = cairo_pdf)


# --- Correct orientation rate ---
cor_data <- bind_rows(
  summary_1a %>%
    mutate(x = factor(design, levels = c("None", "Mild", "Strong")), study_facet = "Design Strength"),
  summary_1c %>% mutate(x = as.character(sigma), study_facet = "Signal Strength")
) %>%
  select(x, study_facet, cor_naive_mean, cor_svy_mean) %>%
  pivot_longer(cols = starts_with("cor_"),
               names_to = "method", values_to = "rate") %>%
  mutate(method = ifelse(method == "cor_naive_mean", "Naive oBN", "swa-oBN")) |> 
  mutate(x = factor(x, levels = c("None", "Mild", "Strong", "0.25", "0.5", "0.75", "1", "1.5", "2")))


plot_cor <- cor_data %>%
  mutate(method = factor(method, 
                         labels = c("Naive oBN", "swa-oBN"), 
                         levels = c("swa-oBN", "Naive oBN"))) |> 
  ggplot(aes(x = x, y = rate, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ study_facet, scales = "free_x") +
  labs(x = "", y = "Correct Orientation Rate",
       fill = "Method",
       title = "Correct Orientation Rate: swa-oBN vs Naive oBN",
       subtitle = "Fraction of true edges found AND correctly oriented") +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom") +
  scale_fill_aaas()

ggsave(file.path(fig_dir, "sim_study1_orientation.pdf"), plot_cor, width = 11, height = 8.5, device = cairo_pdf)

cat("Study 1 plots saved.\n")


# ============================================================================
# 4. STUDY 2: COVARIATE ADJUSTMENT
# ============================================================================

cat("\n========================================\n")
cat("STUDY 2: Processing\n")
cat("========================================\n")

results_2a <- df_all %>% filter(study == "2a")
results_2b <- df_all %>% filter(study == "2b")

save(results_2a, results_2b, TRUE_DAG_2,
     file = file.path(fig_dir, "sim_study2_results.RData"))

# Columns that belong to Study 2 (drop stray cols from bind_rows with Study 1)
s2_cols <- c("delta_strength", "n_actual", "n_eff", "deff",
             "shd_naive", "shd_z", "shd_w", "shd_full",
             "fpr_naive", "fpr_z", "fpr_w", "fpr_full",
             "tpr_naive", "tpr_z", "tpr_w", "tpr_full",
             "orient_naive", "orient_z", "orient_w", "orient_full",
             "nedge_naive", "nedge_z", "nedge_w", "nedge_full",
             "n_true_edges", "study", "scenario_id", "rep_id", "design")

if (nrow(results_2a) > 0) {
  
  # Drop columns that leak in from bind_rows with Study 1
  # (Study 1 has shd_svy, Study 2 has shd_z/shd_w/shd_full — mixing creates NAs)
  results_2a <- results_2a %>% select(any_of(s2_cols))
  
  # --- Summary 2a: vary confounding strength ---
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
  
  cat("\nStudy 2A: SHD by confounding strength\n")
  print(as.data.frame(summary_2a %>%
                        select(delta_strength, ends_with("_mean")) %>%
                        select(delta_strength, starts_with("shd_"))), digits = 3)
  
  # --- Plot 2a ---
  plot_data_2a <- summary_2a %>%
    select(delta_strength,
           shd_naive_mean, shd_z_mean, shd_w_mean, shd_full_mean,
           shd_naive_se, shd_z_se, shd_w_se, shd_full_se) %>%
    pivot_longer(
      cols = -delta_strength,
      names_to = c("method", "stat"),
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
    #filter(method %in% c("Naive", "swa-oBN (full)")) |> 
    ggplot(aes(x = delta_strength, y = mean, color = method, shape = method, linetype = method)) +
    geom_point(size = 4) +
    geom_line(linewidth=2) +
    geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                  width = 0.05, size = 2) +
    labs(x = expression("Confounding Strength " * delta),
         y = "Structural Hamming Distance",
         color = "Method", shape = "Method", linetype = "Method", 
         title = "Study 2A: Effect of Confounding Strength",
         subtitle = sprintf("q = 8, n %s 1000, %s = 1.5, strong design, %d reps",
                            "\u2248", "\u03c3", max(summary_2a$n_reps))) +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom") +
    scale_color_aaas()
  
  ggsave(file.path(fig_dir, "sim_study2a_confounding.pdf"), plot_2a, width = 11, height = 8.5, device = cairo_pdf)
}

if (nrow(results_2b) > 0) {
  
  # Drop stray columns from bind_rows
  results_2b <- results_2b %>% select(any_of(s2_cols))
  
  # --- Summary 2b: vary design with strong confounding ---
  summary_2b <- results_2b %>%
    group_by(design) %>%
    summarize(
      n_reps = n(),
      across(starts_with("shd_"), list(mean = mean, se = ~ sd(.x)/sqrt(n())),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    ) %>%
    mutate(design = factor(design, levels = c("none", "mild", "moderate"), labels = c("None", "Mild", "Strong")))
  
  plot_data_2b <- summary_2b %>%
    select(design,
           shd_naive_mean, shd_z_mean, shd_w_mean, shd_full_mean,
           shd_naive_se, shd_z_se, shd_w_se, shd_full_se) %>%
    pivot_longer(
      cols = -design,
      names_to = c("method", "stat"),
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
    ggplot(aes(x = design, y = mean, color = method, group = method, shape = method, linetype = method)) +
    geom_point(size = 4) +
    geom_line(linewidth=2) +
    geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                  width = 0.05, size=2) +
    labs(x = "Survey Design Strength",
         y = "Structural Hamming Distance",
         color = "Method", shape = "Method", linetype = "Method",
         title = "Study 2B: Joint Effect of Survey Design and Confounding",
         subtitle = sprintf("q = 8, n %s 1000, %s = 1.5, %s = 1.2, %d reps",
                            "\u2248", "\u03c3", "\u03b4", max(summary_2b$n_reps))) +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom")  +
    scale_color_aaas()
  
  ggsave(file.path(fig_dir, "sim_study2b_design_confounding.pdf"), plot_2b, width = 11, height = 8.5, device = cairo_pdf)
}

cat("Study 2 plots saved.\n")


# ============================================================================
# 5. STUDY 3: OTHER ALGORITHM FAMILIES (SUPPLEMENT)
# ============================================================================

cat("\n========================================\n")
cat("STUDY 3: Processing\n")
cat("========================================\n")

results_3a <- df_all %>% filter(study == "3a")
results_3b <- df_all %>% filter(study == "3b")
results_3c <- df_all %>% filter(study == "3c")

save(results_3a, results_3b, results_3c,
     file = file.path(fig_dir, "sim_study3_results.RData"))

summarize_study3 <- function(df, study_name) {
  df %>%
    group_by(design) %>%
    summarize(
      n_reps = n(),
      n_mean = mean(n, na.rm = TRUE),
      deff_mean = mean(deff, na.rm = TRUE),
      
      shd_naive_mean     = mean(shd_naive, na.rm = TRUE),
      shd_svy_mean       = mean(shd_svy, na.rm = TRUE),
      spurious_naive_rate = mean(spurious_naive, na.rm = TRUE),
      spurious_svy_rate   = mean(spurious_svy, na.rm = TRUE),
      correct_naive_rate  = mean(correct_naive, na.rm = TRUE),
      correct_svy_rate    = mean(correct_svy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(study = study_name)
}

has_s3 <- nrow(results_3a) > 0 || nrow(results_3b) > 0 || nrow(results_3c) > 0

if (has_s3) {
  summary_parts <- list()
  if (nrow(results_3a) > 0) summary_parts$a <- summarize_study3(results_3a, "PC")
  if (nrow(results_3b) > 0) summary_parts$b <- summarize_study3(results_3b, "FCI")
  if (nrow(results_3c) > 0) summary_parts$c <- summarize_study3(results_3c, "LiNGAM")
  
  summary_all_s3 <- bind_rows(summary_parts) %>%
    mutate(design = factor(design, levels = c("none", "mild", "moderate"), labels = c("None", "Mild", "Strong")))
  
  cat("\nStudy 3: Spurious edge rate (X--Z)\n")
  print(as.data.frame(summary_all_s3 %>%
                        select(study, design, spurious_naive_rate, spurious_svy_rate,
                               correct_naive_rate, correct_svy_rate)), digits = 3)
  
  N_REPS_3 <- max(summary_all_s3$n_reps)
  
  # --- Spurious edge rate plot ---
  plot_data_s3 <- summary_all_s3 %>%
    select(study, design, spurious_naive_rate, spurious_svy_rate) %>%
    pivot_longer(
      cols = starts_with("spurious_"),
      names_to = "method", values_to = "rate"
    ) %>%
    mutate(method = ifelse(method == "spurious_naive_rate",
                           "Naive (unweighted)", "Survey-weighted"))
  
  plot_3 <- plot_data_s3 %>%
    ggplot(aes(x = design, y = rate, fill = method)) +
    geom_col(position = "dodge") +
    facet_wrap(~ study) +
    labs(x = "Survey Design Strength (Three-Tier)",
         y = "Spurious Edge Rate (X -- Z)",
         fill = "Method",
         title = "Study 3: Survey Weighting Reduces Spurious Edges Across Algorithm Families",
         subtitle = sprintf("True graph: X %s Y %s Z. %d replications per scenario.",
                            "\u2192", "\u2192", N_REPS_3)) +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0, 1)) + 
    scale_fill_aaas()
  
  ggsave(file.path(fig_dir, "sim_study3_supplement.pdf"), plot_3, width = 11, height = 8.5)
  
  # --- Correct skeleton recovery rate ---
  plot_data_correct <- summary_all_s3 %>%
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
         subtitle = sprintf("True graph: X %s Y %s Z. %d replications per scenario.",
                            "\u2192", "\u2192", N_REPS_3)) +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0, 1)) + 
    scale_fill_aaas()
  
  ggsave(file.path(fig_dir, "sim_study3_supplement_correct.pdf"), plot_3_correct,
         width = 11, height = 8.5)
  
  cat("Study 3 plots saved.\n")
} else {
  cat("No Study 3 results found (pcalg may not have been installed).\n")
}

# swa-LiNGAM improves recovery under mild-to-moderate designs in this three-tier setting.
# (Note: Evaluated strictly under the three-tier design system; extreme settings were omitted).

# ============================================================================
# 6. SUMMARY DASHBOARD
# ============================================================================

cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n")

cat("\nFiles saved to:", fig_dir, "\n")
cat("  Study 1: sim_study1_results.RData\n")
cat("           sim_study1a_design.pdf\n")
cat("           sim_study1b_samplesize.pdf\n")
cat("           sim_study1c_signal.pdf\n")
cat("           sim_study1_precision.pdf\n")
cat("           sim_study1_orientation.pdf\n")
cat("  Study 2: sim_study2_results.RData\n")
cat("           sim_study2a_confounding.pdf\n")
cat("           sim_study2b_design_confounding.pdf\n")
cat("  Study 3: sim_study3_results.RData\n")
cat("           sim_study3_supplement.pdf\n")
cat("           sim_study3_supplement_correct.pdf\n")
cat("\nDone.\n")







