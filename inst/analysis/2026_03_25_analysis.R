rm(list = ls())
source(file.path(here::here(), "inst/analysis", "2026_03_25_data.R"))
source(file.path(here::here(), "inst/analysis", "2026_03_25_functions.R"))
source(file.path(here::here(), "inst/analysis", "2026_03_25_plot.R"))
pacman::p_load(ggtext)

need_vars <- colnames(data)[str_detect(colnames(data), "n_")]
outcome_vars <- colnames(data)[str_detect(colnames(data), "o_")]
dag_vars <- c(need_vars, outcome_vars) # needs first, then outcomes
covar_vars <- c("age_v4", "race_sex_v2")

# --- Build analysis objects ---
y_data <- data[, dag_vars]

# Ensure all DAG variables are ordered factors
for (v in names(y_data)) {
  if (!is.ordered(y_data[[v]])) {
    y_data[[v]] <- as.ordered(y_data[[v]])
  }
}

# Covariates (can be any type — continuous, categorical, etc.)
Z_data <- data[, covar_vars]
Z_data$race_sex_v2 <- as.factor(Z_data$race_sex_v2) # ensure factor for dummies
Z_data$age_v4 <- as.factor(Z_data$age_v4)

# Survey weights
w_data <- data$WEIGHT

cat("=== Data Summary ===\n")
cat("Observations:", nrow(y_data), "\n")
cat("DAG nodes:", ncol(y_data), "\n")
cat("  Needs:", length(need_vars), "\n")
cat("  Outcomes:", length(outcome_vars), "\n")
cat("Covariates:", ncol(Z_data), "\n")
cat("Weight range:", round(range(w_data), 2), "\n")
cat("Effective n:", round(sum(w_data)^2 / sum(w_data^2), 1), "\n\n")

cat("Levels per DAG node:\n")
for (v in dag_vars) {
  cat(sprintf(
    "  %-25s %d levels  %s\n",
    v, nlevels(y_data[[v]]),
    ifelse(nlevels(y_data[[v]]) > 2, "[ordinal - identifiable]",
           "[binary - limited]"
    )
  ))
}


# ==============================================================================
# STEP 2: BUILD EXPERT CONSTRAINTS
# ==============================================================================

# Node indices (based on column order in y_data)
node_idx <- setNames(seq_along(dag_vars), dag_vars)

# Blacklist: outcomes cannot cause needs
# Each row is (from, to) meaning from -> to is FORBIDDEN
blacklist <- expand.grid(
  from = node_idx[outcome_vars],
  to   = node_idx[need_vars]
)
blacklist <- as.matrix(blacklist)

cat("\n=== Expert Constraints ===\n")
cat("Blacklisted edges (outcome -> need):", nrow(blacklist), "\n")
cat("Whitelisted edges: 0\n")
cat("Unconstrained: edges among needs, edges among outcomes, edges from needs to outcomes\n\n")

# Print blacklist for verification
cat("Forbidden edges:\n")
for (r in seq_len(nrow(blacklist))) {
  cat(sprintf(
    "  %s -> %s\n",
    dag_vars[blacklist[r, 1]],
    dag_vars[blacklist[r, 2]]
  ))
}

#
# # ==============================================================================
# # STEP 3: RUN swa-oBN (POINT ESTIMATE)
# # ==============================================================================
#
# cat("\n=== Running swa-oBN (point estimate) ===\n")
#
# fit_point <- swa_obn(
#   y         = y_data,
#   Z         = Z_data,
#   weights   = w_data,
#   search    = "greedy",
#   ic        = "bic",
#   link      = "probit",
#   blacklist = blacklist,
#   whitelist = NULL,
#   nstart    = 5,
#   boot      = NULL,
#   verbose   = TRUE,
#   maxit     = 50
# )
#
# cat("\nEstimated adjacency matrix:\n")
# colnames(fit_point$gam) <- dag_vars
# rownames(fit_point$gam) <- dag_vars
# print(fit_point$gam)
# cat(sprintf("BIC: %.2f\n", fit_point$ic_best))
#
#
# # ==============================================================================
# # STEP 4: RUN swa-oBN (BOOTSTRAP FOR UNCERTAINTY)
# # ==============================================================================
#
# #This will take a while...
#
# source(file.path(here::here(), "code", "2026_03_25_bootstrap_base.R"))
#
# boot_result <- run_chunked_bootstrap(
#   y_data, Z_data, w_data,
#   blacklist = blacklist,
#   B_total = 250, chunk_size = 20,
#   nstart = 5, maxit = 25,
#   save_dir = file.path(here::here(), "code", "boot_chunks_PMAX")
# )
#
#
# load_bootstrap_chunks <- function(save_dir = file.path(here::here(), "code", "boot_chunks_PMAX"), var_names = NULL) {
#
#   chunk_files <- sort(list.files(save_dir, pattern = "^chunk_\\d+\\.rds$",
#                                  full.names = TRUE))
#
#   if (length(chunk_files) == 0) stop("No chunk files found in ", save_dir)
#
#   all_gams <- list()
#   total_time <- 0
#
#   for (f in chunk_files) {
#     chunk <- readRDS(f)
#     all_gams <- c(all_gams, chunk$boot_gams)
#     total_time <- total_time + chunk$elapsed_minutes
#     cat(sprintf("  Loaded %s: resamples %d–%d (%.1f min)\n",
#                 basename(f), chunk$b_start, chunk$b_end, chunk$elapsed_minutes))
#   }
#
#   cat(sprintf("Total: %d resamples, %.1f min compute time\n",
#               length(all_gams), total_time))
#
#   q <- nrow(all_gams[[1]])
#   adj_sum <- matrix(0, q, q)
#   n_valid <- 0
#   for (g in all_gams) {
#     if (!any(is.na(g))) {
#       adj_sum <- adj_sum + g
#       n_valid <- n_valid + 1
#     }
#   }
#
#   edge_probs <- adj_sum / n_valid
#   if (!is.null(var_names)) {
#     colnames(edge_probs) <- var_names
#     rownames(edge_probs) <- var_names
#   }
#
#   return(list(
#     edge_probs = edge_probs,
#     boot_gams = all_gams,
#     B_valid = n_valid,
#     B_failed = length(all_gams) - n_valid
#   ))
# }
#
# boot_result <- load_bootstrap_chunks()
# edge_probs <- boot_result$edge_probs
# gam_point <- (edge_probs > 0.50) * 1
#
# # --- Node roles ---
# need_vars    <- dag_vars[grepl("^n_", dag_vars)]
# outcome_vars <- dag_vars[grepl("^o_", dag_vars)]
#
# my_roles <- list(
#   need_ordinal = need_vars,
#   outcome      = outcome_vars
# )
#
# # --- Short display names ---
# my_short_names <- c(
#   n_HelpCaregive      = "Caregiving",
#   n_HelpChildcare     = "Childcare",
#   n_HelpTransport     = "Transport",
#   n_HelpInternet      = "Internet",
#   n_HelpLegal         = "Legal",
#   n_clust_basics      = "Basic\nNeeds",
#   n_clust_work_edu    = "Work/\nEducation",
#   n_clust_discrim_iso = "Discrim./\nIsolation",
#   o_mh                = "Mental\nHealth",
#   o_wb                = "Well-\nbeing"
# )
#
# # # --- Now plot ---
# # p2 <- plot_dag(
# #   gam         = gam_point,
# #   node_names  = dag_vars,
# #   edge_probs  = edge_probs,
# #   node_roles  = my_roles,
# #   short_names = my_short_names,
# #   threshold   = 0.50,
# #   title       = "swa-oBn: Bootstrap (prob > 0.50)",
# #   layout      = "sugiyama"
# # )
# # print(p2)
#
# p2 <- plot_dag(
#   gam         = gam_point,
#   node_names  = dag_vars,
#   edge_probs  = edge_probs,
#   node_roles  = my_roles,
#   short_names = my_short_names,
#   threshold   = 0.50,
#   layout      = "sugiyama",
#   title       = "Consensus DAG: veteran social needs and mental health outcomes",
#   )
#
# print(p2)
#
# ggsave(file.path(here(), "code", "2026_04_18_PMAX.pdf"), p2, height = 10, width=15, device=cairo_pdf)
#
#
#
#
#
# round(boot_result$edge_probs, 3)
# fit_point$gam
# compute_n_eff(w_data)
# summary(w_data)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#


# ==============================================================================
# STEP 6: RUN unwtd-oBN (BOOTSTRAP FOR UNCERTAINTY)
# ==============================================================================

# This will take a while...

#source(file.path(here::here(), "code", "2026_03_25_bootstrap_base.R"))

#boot_result_unwtd <- run_chunked_bootstrap(
#  y_data, Z_data, rep(1, length(w_data)),
#  blacklist = blacklist,
#  B_total = 250, chunk_size = 20,
#  nstart = 5, maxit = 25,
#  save_dir = file.path(here::here(), "code", "boot_chunks_PMAX_UW")
#)


load_bootstrap_chunks <- function(save_dir = file.path(here::here(), "code", "boot_chunks_PMAX_UW"), var_names = NULL) {
  chunk_files <- sort(list.files(save_dir,
                                 pattern = "^chunk_\\d+\\.rds$",
                                 full.names = TRUE
  ))
  
  if (length(chunk_files) == 0) stop("No chunk files found in ", save_dir)
  
  all_gams <- list()
  total_time <- 0
  
  for (f in chunk_files) {
    chunk <- readRDS(f)
    all_gams <- c(all_gams, chunk$boot_gams)
    total_time <- total_time + chunk$elapsed_minutes
    cat(sprintf(
      "  Loaded %s: resamples %d–%d (%.1f min)\n",
      basename(f), chunk$b_start, chunk$b_end, chunk$elapsed_minutes
    ))
  }
  
  cat(sprintf(
    "Total: %d resamples, %.1f min compute time\n",
    length(all_gams), total_time
  ))
  
  q <- nrow(all_gams[[1]])
  adj_sum <- matrix(0, q, q)
  n_valid <- 0
  for (g in all_gams) {
    if (!any(is.na(g))) {
      adj_sum <- adj_sum + g
      n_valid <- n_valid + 1
    }
  }
  
  edge_probs <- adj_sum / n_valid
  if (!is.null(var_names)) {
    colnames(edge_probs) <- var_names
    rownames(edge_probs) <- var_names
  }
  
  return(list(
    edge_probs = edge_probs,
    boot_gams = all_gams,
    B_valid = n_valid,
    B_failed = length(all_gams) - n_valid
  ))
}

boot_result_uw <- load_bootstrap_chunks(save_dir = file.path(here::here(), "inst/analysis", "boot_chunks_PMAX_UW"))
boot_result_sw <- load_bootstrap_chunks(save_dir = file.path(here::here(), "inst/analysis", "boot_chunks_PMAX"))

edge_probs_uw <- boot_result_uw$edge_probs
edge_probs_sw <- boot_result_sw$edge_probs

gam_point_uw <- (edge_probs_uw > 0.50) * 1
gam_point_sw <- (edge_probs_sw > 0.50) * 1

# --- Node roles ---
need_vars <- dag_vars[grepl("^n_", dag_vars)]
outcome_vars <- dag_vars[grepl("^o_", dag_vars)]

my_roles <- list(
  need_ordinal = need_vars,
  outcome      = outcome_vars
)

# --- Short display names ---
my_short_names <- c(
  n_HelpCaregive      = "Caregiving",
  n_HelpChildcare     = "Childcare",
  n_HelpTransport     = "Transport",
  n_HelpInternet      = "Internet",
  n_HelpLegal         = "Legal",
  n_clust_basics      = "Basic\nNeeds",
  n_clust_work_edu    = "Work/\nEducation",
  n_clust_discrim_iso = "I/L/D",
  o_mh                = "Mental\nHealth",
  o_wb                = "Well-\nbeing"
)

p2_uw <- plot_dag(
  gam         = gam_point_uw,
  node_names  = dag_vars,
  edge_probs  = edge_probs_uw,
  node_roles  = my_roles,
  short_names = my_short_names,
  threshold   = 0.50,
  layout      = "sugiyama",
  title       = "Consensus DAG (unweighted): veteran social needs and mental health outcomes",
)
ggsave(file.path(here(), "inst/analysis", "analysis_consensusDAG_uw.pdf"), p2_uw, height = 10, width = 15, device = cairo_pdf)

p2_sw <- plot_dag(
  gam         = gam_point_sw,
  node_names  = dag_vars,
  edge_probs  = edge_probs_sw,
  node_roles  = my_roles,
  short_names = my_short_names,
  threshold   = 0.50,
  layout      = "sugiyama",
  title       = "Consensus DAG (weighted): veteran social needs and mental health outcomes",
)
ggsave(file.path(here(), "inst/analysis", "analysis_consensusDAG_sw.pdf"), p2_sw, height = 10, width = 15, device = cairo_pdf)


layout_df <- data.frame(
  name = c(
    "n_clust_basics", "n_HelpTransport", "n_HelpInternet",
    "n_HelpLegal", "n_clust_work_edu", "n_HelpCaregive",
    "n_clust_discrim_iso", "o_wb", "o_mh",
    "n_HelpChildcare"
  ),
  x = c(0.0, -1.5, -3.0, -3.0, 1.5, 3.0, 0.0, -1.5, 1.5, 4.5),
  y = c(3.0, 1.5, 1.5, 0.3, 1.5, 1.5, -0.5, -2.5, -2.5, 3.0),
  stringsAsFactors = FALSE
)

node_names <- dag_vars

p_w <- plot_dag_aligned(edge_probs_sw, node_names, layout_df,
                        node_roles  = my_roles,
                        short_names = my_short_names,
                        title       = "(a) Weighted"
)
p_u <- plot_dag_aligned(edge_probs_uw, node_names, layout_df,
                        node_roles  = my_roles,
                        short_names = my_short_names,
                        title       = "(b) Unweighted"
)

fig4 <- (p_w | p_u) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


# Drop Childcare
keep <- node_names != "n_HelpChildcare"
node_names_fig4 <- node_names[keep]

layout_df_fig4 <- data.frame(
  name = c(
    "n_clust_basics", "n_HelpTransport", "n_HelpInternet",
    "n_HelpLegal", "n_clust_work_edu", "n_HelpCaregive",
    "n_clust_discrim_iso", "o_wb", "o_mh"
  ),
  x = c(0.0, -1.5, -3.0, -3.0, 1.5, 3.0, 0.0, -1.5, 1.5),
  y = c(3.0, 1.5, 1.5, 0.3, 1.5, 1.5, -0.5, -2.5, -2.5),
  stringsAsFactors = FALSE
)

# Optional: subset the matrices for cleanliness
edge_probs_w_fig4 <- edge_probs_sw[keep, keep]
edge_probs_u_fig4 <- edge_probs_uw[keep, keep]


p_w <- plot_dag_aligned(edge_probs_w_fig4, node_names_fig4, layout_df_fig4,
                        node_roles  = my_roles,
                        short_names = my_short_names,
                        title       = "(a) Weighted"
)
p_u <- plot_dag_aligned(edge_probs_u_fig4, node_names_fig4, layout_df_fig4,
                        node_roles  = my_roles,
                        short_names = my_short_names,
                        title       = "(b) Unweighted"
)

fig4 <- (p_w | p_u) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
factor <- 1.25
ggsave(file.path(here(), "inst/analysis", "analysis_consensusDAG_both.pdf"), fig4,
       height = factor*9, width = factor*21, device = cairo_pdf
)
