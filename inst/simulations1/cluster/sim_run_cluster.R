###############################################################################
# sim_run_cluster.R
#
# Three modes:
#   count   — Print grid size and recommended array settings, then exit
#   <int>   — SLURM array mode: run one chunk of tasks
#   local N — Local parallel mode with N cores
#
# Usage:
#   Rscript sim_run_cluster.R count              # Check before submitting
#   sbatch sim_run_cluster.sh                    # Submit array job
#   Rscript sim_run_cluster.R 1                  # Test: run first chunk
#   Rscript sim_run_cluster.R local 8            # Local test with 8 cores
#   Rscript sim_run_cluster.R local 4 test       # Quick test (5 scenarios)
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

BASE_DIR <- "/ihome/spurkayastha/soumik/2026_surveyCD/simulations"
setwd(BASE_DIR)

# --- Load packages ---
library(MASS)
library(igraph)
library(gRbase)
library(bnlearn)
library(dplyr)

# Conditionally load Study 3 packages
STUDY3_AVAILABLE <- tryCatch({
  library(pcalg)
  library(survey)
  TRUE
}, error = function(e) {
  message("NOTE: pcalg or survey not installed. Study 3 tasks will be skipped.")
  FALSE
})

# --- Source functions ---
source(file.path(BASE_DIR, "sim_utils.R"))
source(file.path(BASE_DIR, "2026_03_11_step2.R"))

# --- Output directory ---
out_dir <- file.path(BASE_DIR, "sim_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


###############################################################################
# CONSTANTS & TRUE DAGs
###############################################################################

# Study 1 settings
Q_NODES_1     <- 8
L_LEVELS      <- 3
EDGE_PROB     <- 0.25
OVERSAMPLE_1  <- c(1, Q_NODES_1)
NSTART_1      <- 5       # hill-climbing restarts (Study 1)
MAXIT_1       <- 25      # max iterations per restart (Study 1)
N_POP_1       <- 20000
SIGMA_1A      <- 1.5
TARGET_N_1A   <- 1000
SIGMA_1B      <- 1.5
SIGMA_1C_TN   <- 1000    # target_n for Study 1c

# Study 2 settings
NSTART_2      <- 10
MAXIT_2       <- 100
N_POP_2       <- 20000
SIGMA_2       <- 1.5
TARGET_N_2    <- 1000
DELTA_STR_2B  <- 1.2     # fixed delta_strength for Study 2b
CONFOUNDED_2  <- c(1, 2, 5, 7, 8)

# Study 3 settings
N_POP_3       <- 20000

# Replication counts
N_REPS_12     <- 100      # Studies 1 and 2
N_REPS_3      <- 100     # Study 3

# Generate TRUE DAGs (deterministic)
set.seed(42)
TRUE_DAG_1 <- generate_random_dag(Q_NODES_1, prob = EDGE_PROB)

TRUE_DAG_2 <- matrix(0L, 8, 8)
TRUE_DAG_2[2, 1] <- 1   # 1 -> 2
TRUE_DAG_2[3, 2] <- 1   # 2 -> 3
TRUE_DAG_2[4, 3] <- 1   # 3 -> 4
TRUE_DAG_2[5, 1] <- 1   # 1 -> 5
TRUE_DAG_2[6, 5] <- 1   # 5 -> 6
TRUE_DAG_2[7, 3] <- 1   # 3 -> 7
TRUE_DAG_2[8, 6] <- 1   # 6 -> 8
TRUE_DAG_2[8, 4] <- 1   # 4 -> 8


###############################################################################
# TASK TABLE
###############################################################################

build_task_table <- function() {
  
  scenario_id <- 0L
  all_tasks <- list()
  
  # --- Study 1a: vary design strength ---
  for (des in c("none", "mild", "moderate", "extreme")) {
    scenario_id <- scenario_id + 1L
    all_tasks[[length(all_tasks) + 1]] <- data.frame(
      study = "1a", scenario_id = scenario_id,
      rep_id = seq_len(N_REPS_12),
      design = des, sigma = SIGMA_1A, target_n = TARGET_N_1A,
      delta_strength = NA_real_, algorithm = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  # --- Study 1b: vary sample size ---
  for (tn in c(300, 500, 1000, 2000, 4000)) {
    scenario_id <- scenario_id + 1L
    all_tasks[[length(all_tasks) + 1]] <- data.frame(
      study = "1b", scenario_id = scenario_id,
      rep_id = seq_len(N_REPS_12),
      design = "moderate", sigma = SIGMA_1B, target_n = tn,
      delta_strength = NA_real_, algorithm = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  # --- Study 1c: vary signal strength ---
  for (sig in c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0)) {
    scenario_id <- scenario_id + 1L
    all_tasks[[length(all_tasks) + 1]] <- data.frame(
      study = "1c", scenario_id = scenario_id,
      rep_id = seq_len(N_REPS_12),
      design = "moderate", sigma = sig, target_n = SIGMA_1C_TN,
      delta_strength = NA_real_, algorithm = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  # --- Study 2a: vary confounding strength ---
  for (ds in c(0, 0.5, 1.0, 1.5, 2.0)) {
    scenario_id <- scenario_id + 1L
    all_tasks[[length(all_tasks) + 1]] <- data.frame(
      study = "2a", scenario_id = scenario_id,
      rep_id = seq_len(N_REPS_12),
      design = "moderate", sigma = NA_real_, target_n = TARGET_N_2,
      delta_strength = ds, algorithm = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  # --- Study 2b: vary design strength (fixed confounding) ---
  for (des in c("none", "mild", "moderate", "extreme")) {
    scenario_id <- scenario_id + 1L
    all_tasks[[length(all_tasks) + 1]] <- data.frame(
      study = "2b", scenario_id = scenario_id,
      rep_id = seq_len(N_REPS_12),
      design = des, sigma = NA_real_, target_n = TARGET_N_2,
      delta_strength = DELTA_STR_2B, algorithm = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  # --- Study 3a/3b/3c: algorithms x designs ---
  algo_study_map <- c(PC = "3a", FCI = "3b", LiNGAM = "3c")
  for (algo in c("PC", "FCI", "LiNGAM")) {
    for (des in c("none", "mild", "moderate", "extreme")) {
      scenario_id <- scenario_id + 1L
      all_tasks[[length(all_tasks) + 1]] <- data.frame(
        study = algo_study_map[algo], scenario_id = scenario_id,
        rep_id = seq_len(N_REPS_3),
        design = des, sigma = NA_real_, target_n = NA_integer_,
        delta_strength = NA_real_, algorithm = algo,
        stringsAsFactors = FALSE
      )
    }
  }
  
  tasks <- bind_rows(all_tasks)
  tasks$task_id <- seq_len(nrow(tasks))
  return(tasks)
}

tasks   <- build_task_table()
n_tasks <- nrow(tasks)


###############################################################################
# TASK DISPATCHER
###############################################################################

run_one_task <- function(params) {
  
  # --- Study 1 ---
  if (params$study %in% c("1a", "1b", "1c")) {
    n_pop <- if (params$study == "1b") {
      max(N_POP_1, params$target_n * 10)
    } else {
      N_POP_1
    }
    
    res <- run_single_replication(
      gam = TRUE_DAG_1, N = n_pop, sigma = params$sigma, L = L_LEVELS,
      design = params$design, target_n = params$target_n,
      nstart = NSTART_1, maxit = MAXIT_1,
      oversample_vars = OVERSAMPLE_1,
      swa_obn_fn = swa_obn
    )
    
    # Attach task metadata
    res$study          <- params$study
    res$scenario_id    <- params$scenario_id
    res$rep_id         <- params$rep_id
    res$design         <- params$design
    res$sigma          <- params$sigma
    res$target_n_param <- params$target_n
    return(res)
  }
  
  # --- Study 2 ---
  if (params$study %in% c("2a", "2b")) {
    res <- run_study2_replication(
      gam = TRUE_DAG_2, N = N_POP_2, sigma = SIGMA_2, L = L_LEVELS,
      design = params$design, target_n = TARGET_N_2,
      confounded_nodes = CONFOUNDED_2,
      delta_strength = params$delta_strength,
      nstart = NSTART_2, maxit = MAXIT_2,
      swa_obn_fn = swa_obn
    )
    
    res$study          <- params$study
    res$scenario_id    <- params$scenario_id
    res$rep_id         <- params$rep_id
    res$design         <- params$design
    res$delta_strength <- params$delta_strength
    return(res)
  }
  
  # --- Study 3 ---
  if (params$study %in% c("3a", "3b", "3c")) {
    if (!STUDY3_AVAILABLE) {
      return(data.frame(
        study = params$study, scenario_id = params$scenario_id,
        rep_id = params$rep_id, design = params$design,
        algorithm = params$algorithm, method = "ERROR",
        error_msg = "pcalg/survey not installed",
        stringsAsFactors = FALSE
      ))
    }
    
    res <- run_study3_replication(
      algorithm = params$algorithm,
      design    = params$design,
      N_POP     = N_POP_3
    )
    
    if (is.null(res)) {
      return(data.frame(
        study = params$study, scenario_id = params$scenario_id,
        rep_id = params$rep_id, design = params$design,
        algorithm = params$algorithm, method = "ERROR",
        error_msg = "algorithm returned NULL",
        stringsAsFactors = FALSE
      ))
    }
    
    res$study       <- params$study
    res$scenario_id <- params$scenario_id
    res$rep_id      <- params$rep_id
    res$design      <- params$design
    res$algorithm   <- params$algorithm
    return(res)
  }
  
  stop("Unknown study: ", params$study)
}


###############################################################################
# MODE: COUNT
###############################################################################

if (length(args) >= 1 && args[1] == "count") {
  
  n_scenarios <- length(unique(tasks$scenario_id))
  
  cat("============================================\n")
  cat("swa-oBN Simulation Grid Summary\n")
  cat("============================================\n")
  cat(sprintf("  Total tasks:         %d\n", n_tasks))
  cat(sprintf("  Unique scenarios:    %d\n", n_scenarios))
  cat("\n  Breakdown by study:\n")
  
  study_summary <- tasks %>%
    group_by(study) %>%
    summarize(
      n_tasks    = n(),
      n_scenarios = n_distinct(scenario_id),
      n_reps     = max(rep_id),
      .groups    = "drop"
    )
  
  for (i in seq_len(nrow(study_summary))) {
    s <- study_summary[i, ]
    # Estimate time per task (seconds)
    est_sec <- switch(substr(s$study, 1, 1),
                      "1" = 10,
                      "2" = 40,
                      "3" = 5
    )
    cat(sprintf("    Study %s: %4d tasks (%2d scenarios x %3d reps) ~%ds/task\n",
                s$study, s$n_tasks, s$n_scenarios, s$n_reps, est_sec))
  }
  
  cat("\n  Study 1: nstart =", NSTART_1, " maxit =", MAXIT_1, "\n")
  cat("  Study 2: nstart =", NSTART_2, " maxit =", MAXIT_2, " (heavier)\n")
  cat("  Study 3: pcalg-based (fast per task)\n")
  
  cat("\n  TRUE DAG 1 edges:", sum(TRUE_DAG_1), "\n")
  cat("  TRUE DAG 2 edges:", sum(TRUE_DAG_2), "\n")
  
  cat("\n  Array size recommendations:\n")
  for (n_arr in c(200, 300, 500)) {
    chunk <- ceiling(n_tasks / n_arr)
    # Rough estimate: average ~15s per task
    est_min <- chunk * 15 / 60
    cat(sprintf("    --array=1-%-4d  ->  %2d tasks/job  ->  ~%.0f min/job\n",
                n_arr, chunk, est_min))
  }
  
  cat("\n  Packages needed:\n")
  cat("    Studies 1-2: MASS, igraph, gRbase, bnlearn, dplyr\n")
  cat("    Study 3:     + pcalg, survey\n")
  
  cat("\n  Test first:  Rscript sim_run_cluster.R 1\n")
  cat("  Submit:      sbatch sim_run_cluster.sh\n")
  cat("  Aggregate:   Rscript sim_summarize.R\n")
  cat("============================================\n")
  quit(save = "no")
}


###############################################################################
# MODE: SLURM ARRAY
###############################################################################

if (length(args) >= 1 && args[1] != "local") {
  
  array_id <- as.integer(args[1])
  n_array  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_MAX",
                                    unset = Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "300")))
  
  chunk_size <- ceiling(n_tasks / n_array)
  task_start <- (array_id - 1) * chunk_size + 1
  task_end   <- min(array_id * chunk_size, n_tasks)
  
  if (task_start > n_tasks) {
    cat(sprintf("Array %d: no tasks (total %d, %d slots). Exiting.\n",
                array_id, n_tasks, n_array))
    quit(save = "no")
  }
  
  my_tasks <- tasks[task_start:task_end, ]
  cat(sprintf("Array %d/%d: tasks %d-%d (%d tasks)\n",
              array_id, n_array, task_start, task_end, nrow(my_tasks)))
  cat(sprintf("Studies in this chunk: %s\n",
              paste(unique(my_tasks$study), collapse = ", ")))
  
  results  <- vector("list", nrow(my_tasks))
  t_start  <- proc.time()
  
  for (i in seq_len(nrow(my_tasks))) {
    params <- my_tasks[i, ]
    
    # Reproducible seed: scenario_id * 10000 + rep_id
    set.seed(params$scenario_id * 10000 + params$rep_id)
    
    t0 <- proc.time()
    results[[i]] <- tryCatch(
      run_one_task(params),
      error = function(e) {
        data.frame(
          study       = params$study,
          scenario_id = params$scenario_id,
          rep_id      = params$rep_id,
          design      = params$design,
          method      = "ERROR",
          error_msg   = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
    elapsed <- (proc.time() - t0)["elapsed"]
    
    if (i %% 5 == 0 || i == nrow(my_tasks)) {
      total_elapsed <- (proc.time() - t_start)["elapsed"]
      rate <- i / total_elapsed * 60
      eta  <- (nrow(my_tasks) - i) / max(rate, 0.01)
      cat(sprintf("  [%d/%d] Study %s | %.1fs | %.1f tasks/min | ETA %.1f min\n",
                  i, nrow(my_tasks), params$study, elapsed, rate, eta))
    }
  }
  
  out <- bind_rows(results)
  out_file <- file.path(out_dir, sprintf("results_%04d.rds", array_id))
  saveRDS(out, out_file)
  
  total_time <- (proc.time() - t_start)["elapsed"]
  n_errors <- sum(grepl("ERROR", out$method, fixed = TRUE) |
                    grepl("ERROR", as.character(out$error_msg), fixed = TRUE),
                  na.rm = TRUE)
  # Count errors more robustly
  if ("method" %in% names(out)) {
    n_errors <- sum(out$method == "ERROR", na.rm = TRUE)
  } else {
    n_errors <- 0
  }
  
  cat(sprintf("\nSaved %d rows -> %s (%.1f min, %d errors)\n",
              nrow(out), out_file, total_time / 60, n_errors))
}


###############################################################################
# MODE: LOCAL PARALLEL
###############################################################################

if (length(args) >= 1 && args[1] == "local") {
  
  n_cores <- if (length(args) >= 2) as.integer(args[2]) else
    max(parallel::detectCores() - 1, 1)
  
  if (length(args) >= 3 && args[3] == "test") {
    tasks <- tasks %>% filter(scenario_id <= 5)
    cat(sprintf("TEST MODE: %d tasks only\n", nrow(tasks)))
  }
  
  cat(sprintf("Local: %d tasks on %d cores\n", nrow(tasks), n_cores))
  
  scenario_groups <- split(tasks, tasks$scenario_id)
  
  run_batch <- function(batch) {
    results <- vector("list", nrow(batch))
    for (i in seq_len(nrow(batch))) {
      params <- batch[i, ]
      set.seed(params$scenario_id * 10000 + params$rep_id)
      results[[i]] <- tryCatch(
        run_one_task(params),
        error = function(e) {
          data.frame(study = params$study, scenario_id = params$scenario_id,
                     rep_id = params$rep_id, method = "ERROR",
                     error_msg = conditionMessage(e),
                     stringsAsFactors = FALSE)
        }
      )
    }
    bind_rows(results)
  }
  
  library(parallel)
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {
    library(MASS); library(igraph); library(gRbase); library(bnlearn)
    library(dplyr)
    tryCatch({ library(pcalg); library(survey) }, error = function(e) NULL)
  })
  
  # Export everything the workers need
  clusterExport(cl, c(
    # Functions from sim_utils.R
    "generate_random_dag", "generate_obn_data", "apply_survey_sampling",
    "compute_shd", "compute_edge_metrics", "run_single_replication",
    "generate_confounded_population", "run_study2_replication",
    "sample_with_design", "svy_ci_test", "svy_lingam",
    "skeleton_metrics_3node", "run_study3_replication",
    "DESIGN_PARAMS_S3",
    # Functions from 2026_03_11_step2.R
    "swa_obn", "greedy_search", "multi_start_search", "exhaustive_search",
    "node_score", "cached_node_score", "compute_n_eff",
    "is_dag_after_add", "is_dag_after_rev", "is_allowed",
    # Constants & globals
    "TRUE_DAG_1", "TRUE_DAG_2",
    "Q_NODES_1", "L_LEVELS", "EDGE_PROB", "OVERSAMPLE_1",
    "NSTART_1", "MAXIT_1", "NSTART_2", "MAXIT_2",
    "N_POP_1", "N_POP_2", "N_POP_3",
    "SIGMA_1A", "SIGMA_1B", "SIGMA_1C_TN", "TARGET_N_1A",
    "SIGMA_2", "TARGET_N_2", "DELTA_STR_2B", "CONFOUNDED_2",
    "STUDY3_AVAILABLE",
    # Dispatcher
    "run_one_task"
  ))
  
  t0 <- proc.time()
  all_results <- parLapplyLB(cl, scenario_groups, run_batch)
  stopCluster(cl)
  
  out <- bind_rows(all_results)
  out_file <- file.path(out_dir, "results_all.rds")
  saveRDS(out, out_file)
  cat(sprintf("Done: %d rows in %.1f min -> %s\n",
              nrow(out), (proc.time() - t0)["elapsed"] / 60, out_file))
}


###############################################################################
# NO ARGS
###############################################################################

if (length(args) == 0) {
  cat("Usage:\n")
  cat("  Rscript sim_run_cluster.R count              # Grid info\n")
  cat("  Rscript sim_run_cluster.R <array_id>         # SLURM task\n")
  cat("  Rscript sim_run_cluster.R local <cores>      # Local parallel\n")
  cat("  Rscript sim_run_cluster.R local 4 test       # Quick test\n")
}