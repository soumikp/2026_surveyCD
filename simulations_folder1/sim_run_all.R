# ==============================================================================
# sim_run_all.R
# Master script: runs all simulation studies in sequence
#
# BEFORE RUNNING:
#   1. Place these files in the same directory:
#      - sim_utils.R           (data generation, sampling, metrics)
#      - 2026_03_11_step2.R    (swa_obn function)
#      - sim_study1_main.R     (Study 1: swa-oBN vs naive oBN)
#      - sim_study2_covariates.R (Study 2: covariate adjustment)
#      - sim_study3_supplement.R (Study 3: other algorithms)
#
#   2. Install required packages:
#      install.packages(c("MASS", "igraph", "gRbase", "bnlearn",
#                         "ggplot2", "dplyr", "tidyr",
#                         "pcalg", "survey"))
#
# RUNTIME ESTIMATES (on a modern laptop):
#   Study 1: ~4-8 hours  (50 reps x 15 scenarios x ~10s each)
#   Study 2: ~3-6 hours  (50 reps x 9 scenarios x ~10s each)
#   Study 3: ~1-2 hours  (100 reps x 12 scenarios, fast per rep)
#
# To run a subset, source individual study files instead.
#
# OUTPUT:
#   sim_study1_results.RData, sim_study2_results.RData, sim_study3_results.RData
#   sim_study1a_design.pdf, sim_study1b_samplesize.pdf, sim_study1c_signal.pdf
#   sim_study1_orientation.pdf
#   sim_study2a_confounding.pdf, sim_study2b_design_confounding.pdf
#   sim_study3_supplement.pdf, sim_study3_supplement_correct.pdf
# ==============================================================================

t0 <- Sys.time()
cat("========================================\n")
cat("Starting simulation suite:", format(t0), "\n")
cat("========================================\n\n")

# --- Study 1 ---
cat(">>> STUDY 1: swa-oBN vs naive oBN <<<\n")
t1 <- Sys.time()
source("sim_study1_main.R")
cat(sprintf("Study 1 completed in %.1f minutes\n\n", 
            difftime(Sys.time(), t1, units = "mins")))

# --- Study 2 ---
cat(">>> STUDY 2: Covariate adjustment <<<\n")
t2 <- Sys.time()
source("sim_study2_covariates.R")
cat(sprintf("Study 2 completed in %.1f minutes\n\n",
            difftime(Sys.time(), t2, units = "mins")))

# --- Study 3 ---
cat(">>> STUDY 3: Other algorithm families (supplement) <<<\n")
t3 <- Sys.time()
source("sim_study3_supplement.R")
cat(sprintf("Study 3 completed in %.1f minutes\n\n",
            difftime(Sys.time(), t3, units = "mins")))

# --- Done ---
total <- difftime(Sys.time(), t0, units = "hours")
cat("========================================\n")
cat(sprintf("All simulations complete. Total time: %.1f hours\n", total))
cat("========================================\n")
