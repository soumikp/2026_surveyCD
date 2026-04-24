#!/bin/bash
#
# swa-oBN Simulation — Pitt CRC SLURM submission
#
# Usage:
#   Step 1: Check grid size
#           Rscript sim_run_cluster.R count
#
#   Step 2: Submit
#           sbatch sim_run_cluster.sh
#
#   Step 3: Monitor
#           squeue -u $USER
#           sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS
#
#   Step 4: Aggregate
#           Rscript sim_summarize.R
#
###############################################################################

#SBATCH --job-name=swa_obn_sim
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --array=1-300
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-01:00:00
#SBATCH --output=/ihome/spurkayastha/soumik/2026_surveyCD/simulations/log/sim_%A_%a.out
#SBATCH --error=/ihome/spurkayastha/soumik/2026_surveyCD/simulations/log/sim_%A_%a.err
#SBATCH --mail-user=soumik@pitt.edu
#SBATCH --mail-type=END,FAIL

###############################################################################
# Notes:
#
# --array=1-300       300 array tasks. Each handles ~8 tasks from the grid
#                     (2400 total / 300 = 8 tasks per job).
#                     Run `Rscript sim_run_cluster.R count` to verify.
#
# --time=0-01:00:00   1 hour per job. Most will finish in 5-15 min.
#                     Study 2 tasks are heaviest (nstart=10, maxit=100,
#                     4 models per rep). Adjust upward if jobs time out.
#
# --mem=8G            Generous for R + polr + 20k population generation.
#                     Monitor with sacct and reduce to 4G if wasteful.
#
# To run only Studies 1-2 (skip pcalg dependency), filter in the R script
# or submit a smaller array covering only tasks 1-1200.
###############################################################################

# Create log directory if it doesn't exist
mkdir -p /ihome/spurkayastha/soumik/2026_surveyCD/simulations/log

# Load R module (check available versions with: module spider R)
module load gcc
module load r

# Print job info for debugging
echo "============================================"
echo "SLURM Job ID:    ${SLURM_JOB_ID}"
echo "Array Task ID:   ${SLURM_ARRAY_TASK_ID}"
echo "Array Task Max:  ${SLURM_ARRAY_TASK_MAX}"
echo "Node:            ${SLMD_NODENAME}"
echo "Working Dir:     $(pwd)"
echo "Start Time:      $(date)"
echo "============================================"

# Set working directory
cd /ihome/spurkayastha/soumik/2026_surveyCD/simulations

# Run the simulation
Rscript sim_run_cluster.R ${SLURM_ARRAY_TASK_ID}

echo "End Time: $(date)"
echo "Exit Code: $?"