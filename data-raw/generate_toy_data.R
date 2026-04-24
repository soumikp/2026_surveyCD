# Generate synthetic toy data for swaoBN examples/vignettes
# This script uses the simulation DGP to create an ordinal dataset
# that loosely resembles the structure of the SHEP data (but is purely synthetic).

set.seed(2026)

# Provide path to the simulation utility file to use the data generator
sim_utils_path <- system.file("simulations", "laptop", "sim_utils.R", package = "swaoBN")
if (sim_utils_path == "") {
  # For local execution before package installation
  sim_utils_path <- file.path("inst", "simulations", "laptop", "sim_utils.R")
}

source(sim_utils_path)

# Simulate a DAG for 13 nodes (similar to the 13 need variables)
q <- 13
# Let's create a sparse Erdos-Renyi DAG
gam_true <- generate_random_dag(q, prob = 0.2)

# Generate ordinal data (n = 500, mostly 3 levels)
L <- rep(3, q)
L[c(1, 4, 7)] <- 2 # Some binary variables
pop <- generate_obn_data(n = 500, gam = gam_true, sigma = 1.0, L = L)

# Create some informative sampling weights (simulating survey design)
svy <- apply_survey_sampling(
  pop_data = pop$data,
  design = "moderate",
  oversample_vars = c(1, q),
  target_n = 300
)

# Extract final dataset
toy_data <- svy$sample_data
toy_weights <- svy$weights

# Add weights as a column for convenience in the example, or store separately.
# Usually weights are passed separately. Let's return a list.
swaoBN_toy <- list(
  data = toy_data,
  weights = toy_weights,
  true_gam = gam_true
)

# Save to data/
usethis::use_data(swaoBN_toy, overwrite = TRUE)
