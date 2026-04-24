# ==============================================================================
# swa-oBN Post-Analysis Diagnostics
#
# Runs AFTER:
#   - 2026_03_25_data.R        (loads `data`)
#   - 2026_03_25_analysis.R    (bootstrap chunks saved to boot_chunks_PMAX/)
#
# Self-contained: defines run_brant and brant_post_analysis internally.
#
# Checks:
#   Brant test for proportional odds on the bootstrap MPM consensus DAG
#   (each node regressed on its discovered parents + covariates Z).
# ==============================================================================

rm(list = ls())

pacman::p_load(here, tidyverse, brant, MASS)

source(file.path(here::here(), "code", "2026_03_25_data.R"))


# ==============================================================================
# 1. HELPER FUNCTIONS (copied verbatim from preDiagnostics.R)
# ==============================================================================

# --- Helper function to run one Brant test ---
run_brant <- function(response_var, predictor_df, weights, label = "") {
  test_df <- data.frame(Y = response_var)
  test_df <- cbind(test_df, predictor_df)
  
  result <- list(label = label, omnibus_X2 = NA, omnibus_p = NA,
                 conclusion = "NOT TESTED", per_var = NULL)
  
  tryCatch({
    mod <- MASS::polr(Y ~ ., data = test_df, weights = weights,
                      method = "probit", Hess = TRUE)
    bt <- brant::brant(mod)
    result$omnibus_X2 <- bt[1, "X2"]
    result$omnibus_p  <- bt[1, "probability"]
    result$conclusion <- if (result$omnibus_p < 0.001) {
      "REJECT (p < 0.001)"
    } else if (result$omnibus_p < 0.05) {
      "REJECT (p < 0.05)"
    } else {
      "PASS"
    }
    # Store per-variable breakdown
    if (nrow(bt) > 1) {
      result$per_var <- bt[2:nrow(bt), , drop = FALSE]
    }
  }, error = function(e) {
    result$conclusion <<- paste("ERROR:", conditionMessage(e))
  })
  
  return(result)
}


#' Brant test for the estimated DAG (post-analysis)
#'
#' @param gam Estimated adjacency matrix (q x q). gam[i,j] = 1 means j -> i.
#' @param y_data Data frame of ordinal DAG variables
#' @param Z_data Data frame of covariates
#' @param weights Normalized weight vector
#' @param var_names Character vector of DAG variable names
#' @return Data frame of Brant test results per node
brant_post_analysis <- function(gam, y_data, Z_data, weights, var_names) {
  
  cat("\n==============================================================\n")
  cat("  POST-ANALYSIS BRANT TEST (discovered parents + Z)\n")
  cat("==============================================================\n\n")
  cat("  Testing proportional odds for each ordinal node regressed on\n")
  cat("  its DISCOVERED parent set + covariates Z. This is the definitive\n")
  cat("  check of Assumption 4 for the estimated DAG.\n\n")
  
  q <- ncol(y_data)
  results <- data.frame(
    node = character(), n_parents = integer(), parents = character(),
    omnibus_X2 = numeric(), omnibus_p = numeric(),
    conclusion = character(), stringsAsFactors = FALSE
  )
  
  for (i in seq_len(q)) {
    v <- var_names[i]
    nl <- nlevels(y_data[[v]])
    if (nl < 3) {
      cat(sprintf("  %-25s  SKIP (binary, L=%d)\n", v, nl))
      next
    }
    
    # Find parents from adjacency matrix: gam[i, j] = 1 means j -> i
    parent_idx <- which(gam[i, ] == 1)
    parent_names <- var_names[parent_idx]
    n_pa <- length(parent_idx)
    
    # Build predictor set: parents (as numeric) + covariates Z
    pred_df <- Z_data
    if (n_pa > 0) {
      for (pv in parent_names) {
        pred_df[[pv]] <- as.numeric(y_data[[pv]])
      }
    }
    
    parent_str <- if (n_pa > 0) paste(parent_names, collapse = ", ") else "(none)"
    label <- sprintf("%s ~ %s + Z", v, parent_str)
    
    res <- run_brant(y_data[[v]], pred_df, weights, label = label)
    
    results <- rbind(results, data.frame(
      node = v, n_parents = n_pa, parents = parent_str,
      omnibus_X2 = round(res$omnibus_X2, 2),
      omnibus_p = round(res$omnibus_p, 4),
      conclusion = res$conclusion,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("  %-25s  pa={%s}  X2 = %7.2f  p = %.4f  %s\n",
                v, parent_str, res$omnibus_X2, res$omnibus_p, res$conclusion))
    
    # Per-predictor breakdown for rejections
    if (!is.na(res$omnibus_p) && res$omnibus_p < 0.05 && !is.null(res$per_var)) {
      cat("    Per-predictor breakdown:\n")
      for (r in seq_len(nrow(res$per_var))) {
        pname <- rownames(res$per_var)[r]
        # Flag if the violation comes from a parent (not a covariate)
        is_parent <- pname %in% parent_names
        tag <- if (is_parent) " [PARENT]" else ""
        cat(sprintf("      %-20s  X2 = %6.2f  p = %.4f%s\n",
                    pname, res$per_var[r, "X2"], res$per_var[r, "pr(>X2)"], tag))
      }
    }
  }
  
  n_reject_post <- sum(results$omnibus_p < 0.05, na.rm = TRUE)
  n_tested_post <- sum(!is.na(results$omnibus_p))
  
  cat(sprintf("\n  Post-analysis summary: %d / %d ordinal nodes reject at alpha = 0.05\n",
              n_reject_post, n_tested_post))
  
  cat("\n  KEY QUESTION: Are violations driven by parent variables or\n")
  cat("  covariates? If only covariates violate PO, the ordinal\n")
  cat("  identifiability argument (which operates through parent-child\n")
  cat("  regression structure) may still be sound. If parent variables\n")
  cat("  drive the violation, that directly threatens Assumption 4.\n\n")
  
  return(results)
}


# ==============================================================================
# 2. RECONSTRUCT DATA OBJECTS
# ==============================================================================

need_vars    <- colnames(data)[str_detect(colnames(data), "n_")]
outcome_vars <- colnames(data)[str_detect(colnames(data), "o_")]
dag_vars     <- c(need_vars, outcome_vars)
covar_vars   <- c("age_v4", "race_sex_v2")

y_data <- data[, dag_vars]
for (v in names(y_data)) {
  if (!is.ordered(y_data[[v]])) y_data[[v]] <- as.ordered(y_data[[v]])
}

Z_data <- data[, covar_vars]
Z_data$race_sex_v2 <- as.factor(Z_data$race_sex_v2)
Z_data$age_v4      <- as.factor(Z_data$age_v4)

w_data <- data$WEIGHT
n      <- nrow(y_data)
w_norm <- w_data * (n / sum(w_data))


# ==============================================================================
# 3. LOAD BOOTSTRAP CHUNKS → MPM CONSENSUS DAG
# ==============================================================================

load_bootstrap_chunks <- function(save_dir, var_names = NULL) {
  chunk_files <- sort(list.files(save_dir, pattern = "^chunk_\\d+\\.rds$",
                                 full.names = TRUE))
  if (length(chunk_files) == 0) stop("No chunk files found in ", save_dir)
  
  all_gams <- list()
  for (f in chunk_files) {
    chunk <- readRDS(f)
    all_gams <- c(all_gams, chunk$boot_gams)
  }
  
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
  
  list(edge_probs = edge_probs, B_valid = n_valid,
       B_failed = length(all_gams) - n_valid)
}

boot_result <- load_bootstrap_chunks(
  save_dir  = file.path(here::here(), "code", "boot_chunks_PMAX"),
  var_names = dag_vars
)

edge_probs <- boot_result$edge_probs
gam_point  <- (edge_probs > 0.50) * 1

cat(sprintf("Bootstrap loaded: B_valid = %d, B_failed = %d\n\n",
            boot_result$B_valid, boot_result$B_failed))


# ==============================================================================
# 4. POST-ANALYSIS BRANT — MPM CONSENSUS DAG
# ==============================================================================

brant_post_mpm <- brant_post_analysis(
  gam       = gam_point,
  y_data    = y_data,
  Z_data    = Z_data,
  weights   = w_norm,
  var_names = dag_vars
)


# ==============================================================================
# 5. SAVE
# ==============================================================================

out_dir <- file.path(here::here(), "code", "postDiagnostics_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

saveRDS(list(
  brant_mpm  = brant_post_mpm,
  edge_probs = edge_probs,
  gam_mpm    = gam_point,
  B_valid    = boot_result$B_valid,
  timestamp  = Sys.time()
), file = file.path(out_dir, "post_diagnostics.rds"))

cat(sprintf("\nSaved to %s\n", file.path(out_dir, "post_diagnostics.rds")))