# ==============================================================================
# swa-oBN Pre-Analysis Diagnostics
# Run AFTER step1.R (which loads and prepares the CAHPS data)
# Run BEFORE step3.R (which fits the swa-oBN)
#
# Checks:
#   1. Level distributions and identifiability tiers (Assumption 6)
#   2. Survey weight diagnostics (Assumption 8)
#   3. Brant test for proportional odds (Assumption 4)
#   4. Sparse cell diagnostics (practical convergence)
#   5. Composite variable validation
# ==============================================================================

# --- Prerequisites ---
# Assumes step1.R has been sourced and `data` exists with:
#   - ordinal DAG variables (n_* and o_* columns)
#   - WEIGHT column
#   - age_v4, race_sex_v2 covariates
rm(list = ls())
source(file.path(here::here(), "code", "2026_03_25_data.R"))
pacman::p_load(here, tidyverse, brant, MASS)

# If step1 hasn't been sourced yet, uncomment:
# source(file.path(here::here(), "analyses", "finals", "2026_03_11_step1.R"))

# ==============================================================================
# 0. IDENTIFY VARIABLES
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


cat("##################################################################\n")
cat("#  swa-oBN Pre-Analysis Diagnostics                              #\n")
cat("##################################################################\n\n")

# ==============================================================================
# 1. LEVEL DISTRIBUTIONS & IDENTIFIABILITY (Assumption 6)
# ==============================================================================

cat("==============================================================\n")
cat("  1. LEVEL DISTRIBUTIONS & IDENTIFIABILITY TIERS\n")
cat("==============================================================\n\n")

level_summary <- data.frame(
  variable   = character(),
  n_levels   = integer(),
  min_freq   = integer(),
  min_pct    = numeric(),
  level_dist = character(),
  ident_tier = character(),
  stringsAsFactors = FALSE
)

for (v in dag_vars) {
  x <- y_data[[v]]
  tab <- table(x, useNA = "ifany")
  nl <- nlevels(x)
  min_f <- min(tab[names(tab) != "NA"])
  min_p <- round(100 * min_f / sum(!is.na(x)), 2)
  dist_str <- paste(sprintf("%s:%d", names(tab), tab), collapse = "  ")
  
  tier <- if (nl >= 3) {
    "Tier 1: Fully identifiable (L >= 3)"
  } else if (nl == 2) {
    "Tier 2: Binary (direction not identified by ordinal machinery)"
  } else {
    "PROBLEM: Constant or single-level variable"
  }
  
  level_summary <- rbind(level_summary, data.frame(
    variable = v, n_levels = nl, min_freq = min_f,
    min_pct = min_p, level_dist = dist_str, ident_tier = tier,
    stringsAsFactors = FALSE
  ))
}

# Print level distributions
for (i in seq_len(nrow(level_summary))) {
  row <- level_summary[i, ]
  flag <- ""
  if (row$min_freq < 30)  flag <- " ** SPARSE **"
  if (row$min_freq < 10)  flag <- " ** VERY SPARSE — convergence risk **"
  if (row$n_levels < 2)   flag <- " ** DEGENERATE — exclude from DAG **"
  
  cat(sprintf("  %-25s  L=%d  min_cell=%d (%.1f%%)%s\n",
              row$variable, row$n_levels, row$min_freq, row$min_pct, flag))
  cat(sprintf("    %s\n", row$level_dist))
}

cat("\nIdentifiability summary:\n")
tier_tab <- table(level_summary$ident_tier)
for (t in names(tier_tab)) {
  cat(sprintf("  %s: %d variables\n", t, tier_tab[t]))
}

# Flag non-degeneracy concerns
sparse_vars <- level_summary$variable[level_summary$min_freq < 30]
if (length(sparse_vars) > 0) {
  cat(sprintf("\n  WARNING: %d variable(s) have levels with < 30 obs:\n",
              length(sparse_vars)))
  cat(sprintf("    %s\n", paste(sparse_vars, collapse = ", ")))
  cat("  Consider collapsing levels or excluding these variables.\n")
}

# Flag binary variables
binary_vars <- level_summary$variable[level_summary$n_levels == 2]
if (length(binary_vars) > 0) {
  cat(sprintf("\n  NOTE: %d binary variable(s) — edge direction identified only\n",
              length(binary_vars)))
  cat("  by BIC penalty or domain constraints, not ordinal machinery:\n")
  cat(sprintf("    %s\n", paste(binary_vars, collapse = ", ")))
}

cat("\n")


# ==============================================================================
# 2. SURVEY WEIGHT DIAGNOSTICS (Assumption 8)
# ==============================================================================

cat("==============================================================\n")
cat("  2. SURVEY WEIGHT DIAGNOSTICS\n")
cat("==============================================================\n\n")

n <- length(w_data)
n_eff <- sum(w_data)^2 / sum(w_data^2)
deff <- n / n_eff
cv_w <- sd(w_data) / mean(w_data)

# Normalized weights (as the algorithm will use them)
w_norm <- w_data * (n / sum(w_data))

cat(sprintf("  Sample size (n):             %d\n", n))
cat(sprintf("  Effective sample size:       %.1f\n", n_eff))
cat(sprintf("  Design effect (n / n_eff):   %.2f\n", deff))
cat(sprintf("  CV of weights:               %.3f\n", cv_w))
cat(sprintf("\n  Raw weight distribution:\n"))
cat(sprintf("    Min:    %.4f\n", min(w_data)))
cat(sprintf("    Q1:     %.4f\n", quantile(w_data, 0.25)))
cat(sprintf("    Median: %.4f\n", median(w_data)))
cat(sprintf("    Q3:     %.4f\n", quantile(w_data, 0.75)))
cat(sprintf("    Max:    %.4f\n", max(w_data)))
cat(sprintf("    Ratio (max/min): %.1f\n", max(w_data) / min(w_data)))

cat(sprintf("\n  Normalized weight distribution (sum to n = %d):\n", n))
cat(sprintf("    Min:    %.4f\n", min(w_norm)))
cat(sprintf("    Max:    %.4f\n", max(w_norm)))
cat(sprintf("    Ratio:  %.1f\n", max(w_norm) / min(w_norm)))

# Check Assumption 8: bounded weights
ratio <- max(w_data) / min(w_data)
if (ratio > 100) {
  cat("\n  WARNING: Weight ratio > 100. Extreme weights detected.\n")
  cat("  Consider weight trimming (e.g., at 1st/99th percentile).\n")
  cat("  Extreme weights can destabilize polr convergence and inflate\n")
  cat("  the variance of the pseudo-MLE (see Assumption 8).\n")
} else if (ratio > 20) {
  cat("\n  NOTE: Weight ratio moderately large (> 20).\n")
  cat("  Kish n_eff is a reasonable approximation, but node-specific\n")
  cat("  design effects may vary. Report this ratio in the paper.\n")
} else {
  cat("\n  Weights are well-behaved. Kish n_eff is a good approximation.\n")
}

# Weight quantile plot data
cat("\n  Top 10 largest weights (check for dominant observations):\n")
top_w <- sort(w_data, decreasing = TRUE)[1:10]
for (k in seq_along(top_w)) {
  cat(sprintf("    #%d: %.4f  (%.1f%% of total weight)\n",
              k, top_w[k], 100 * top_w[k] / sum(w_data)))
}

cat("\n")


# ==============================================================================
# 3. BRANT TEST FOR PROPORTIONAL ODDS (Assumption 4)
# ==============================================================================

cat("==============================================================\n")
cat("  3. BRANT TEST FOR PROPORTIONAL ODDS\n")
cat("==============================================================\n\n")

cat("  Assumption 4 requires proportional odds for the FULL node-level\n")
cat("  regression: P(X_j <= l | X_{pa(j)}, Z). Since pa(j) is unknown\n")
cat("  before running swa-oBN, we test in two stages:\n\n")
cat("  3a. PRE-ANALYSIS: Each ordinal node regressed on Z only, and on\n")
cat("      each possible single parent + Z. This screens for variables\n")
cat("      where proportional odds is untenable under any parent config.\n\n")
cat("  3b. POST-ANALYSIS (run after step3): Each node regressed on its\n")
cat("      discovered parents + Z. This is the definitive test.\n")
cat("      -> See brant_post_analysis() function defined at end.\n\n")

ordinal_vars <- level_summary$variable[level_summary$n_levels >= 3]

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


# ---- 3a. PRE-ANALYSIS: Z-only and single-parent screening ----

cat("  --- 3a. Z-only regressions ---\n\n")

brant_z_only <- data.frame(
  variable = character(), omnibus_X2 = numeric(),
  omnibus_p = numeric(), conclusion = character(),
  stringsAsFactors = FALSE
)

for (v in ordinal_vars) {
  res <- run_brant(y_data[[v]], Z_data, w_norm, label = paste(v, "~ Z"))
  
  brant_z_only <- rbind(brant_z_only, data.frame(
    variable = v, omnibus_X2 = round(res$omnibus_X2, 2),
    omnibus_p = round(res$omnibus_p, 4), conclusion = res$conclusion,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("  %-25s  X2 = %7.2f  p = %.4f  %s\n",
              v, res$omnibus_X2, res$omnibus_p, res$conclusion))
  
  if (!is.na(res$omnibus_p) && res$omnibus_p < 0.05 && !is.null(res$per_var)) {
    cat("    Per-covariate breakdown:\n")
    for (r in seq_len(nrow(res$per_var))) {
      cat(sprintf("      %-20s  X2 = %6.2f  p = %.4f\n",
                  rownames(res$per_var)[r],
                  res$per_var[r, "X2"], res$per_var[r, "probability"]))
    }
  }
}

n_reject_z <- sum(brant_z_only$omnibus_p < 0.05, na.rm = TRUE)
n_tested_z <- sum(!is.na(brant_z_only$omnibus_p))
cat(sprintf("\n  Z-only summary: %d / %d reject at alpha = 0.05\n\n",
            n_reject_z, n_tested_z))


cat("  --- 3a. Single-parent + Z regressions ---\n\n")
cat("  For each ordinal child node, regress on each possible single\n")
cat("  parent (treated as numeric) + covariates Z. This screens whether\n")
cat("  adding any DAG parent breaks proportional odds.\n\n")

# We'll track: for each child, how many potential parents cause rejection?
brant_single_parent <- data.frame(
  child = character(), parent = character(),
  omnibus_X2 = numeric(), omnibus_p = numeric(),
  conclusion = character(), stringsAsFactors = FALSE
)

for (child_v in ordinal_vars) {
  # Potential parents = all other DAG variables
  potential_parents <- setdiff(dag_vars, child_v)
  
  for (parent_v in potential_parents) {
    # Build predictor set: parent (as numeric) + covariates Z
    pred_df <- Z_data
    pred_df[[parent_v]] <- as.numeric(y_data[[parent_v]])
    
    res <- run_brant(y_data[[child_v]], pred_df, w_norm,
                     label = paste(child_v, "~", parent_v, "+ Z"))
    
    brant_single_parent <- rbind(brant_single_parent, data.frame(
      child = child_v, parent = parent_v,
      omnibus_X2 = round(res$omnibus_X2, 2),
      omnibus_p = round(res$omnibus_p, 4),
      conclusion = res$conclusion,
      stringsAsFactors = FALSE
    ))
  }
}

# Summarize by child: how many single-parent configs reject?
cat("  Per-child summary (single-parent + Z regressions):\n\n")
for (child_v in ordinal_vars) {
  sub <- brant_single_parent[brant_single_parent$child == child_v, ]
  n_tested_sp <- sum(!is.na(sub$omnibus_p))
  n_reject_sp <- sum(sub$omnibus_p < 0.05, na.rm = TRUE)
  
  status <- if (n_reject_sp == 0) {
    "all pass"
  } else if (n_reject_sp < n_tested_sp / 2) {
    "some rejections"
  } else {
    "MAJORITY REJECT"
  }
  
  cat(sprintf("  %-25s  %d / %d reject  [%s]\n",
              child_v, n_reject_sp, n_tested_sp, status))
  
  # Print the rejecting parents if any
  if (n_reject_sp > 0 && n_reject_sp <= 5) {
    rejects <- sub[!is.na(sub$omnibus_p) & sub$omnibus_p < 0.05, ]
    for (r in seq_len(nrow(rejects))) {
      cat(sprintf("    -> parent = %-20s  X2 = %6.2f  p = %.4f\n",
                  rejects$parent[r], rejects$omnibus_X2[r], rejects$omnibus_p[r]))
    }
  }
}

n_reject_sp_total <- sum(brant_single_parent$omnibus_p < 0.05, na.rm = TRUE)
n_tested_sp_total <- sum(!is.na(brant_single_parent$omnibus_p))
cat(sprintf("\n  Single-parent overall: %d / %d regressions reject at alpha = 0.05\n",
            n_reject_sp_total, n_tested_sp_total))

# Combined interpretation
n_reject <- n_reject_z  # use Z-only count for summary section later
n_tested <- n_tested_z

cat("\n  INTERPRETATION:\n")
cat("  - Z-only rejections indicate baseline PO violations from demographics.\n")
cat("  - Single-parent rejections indicate which parent-child relationships\n")
cat("    are problematic for the ordinal regression specification.\n")
cat("  - Brant rejections are common in large samples (n > 1000) due to\n")
cat("    high power. The key question is whether violations are large enough\n")
cat("    to compromise the functional-form asymmetry driving Theorem 1.\n")
cat("  - The DEFINITIVE test runs after the DAG is estimated (see 3b below).\n")

if (n_reject_z == 0 && n_reject_sp_total == 0) {
  cat("\n  All pre-analysis Brant tests pass. Assumption 4 is well-supported.\n")
}


# ---- 3b. POST-ANALYSIS FUNCTION (to be called after step3) ----
# Define a function that takes the estimated DAG and runs Brant on
# the actual discovered parent sets.

#' Brant test for the estimated DAG (post-analysis)
#'
#' Call this AFTER running swa_obn in step3.R:
#'   brant_post <- brant_post_analysis(fit_point$gam, y_data, Z_data, w_norm, dag_vars)
#'
#' @param gam Estimated adjacency matrix (q x q). gam[i,j] = 1 means j -> i.
#' @param y_data Data frame of ordinal DAG variables
#' @param Z_data Data frame of covariates
#' @param weights Normalized weight vector
#' @param var_names Character vector of DAG variable names
#' @return Data frame of Brant test results per node
brant_post_analysis <- function(gam, y_data, Z_data, weights, var_names) {
  
  cat("\n==============================================================\n")
  cat("  3b. POST-ANALYSIS BRANT TEST (discovered parents + Z)\n")
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
                    pname, res$per_var[r, "X2"], res$per_var[r, "probability"], tag))
      }
    }
  }
  
  n_reject_post <- sum(results$omnibus_p < 0.05, na.rm = TRUE)
  n_tested_post <- sum(!is.na(results$omnibus_p))
  
  cat(sprintf("\n  Post-analysis summary: %d / %d ordinal nodes reject at alpha = 0.05\n",
              n_reject_post, n_tested_post))
  
  # Check if rejections are driven by parents vs covariates
  cat("\n  KEY QUESTION: Are violations driven by parent variables or\n")
  cat("  covariates? If only covariates violate PO, the ordinal\n")
  cat("  identifiability argument (which operates through parent-child\n")
  cat("  regression structure) may still be sound. If parent variables\n")
  cat("  drive the violation, that directly threatens Assumption 4.\n\n")
  
  return(results)
}

cat("\n")


# ==============================================================================
# 4. SPARSE CELL DIAGNOSTICS (Practical Convergence)
# ==============================================================================

cat("==============================================================\n")
cat("  4. SPARSE CELL DIAGNOSTICS\n")
cat("==============================================================\n\n")

cat("  Checking pairwise cross-tabulations for empty/sparse cells.\n")
cat("  Empty cells can cause polr to fail during hill-climbing.\n\n")

sparse_pairs <- data.frame(
  var1 = character(), var2 = character(),
  n_empty_cells = integer(), min_cell = integer(),
  total_cells = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_along(dag_vars)) {
  for (j in seq_along(dag_vars)) {
    if (i >= j) next
    v1 <- dag_vars[i]; v2 <- dag_vars[j]
    tab <- table(y_data[[v1]], y_data[[v2]])
    n_empty <- sum(tab == 0)
    min_cell <- min(tab)
    total_cells <- prod(dim(tab))
    
    if (n_empty > 0 || min_cell < 10) {
      sparse_pairs <- rbind(sparse_pairs, data.frame(
        var1 = v1, var2 = v2,
        n_empty_cells = n_empty, min_cell = min_cell,
        total_cells = total_cells,
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(sparse_pairs) > 0) {
  cat(sprintf("  %d variable pairs have empty or sparse cells:\n\n", nrow(sparse_pairs)))
  # Sort by severity
  sparse_pairs <- sparse_pairs[order(-sparse_pairs$n_empty_cells, sparse_pairs$min_cell), ]
  for (r in seq_len(min(nrow(sparse_pairs), 20))) {
    row <- sparse_pairs[r, ]
    flag <- if (row$n_empty_cells > 0) "EMPTY CELLS" else "sparse"
    cat(sprintf("    %-22s x %-22s  empty=%d  min=%d  (%d cells)  [%s]\n",
                row$var1, row$var2, row$n_empty_cells, row$min_cell,
                row$total_cells, flag))
  }
  if (nrow(sparse_pairs) > 20) {
    cat(sprintf("    ... and %d more pairs\n", nrow(sparse_pairs) - 20))
  }
  cat("\n  NOTE: Empty cells may cause polr convergence failures during\n")
  cat("  hill-climbing. The cached_node_score function returns Inf for\n")
  cat("  failed fits, which effectively prevents the algorithm from\n")
  cat("  adding edges involving sparse parent configurations. This is\n")
  cat("  conservative (biases toward sparser graphs).\n")
} else {
  cat("  No empty or very sparse cells detected. Good.\n")
}

cat("\n")


# ==============================================================================
# 5. COMPOSITE VARIABLE VALIDATION
# ==============================================================================

cat("==============================================================\n")
cat("  5. COMPOSITE VARIABLE VALIDATION\n")
cat("==============================================================\n\n")

cat("  The following DAG variables are composites created by rounding\n")
cat("  the mean of constituent items. Checking whether rounding\n")
cat("  preserves a meaningful ordinal structure.\n\n")

# Identify composites (those with 'clust' in the name from step1)
composite_vars <- dag_vars[str_detect(dag_vars, "clust")]
non_composite_vars <- setdiff(dag_vars, composite_vars)

if (length(composite_vars) > 0) {
  for (v in composite_vars) {
    x <- y_data[[v]]
    tab <- table(x)
    nl <- nlevels(x)
    entropy <- -sum((tab/sum(tab)) * log(tab/sum(tab) + 1e-10))
    max_entropy <- log(nl)
    evenness <- entropy / max_entropy  # 1 = perfectly even, 0 = degenerate
    
    cat(sprintf("  %-25s  L=%d  evenness=%.2f\n", v, nl, evenness))
    cat(sprintf("    Distribution: %s\n",
                paste(sprintf("%s:%.1f%%", names(tab), 100*tab/sum(tab)),
                      collapse = "  ")))
    
    if (evenness < 0.5) {
      cat("    WARNING: Highly uneven distribution. Rounding may have\n")
      cat("    collapsed meaningful variation. Consider alternative\n")
      cat("    aggregation (sum score, latent variable).\n")
    }
    if (nl < 3) {
      cat("    WARNING: Composite has only 2 levels — ordinal identifiability\n")
      cat("    does not apply. Consider keeping constituent items separate.\n")
    }
    cat("\n")
  }
} else {
  cat("  No composite variables detected.\n\n")
}

cat("  Non-composite DAG variables (original survey items):\n")
for (v in non_composite_vars) {
  x <- y_data[[v]]
  tab <- table(x)
  cat(sprintf("    %-25s  L=%d  %s\n", v, nlevels(x),
              paste(sprintf("%s:%.1f%%", names(tab), 100*tab/sum(tab)),
                    collapse = "  ")))
}


cat("\n")

# ==============================================================================
# 6. SUMMARY AND RECOMMENDATIONS
# ==============================================================================

cat("==============================================================\n")
cat("  6. SUMMARY AND RECOMMENDATIONS\n")
cat("==============================================================\n\n")

# Gather flags
flags <- character(0)

# Identifiability
n_binary <- sum(level_summary$n_levels == 2)
n_ordinal <- sum(level_summary$n_levels >= 3)
if (n_binary > 0) {
  flags <- c(flags, sprintf(
    "- %d/%d DAG variables are binary (direction not identified by ordinal machinery)",
    n_binary, nrow(level_summary)))
}

# Weights
if (ratio > 100) {
  flags <- c(flags, "- Extreme weight ratio (> 100:1). Consider trimming.")
} else if (ratio > 20) {
  flags <- c(flags, sprintf("- Moderate weight ratio (%.0f:1). Report n_eff = %.0f in paper.", ratio, n_eff))
}

# Brant (pre-analysis: Z-only)
if (n_reject_z > 0) {
  flags <- c(flags, sprintf(
    "- %d/%d ordinal variables reject Brant (Z-only). Report as limitation.",
    n_reject_z, n_tested_z))
}
# Brant (pre-analysis: single-parent + Z)
if (n_reject_sp_total > 0) {
  pct_reject_sp <- round(100 * n_reject_sp_total / n_tested_sp_total, 0)
  flags <- c(flags, sprintf(
    "- %d/%d single-parent regressions reject Brant (%d%%). Run post-analysis after step3.",
    n_reject_sp_total, n_tested_sp_total, pct_reject_sp))
}

# Sparsity
if (nrow(sparse_pairs) > 0) {
  n_empty_pairs <- sum(sparse_pairs$n_empty_cells > 0)
  flags <- c(flags, sprintf("- %d variable pairs have empty cells. May cause polr failures.", n_empty_pairs))
}

# Composites
n_composite_binary <- sum(level_summary$n_levels[level_summary$variable %in% composite_vars] == 2)
if (n_composite_binary > 0) {
  flags <- c(flags, sprintf(
    "- %d composite variable(s) collapsed to binary. Consider alternative aggregation.", n_composite_binary))
}

if (length(flags) == 0) {
  cat("  All diagnostics pass. Data are appropriate for swa-oBN.\n")
} else {
  cat("  Issues to address before running swa-oBN:\n\n")
  for (f in flags) cat(sprintf("  %s\n", f))
}

cat(sprintf("\n  Variables entering the DAG: %d (%d ordinal, %d binary)\n",
            nrow(level_summary), n_ordinal, n_binary))
cat(sprintf("  Effective sample size: %.0f (design effect = %.2f)\n", n_eff, deff))
cat(sprintf("  Covariates in Z matrix: %s\n", paste(covar_vars, collapse = ", ")))

cat("\n  NEXT STEPS:\n")
cat("  1. If diagnostics are acceptable, proceed to step3.R to run swa-oBN.\n")
cat("  2. After step3 completes, run the post-analysis Brant test:\n")
cat("       brant_post <- brant_post_analysis(\n")
cat("         fit_point$gam, y_data, Z_data, w_norm, dag_vars)\n")
cat("     This is the DEFINITIVE check of Assumption 4 on the discovered DAG.\n")

cat("\n##################################################################\n")
cat("#  End of Pre-Analysis Diagnostics                               #\n")
cat("##################################################################\n")