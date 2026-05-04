# ==============================================================================
# 2026_04_27_effects.R
#
# Interventional effect estimation on the recovered swa-oBN.
# Computes do-calculus quantities by reusing the existing bootstrap chunks
# (no re-running of structure search needed):
#
#   - Total effect:      P(MH | do(BN = b))
#   - Controlled direct: P(MH | do(BN = b), do(D/I = d_ref))
#   - Implied indirect:  Total - CDE  (loosely, the mediated portion)
#
# Identification (g-formula on the recovered DAG, with Z as adjustment set):
#
#   P(MH | do(BN=b)) =
#     sum_{di, z} P(MH | BN=b, D/I=di, pa(MH)\{BN,DI}, Z=z) *
#                 P(D/I | BN=b, pa(D/I)\BN,        Z=z) *
#                 P(Z=z)
#
# Each conditional is fit by survey::svyolr on each bootstrap resample,
# with parents read from boot_gams[[b]]. Marginal of Z is the HT estimate.
#
# REQUIREMENTS (assumed already in workspace from the main analysis script):
#   y_data, Z_data, w_data, dag_vars, blacklist
#   boot_result_sw  (loaded by load_bootstrap_chunks for the weighted run)
#
# RUNTIME: a few minutes on a laptop for B = 250.
# ==============================================================================

suppressPackageStartupMessages({
  library(survey)
  library(MASS)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# Configuration: which intervention question?
# ------------------------------------------------------------------------------
# Treatment node, mediator node, and outcome node (as column names of y_data).
# Levels: 1 = no need, 2 = needed and got, 3 = needed and didn't get.
# We frame the contrast as "worst" (3) vs "best" (1) on BN.

TX_NODE  <- "n_clust_basics"
MED_NODE <- "n_clust_discrim_iso"
OUT_NODE <- "o_mh"           # also rerun with "o_wb" if desired

TX_HIGH <- 3                 # "needed support but did not get it"
TX_LOW  <- 1                 # "no support needed"
MED_REF <- 1                 # reference level of mediator for CDE

# ------------------------------------------------------------------------------
# Helper: extract parent set from an adjacency matrix
# Convention used in this codebase: gam[i, j] = 1 means j -> i
# ------------------------------------------------------------------------------
parents_of <- function(node_name, gam, var_names) {
  i <- match(node_name, var_names)
  if (is.na(i)) stop("Node not found: ", node_name)
  par_idx <- which(gam[i, ] == 1)
  var_names[par_idx]
}

# ------------------------------------------------------------------------------
# Helper: fit a survey-weighted ordinal regression with graceful fallback.
# Returns a function predict_probs(newdata) that returns a matrix [n x L]
# of predicted category probabilities.
# ------------------------------------------------------------------------------
fit_ordinal_model <- function(response_name, parent_names, covar_names,
                              df_full, weights_vec) {
  # Build formula: Y ~ pa(Y) + Z
  rhs_vars <- c(parent_names, covar_names)
  if (length(rhs_vars) == 0) rhs_vars <- "1"
  form <- as.formula(paste(response_name, "~", paste(rhs_vars, collapse = " + ")))
  
  df_full$.w <- weights_vec
  des <- tryCatch(
    svydesign(ids = ~1, weights = ~.w, data = df_full),
    error = function(e) NULL
  )
  
  fit <- NULL
  if (!is.null(des)) {
    fit <- tryCatch(
      survey::svyolr(form, design = des, method = "probit"),
      error = function(e) NULL
    )
  }
  if (is.null(fit)) {
    # Fallback: weighted polr (point estimates equivalent; SEs naive)
    fit <- tryCatch(
      MASS::polr(form, data = df_full, weights = .w, method = "probit"),
      error = function(e) NULL
    )
  }
  if (is.null(fit)) return(NULL)
  
  predict_probs <- function(newdata) {
    suppressWarnings(predict(fit, newdata = newdata, type = "probs"))
  }
  predict_probs
}

# ==============================================================================
# Patched functions for 2026_04_27_effects.R
# Fix: complete-case within each resample to handle 82 rows with race_sex_v2 = NA.
# Drop-in replacements for compute_total_effect_one and compute_cde_one.
# ==============================================================================

compute_total_effect_one <- function(y_b, Z_b, w_b, gam_b, var_names,
                                     tx_node, out_node, med_node,
                                     tx_val, covar_names) {
  
  df_full <- cbind(as.data.frame(y_b), as.data.frame(Z_b))
  
  # --- Complete-case within this resample ---
  # 82 rows (~1.5%) of the analytic sample have race_sex_v2 = NA.
  # MASS::polr in swa_obn handled these silently per-fit; survey::svyolr does not,
  # so predictions returned NA across the board. We complete-case the resample
  # for effect estimation; structure recovery used the original (NA-tolerant) fit.
  keep    <- complete.cases(df_full)
  df_full <- df_full[keep, , drop = FALSE]
  w_b     <- w_b[keep]
  if (nrow(df_full) < 100) return(NULL)  # too small after dropping; skip
  
  pa_med <- parents_of(med_node, gam_b, var_names)
  pa_out <- parents_of(out_node, gam_b, var_names)
  
  # Existence check: TX must be in some ancestor chain to OUT.
  # Direct-parent check is sufficient for our application; loosen if needed.
  if (!(tx_node %in% pa_med) && !(tx_node %in% pa_out)) return(NULL)
  
  fit_med_predict <- fit_ordinal_model(med_node, pa_med, covar_names, df_full, w_b)
  fit_out_predict <- fit_ordinal_model(out_node, pa_out, covar_names, df_full, w_b)
  if (is.null(fit_med_predict) || is.null(fit_out_predict)) return(NULL)
  
  n_b        <- nrow(df_full)
  med_levels <- levels(df_full[[med_node]])
  out_levels <- levels(df_full[[out_node]])
  L_med      <- length(med_levels)
  L_out      <- length(out_levels)
  
  df_do <- df_full
  df_do[[tx_node]] <- factor(rep(tx_val, n_b),
                             levels  = levels(df_full[[tx_node]]),
                             ordered = is.ordered(df_full[[tx_node]]))
  
  P_med <- fit_med_predict(df_do)
  if (is.null(dim(P_med))) P_med <- matrix(P_med, nrow = n_b)
  
  out_probs <- matrix(0, nrow = n_b, ncol = L_out)
  for (k in seq_len(L_med)) {
    df_dok <- df_do
    df_dok[[med_node]] <- factor(rep(med_levels[k], n_b),
                                 levels  = med_levels,
                                 ordered = is.ordered(df_full[[med_node]]))
    P_out_k <- fit_out_predict(df_dok)
    if (is.null(dim(P_out_k))) P_out_k <- matrix(P_out_k, nrow = n_b)
    out_probs <- out_probs + P_med[, k] * P_out_k
  }
  
  # HT-weighted marginalization with na.rm guard (belt-and-suspenders)
  w_norm   <- w_b / sum(w_b)
  marginal <- colSums(w_norm * out_probs, na.rm = TRUE)
  
  # If any column was entirely NA, na.rm = TRUE returns 0 for that column;
  # detect that case and bail out instead of returning a misleading 0.
  any_all_na <- apply(out_probs, 2, function(col) all(is.na(col)))
  if (any(any_all_na)) return(NULL)
  
  names(marginal) <- out_levels
  marginal
}


compute_cde_one <- function(y_b, Z_b, w_b, gam_b, var_names,
                            tx_node, out_node, med_node,
                            tx_val, med_ref_val, covar_names) {
  
  df_full <- cbind(as.data.frame(y_b), as.data.frame(Z_b))
  
  keep    <- complete.cases(df_full)
  df_full <- df_full[keep, , drop = FALSE]
  w_b     <- w_b[keep]
  if (nrow(df_full) < 100) return(NULL)
  
  pa_out <- parents_of(out_node, gam_b, var_names)
  
  fit_out_predict <- fit_ordinal_model(out_node, pa_out, covar_names, df_full, w_b)
  if (is.null(fit_out_predict)) return(NULL)
  
  n_b        <- nrow(df_full)
  out_levels <- levels(df_full[[out_node]])
  
  df_do <- df_full
  df_do[[tx_node]]  <- factor(rep(tx_val, n_b),
                              levels  = levels(df_full[[tx_node]]),
                              ordered = is.ordered(df_full[[tx_node]]))
  df_do[[med_node]] <- factor(rep(med_ref_val, n_b),
                              levels  = levels(df_full[[med_node]]),
                              ordered = is.ordered(df_full[[med_node]]))
  
  P_out <- fit_out_predict(df_do)
  if (is.null(dim(P_out))) P_out <- matrix(P_out, nrow = n_b)
  
  w_norm   <- w_b / sum(w_b)
  marginal <- colSums(w_norm * P_out, na.rm = TRUE)
  
  any_all_na <- apply(P_out, 2, function(col) all(is.na(col)))
  if (any(any_all_na)) return(NULL)
  
  names(marginal) <- out_levels
  marginal
}

# ==============================================================================
# Main loop over bootstrap resamples
# ==============================================================================

run_intervention_bootstrap <- function(y_data, Z_data, w_data, boot_gams,
                                       var_names, covar_names,
                                       tx_node, med_node, out_node,
                                       tx_high = 3, tx_low = 1, med_ref = 1,
                                       seed = 42) {
  
  n        <- nrow(y_data)
  B        <- length(boot_gams)
  w_norm0  <- as.numeric(w_data) * (n / sum(as.numeric(w_data)))
  
  # Storage: total effect under TX_HIGH, TX_LOW, and CDE under TX_HIGH, TX_LOW
  out_levels <- levels(y_data[[out_node]])
  L_out      <- length(out_levels)
  
  total_high <- matrix(NA_real_, nrow = B, ncol = L_out, dimnames = list(NULL, out_levels))
  total_low  <- matrix(NA_real_, nrow = B, ncol = L_out, dimnames = list(NULL, out_levels))
  cde_high   <- matrix(NA_real_, nrow = B, ncol = L_out, dimnames = list(NULL, out_levels))
  cde_low    <- matrix(NA_real_, nrow = B, ncol = L_out, dimnames = list(NULL, out_levels))
  
  # Replay the same RNG sequence used by run_chunked_bootstrap (seed = 42, sequential)
  set.seed(seed)
  
  for (b in seq_len(B)) {
    idx <- sample(n, replace = TRUE)
    y_b <- y_data[idx, , drop = FALSE]
    Z_b <- Z_data[idx, , drop = FALSE]
    w_b <- w_norm0[idx]
    gam_b <- boot_gams[[b]]
    
    if (any(is.na(gam_b))) next  # failed resample, skip
    
    # Total effects
    th <- tryCatch(
      compute_total_effect_one(y_b, Z_b, w_b, gam_b, var_names,
                               tx_node, out_node, med_node,
                               tx_high, covar_names),
      error = function(e) NULL
    )
    tl <- tryCatch(
      compute_total_effect_one(y_b, Z_b, w_b, gam_b, var_names,
                               tx_node, out_node, med_node,
                               tx_low, covar_names),
      error = function(e) NULL
    )
    # CDEs
    ch <- tryCatch(
      compute_cde_one(y_b, Z_b, w_b, gam_b, var_names,
                      tx_node, out_node, med_node,
                      tx_high, med_ref, covar_names),
      error = function(e) NULL
    )
    cl <- tryCatch(
      compute_cde_one(y_b, Z_b, w_b, gam_b, var_names,
                      tx_node, out_node, med_node,
                      tx_low, med_ref, covar_names),
      error = function(e) NULL
    )
    
    if (!is.null(th)) total_high[b, ] <- th
    if (!is.null(tl)) total_low [b, ] <- tl
    if (!is.null(ch)) cde_high  [b, ] <- ch
    if (!is.null(cl)) cde_low   [b, ] <- cl
    
    if (b %% 25 == 0) cat(sprintf("  Effect bootstrap %d / %d\n", b, B))
  }
  
  list(
    total_high = total_high,
    total_low  = total_low,
    cde_high   = cde_high,
    cde_low    = cde_low,
    out_levels = out_levels,
    tx_node    = tx_node,
    med_node   = med_node,
    out_node   = out_node
  )
}

# ==============================================================================
# Summarize: contrasts and 95% CIs
# ==============================================================================

summarize_effects <- function(eff, focus_level = NULL) {
  # focus_level: outcome level to focus on (e.g., the worst MH category).
  # If NULL, summarize the expected level (treating ordinal as numeric 1..L).
  
  out_levels <- eff$out_levels
  L_out      <- length(out_levels)
  
  if (is.null(focus_level)) {
    # Expected level under each interventional distribution (1..L)
    lev_num <- seq_len(L_out)
    th_summary <- eff$total_high %*% lev_num
    tl_summary <- eff$total_low  %*% lev_num
    ch_summary <- eff$cde_high   %*% lev_num
    cl_summary <- eff$cde_low    %*% lev_num
    label <- "E[OUT | do(.)]"
  } else {
    k <- match(focus_level, out_levels)
    if (is.na(k)) stop("focus_level not in outcome levels: ", focus_level)
    th_summary <- eff$total_high[, k]
    tl_summary <- eff$total_low [, k]
    ch_summary <- eff$cde_high  [, k]
    cl_summary <- eff$cde_low   [, k]
    label <- sprintf("P(OUT = %s | do(.))", focus_level)
  }
  
  total_contrast <- as.numeric(th_summary) - as.numeric(tl_summary)
  cde_contrast   <- as.numeric(ch_summary) - as.numeric(cl_summary)
  indirect       <- total_contrast - cde_contrast  # implied through-mediator
  
  fmt <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(c(point = NA, lo = NA, hi = NA))
    c(point = mean(x),
      lo    = quantile(x, 0.025, names = FALSE),
      hi    = quantile(x, 0.975, names = FALSE))
  }
  
  out <- rbind(
    Total_effect    = fmt(total_contrast),
    CDE_at_med_ref  = fmt(cde_contrast),
    Implied_through_mediator = fmt(indirect)
  )
  list(label = label, contrasts = out, raw = list(
    total = total_contrast, cde = cde_contrast, indirect = indirect
  ))
}

# ==============================================================================
# RUN
# ==============================================================================
# Assumes the main analysis script has already been sourced and the workspace
# contains: y_data, Z_data, w_data, dag_vars, boot_result_sw

cat("\n=== Effect estimation: do(", TX_NODE, ") on ", OUT_NODE,
    " via ", MED_NODE, " ===\n", sep = "")

eff <- run_intervention_bootstrap(
  y_data       = y_data,
  Z_data       = Z_data,
  w_data       = w_data,
  boot_gams    = boot_result_sw$boot_gams,
  var_names    = dag_vars,
  covar_names  = colnames(Z_data),
  tx_node      = TX_NODE,
  med_node     = MED_NODE,
  out_node     = OUT_NODE,
  tx_high      = TX_HIGH,
  tx_low       = TX_LOW,
  med_ref      = MED_REF,
  seed         = 42
)

# Save raw bootstrap distributions
saveRDS(eff,
        file = file.path(here::here(), "inst/analysis",
                         sprintf("boot_effects.rds",
                                 TX_NODE, MED_NODE, OUT_NODE)))

# Summary 1: contrast on E[OUT | do(.)]
cat("\n--- Contrast on expected outcome level (TX = 3 vs TX = 1) ---\n")
s_mean <- summarize_effects(eff, focus_level = NULL)
print(round(s_mean$contrasts, 4))

# Summary 2: contrast on P(OUT = worst category)
worst <- tail(eff$out_levels, 1)
cat(sprintf("\n--- Contrast on P(%s = %s | do(.)) ---\n", OUT_NODE, worst))
s_worst <- summarize_effects(eff, focus_level = worst)
print(round(s_worst$contrasts, 4))

# Mediation share (only meaningful if total contrast is sizable)
total_pt <- s_worst$contrasts["Total_effect", "point"]
indir_pt <- s_worst$contrasts["Implied_through_mediator", "point"]
if (!is.na(total_pt) && abs(total_pt) > 1e-3) {
  cat(sprintf("\nImplied mediated share (point estimate): %.1f%%\n",
              100 * indir_pt / total_pt))
}
