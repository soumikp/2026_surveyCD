###############################################################################
# verify_application_v2.R
#
# Same as v1 but with the correct adjacency convention: M[i, j] = 1 means
# there is an edge from j to i (i.e., M[child, parent]). Therefore
#   ep[from, to] is read as ep[child=from_var, parent=to_var]  -- WRONG
#   ep[to, from] is read as ep[child=to_var, parent=from_var]  -- RIGHT
# So to retrieve P(A -> B), use ep[B, A].
###############################################################################

suppressPackageStartupMessages({
    library(dplyr)
})

ANALYSIS_DIR <- "/Users/soumikp/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research/2026_surveyCD/inst/analysis"
setwd(ANALYSIS_DIR)

dag_vars <- c(
    "n_HelpCaregive", "n_HelpChildcare", "n_HelpTransport", "n_HelpInternet",
    "n_HelpLegal", "n_clust_basics", "n_clust_work_edu", "n_clust_discrim_iso",
    "o_mh", "o_wb"
)

load_bootstrap <- function(save_dir, var_names = NULL) {
    chunk_files <- sort(list.files(save_dir,
        pattern = "^chunk_\\d+\\.rds$",
        full.names = TRUE
    ))
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
    list(edge_probs = edge_probs, B_valid = n_valid, B_total = length(all_gams))
}

# Helper: P(from -> to) given the M[child, parent] convention
P_edge <- function(ep, from, to) ep[to, from]

# Helper: out-edges (where 'from' is the parent) at threshold
out_edges <- function(ep, from, thr = 0.50) {
    # ep[, from] = column 'from' = P(from -> i) for each row i
    which(ep[, from] >= thr)
}

# Helper: in-edges (where 'to' is the child) at threshold
in_edges <- function(ep, to, thr = 0.50) {
    # ep[to, ] = row 'to' = P(j -> to) for each column j
    which(ep[to, ] >= thr)
}

chk <- function(label, computed, claimed) {
    match_flag <- if (isTRUE(all.equal(computed, claimed,
        tolerance = 0.005
    ))) {
        "OK   "
    } else {
        "MISMATCH"
    }
    cat(sprintf(
        "  [%s] %-50s computed=%-12s claimed=%s\n",
        match_flag, label,
        format(computed, nsmall = 2), format(claimed, nsmall = 2)
    ))
}

cat("\n=========================================================\n")
cat(" VERIFICATION v2: §7 APPLICATION + APPENDIX F\n")
cat(" Convention: ep[child, parent] = P(parent -> child)\n")
cat("=========================================================\n")

boot_sw <- load_bootstrap("boot_chunks_PMAX", var_names = dag_vars)
boot_uw <- load_bootstrap("boot_chunks_PMAX_UW", var_names = dag_vars)
ep_sw <- boot_sw$edge_probs
ep_uw <- boot_uw$edge_probs

cat(sprintf(
    "\n  Weighted:   B_valid=%d / B_total=%d\n",
    boot_sw$B_valid, boot_sw$B_total
))
cat(sprintf(
    "  Unweighted: B_valid=%d / B_total=%d\n",
    boot_uw$B_valid, boot_uw$B_total
))

###############################################################################
# Build directed edge tables under correct convention:
# row index = child, column index = parent
###############################################################################
edges_at <- function(ep, thr = 0.50) {
    idx <- which(ep >= thr, arr.ind = TRUE)
    data.frame(
        from = colnames(ep)[idx[, 2]], # parent
        to   = rownames(ep)[idx[, 1]], # child
        p    = ep[idx]
    ) %>% arrange(desc(p))
}

mpm_sw <- edges_at(ep_sw, 0.50)
mpm_uw <- edges_at(ep_uw, 0.50)

cat("\n--- Consensus DAG (weighted, p >= 0.50) ---\n")
chk("Total consensus edges (>=0.50)", nrow(mpm_sw), 14)
chk("Robust edges (>=0.75)", sum(mpm_sw$p >= 0.75), 8)
chk("Suggestive edges [0.50,0.75)", sum(mpm_sw$p >= 0.50 & mpm_sw$p < 0.75), 6)

cat("\n  Weighted edges (sorted):\n")
print(mpm_sw, row.names = FALSE, digits = 3)

cat("\n--- Consensus DAG (unweighted, p >= 0.50) ---\n")
cat(sprintf("  Total unweighted edges: %d (manuscript claims 18)\n", nrow(mpm_uw)))
cat("\n  Unweighted edges (sorted):\n")
print(mpm_uw, row.names = FALSE, digits = 3)

###############################################################################
# Childcare
###############################################################################
cat("\n--- Childcare connectivity ---\n")
cc_in <- ep_sw["n_HelpChildcare", ]
cc_out <- ep_sw[, "n_HelpChildcare"]
chk("Childcare out-degree at p>=0.50", sum(cc_out >= 0.50), 0)
chk("Childcare in-degree  at p>=0.50", sum(cc_in >= 0.50), 0)
cat(sprintf("  Max childcare in-prob:  %.3f\n", max(cc_in)))
cat(sprintf("  Max childcare out-prob: %.3f\n", max(cc_out)))

###############################################################################
# MH<->WB
###############################################################################
cat("\n--- MH<->WB orientation ---\n")
chk("p(WB -> MH) weighted", round(P_edge(ep_sw, "o_wb", "o_mh"), 2), 0.54)
chk("p(MH -> WB) weighted", round(P_edge(ep_sw, "o_mh", "o_wb"), 2), 0.46)

###############################################################################
# Blacklist
###############################################################################
cat("\n--- Blacklist verification ---\n")
need_vars <- dag_vars[grepl("^n_", dag_vars)]
outcome_vars <- dag_vars[grepl("^o_", dag_vars)]
# Blacklist: outcome -> need; in M[child, parent], that's M[need, outcome]
blk_sw <- ep_sw[need_vars, outcome_vars]
blk_uw <- ep_uw[need_vars, outcome_vars]
chk("Blacklist cells count", length(blk_sw), 16)
chk("Max p in blacklist cells (weighted)", max(blk_sw), 0)
chk("Max p in blacklist cells (unweighted)", max(blk_uw), 0)

###############################################################################
# Weighted vs unweighted edge sets
###############################################################################
cat("\n--- Weighted vs Unweighted ---\n")
sw_set <- paste(mpm_sw$from, mpm_sw$to, sep = "->")
uw_set <- paste(mpm_uw$from, mpm_uw$to, sep = "->")
shared <- intersect(sw_set, uw_set)
only_sw <- setdiff(sw_set, uw_set)
only_uw <- setdiff(uw_set, sw_set)

reversed <- character(0)
for (e in only_sw) {
    parts <- strsplit(e, "->")[[1]]
    rev_e <- paste(parts[2], parts[1], sep = "->")
    if (rev_e %in% only_uw) reversed <- c(reversed, e)
}

chk("Total weighted edges", length(sw_set), 14)
chk("Total unweighted edges", length(uw_set), 18)
chk("Shared", length(shared), 11)
chk("Only-weighted", length(only_sw), 3)
chk("Only-unweighted", length(only_uw), 7)
chk("Reversed", length(reversed), 1)

cat("\n  Reversed: ")
print(reversed)
cat("\n  Only-weighted: ")
print(only_sw)
cat("\n  Only-unweighted: ")
print(only_uw)

###############################################################################
# Discrim out-degree
###############################################################################
cat("\n--- Discrim/Iso out-degree ---\n")
disc_out_sw <- sum(ep_sw[, "n_clust_discrim_iso"] >= 0.50)
disc_out_uw <- sum(ep_uw[, "n_clust_discrim_iso"] >= 0.50)
chk("Discrim out-degree (weighted)", disc_out_sw, 4)
chk("Discrim out-degree (unweighted)", disc_out_uw, 8)

###############################################################################
# App F three severing comparisons
###############################################################################
cat("\n--- App F: severing comparisons ---\n")
cat(sprintf(
    "  basic -> discrim:    weighted=%.2f  unweighted=%.2f\n",
    P_edge(ep_sw, "n_clust_basics", "n_clust_discrim_iso"),
    P_edge(ep_uw, "n_clust_basics", "n_clust_discrim_iso")
))
cat(sprintf(
    "  discrim -> basic:    weighted=%.2f  unweighted=%.2f\n",
    P_edge(ep_sw, "n_clust_discrim_iso", "n_clust_basics"),
    P_edge(ep_uw, "n_clust_discrim_iso", "n_clust_basics")
))
cat(sprintf(
    "  basic -> WB:         weighted=%.2f  unweighted=%.2f\n",
    P_edge(ep_sw, "n_clust_basics", "o_wb"),
    P_edge(ep_uw, "n_clust_basics", "o_wb")
))
cat(sprintf(
    "  internet -> discrim: weighted=%.2f  unweighted=%.2f\n",
    P_edge(ep_sw, "n_HelpInternet", "n_clust_discrim_iso"),
    P_edge(ep_uw, "n_HelpInternet", "n_clust_discrim_iso")
))

###############################################################################
# App E: cells flipping across 0.50 (excluding diagonal and blacklist)
###############################################################################
cat("\n--- App E: flips across 0.50 ---\n")
flip <- ((ep_sw >= 0.50) != (ep_uw >= 0.50))
diag(flip) <- FALSE
flip[need_vars, outcome_vars] <- FALSE # exclude blacklist cells
chk("Flipping cells", sum(flip), 10)

flip_idx <- which(flip, arr.ind = TRUE)
flip_df <- data.frame(
    from = colnames(ep_sw)[flip_idx[, 2]], # parent
    to   = rownames(ep_sw)[flip_idx[, 1]], # child
    p_sw = ep_sw[flip],
    p_uw = ep_uw[flip]
)
print(flip_df, row.names = FALSE, digits = 3)

cat("\n=========================================================\n")
cat(" END VERIFICATION v2\n")
cat("=========================================================\n\n")
