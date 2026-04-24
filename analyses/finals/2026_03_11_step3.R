rm(list = ls())
source(file.path(here::here(), "analyses", "finals", "2026_03_11_step1.R"))
source(file.path(here::here(), "analyses", "finals", "2026_03_11_step2.R"))


need_vars    <- colnames(data)[str_detect(colnames(data), "n_")]
outcome_vars <- colnames(data)[str_detect(colnames(data), "o_")]
dag_vars     <- c(need_vars, outcome_vars)  # needs first, then outcomes
covar_vars   <- c("age_v4", "race_sex_v2")

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
Z_data$race_sex_v2 <- as.factor(Z_data$race_sex_v2)  # ensure factor for dummies
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
  cat(sprintf("  %-25s %d levels  %s\n",
              v, nlevels(y_data[[v]]),
              ifelse(nlevels(y_data[[v]]) > 2, "[ordinal - identifiable]",
                     "[binary - limited]")))
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
cat("Unconstrained: edges among needs, edges among outcomes,\n")
cat("               edges from needs to outcomes\n\n")

# Print blacklist for verification
cat("Forbidden edges:\n")
for (r in seq_len(nrow(blacklist))) {
  cat(sprintf("  %s -> %s\n",
              dag_vars[blacklist[r, 1]],
              dag_vars[blacklist[r, 2]]))
}


# ==============================================================================
# STEP 3: RUN swa-oBN (POINT ESTIMATE)
# ==============================================================================

cat("\n=== Running swa-oBN (point estimate) ===\n")

fit_point <- swa_obn(
  y         = y_data,
  Z         = Z_data,
  weights   = w_data,
  search    = "greedy",
  ic        = "bic",
  link      = "probit",
  blacklist = blacklist,
  whitelist = NULL,
  nstart    = 10,
  boot      = NULL,
  verbose   = TRUE,
  maxit     = 100
)

cat("\nEstimated adjacency matrix:\n")
colnames(fit_point$gam) <- dag_vars
rownames(fit_point$gam) <- dag_vars
print(fit_point$gam)
cat(sprintf("BIC: %.2f\n", fit_point$ic_best))


# ==============================================================================
# STEP 4: RUN swa-oBN (BOOTSTRAP FOR UNCERTAINTY)
# ==============================================================================

#This will take a while...

fit_boot <- swa_obn(
  y         = y_data,
  Z         = Z_data,
  weights   = w_data,
  search    = "greedy",
  ic        = "bic",
  link      = "probit",
  blacklist = blacklist,
  whitelist = NULL,
  nstart    = 5,       # fewer restarts per bootstrap for speed
  boot      = 100,
  verbose   = TRUE,
  maxit     = 25
)

cat("\nEdge inclusion probabilities:\n")
colnames(fit_boot$edge_probs) <- dag_vars
rownames(fit_boot$edge_probs) <- dag_vars
print(round(fit_boot$edge_probs, 3))


# ==============================================================================
# STEP 5: RESULTS TABLE
# ==============================================================================

# Extract all edges with inclusion probability > 0.5
extract_edges <- function(gam_point, edge_probs, var_names) {
  q <- nrow(gam_point)
  edges <- data.frame(
    from      = character(),
    to        = character(),
    present   = logical(),
    prob      = numeric(),
    edge_type = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(q)) {
    for (j in seq_len(q)) {
      if (i == j) next
      if (edge_probs[i, j] < 0.10) next  # skip very rare edges
      
      # Classify edge type
      from_name <- var_names[j]
      to_name   <- var_names[i]
      from_is_need <- from_name %in% need_vars
      to_is_need   <- to_name %in% need_vars
      
      if (from_is_need && to_is_need) {
        etype <- "need -> need"
      } else if (from_is_need && !to_is_need) {
        etype <- "need -> outcome"
      } else if (!from_is_need && !to_is_need) {
        etype <- "outcome -> outcome"
      } else {
        etype <- "other"
      }
      
      edges <- rbind(edges, data.frame(
        from      = from_name,
        to        = to_name,
        present   = as.logical(gam_point[i, j]),
        prob      = edge_probs[i, j],
        edge_type = etype,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  edges <- edges[order(-edges$prob), ]
  return(edges)
}

edge_table <- extract_edges(fit_boot$gam, fit_boot$edge_probs, dag_vars)

cat("\n=== Significant Edges (inclusion prob > 0.50) ===\n")
sig_edges <- edge_table[edge_table$prob > 0.50, ]
print(sig_edges, row.names = FALSE)

cat("\n=== Suggestive Edges (0.30 < prob < 0.50) ===\n")
sug_edges <- edge_table[edge_table$prob > 0.30 & edge_table$prob <= 0.50, ]
print(sug_edges, row.names = FALSE)


# ==============================================================================
# STEP 6: IDENTIFIABILITY TIER ANNOTATION
# ==============================================================================

annotate_identifiability <- function(edge_table, y_data) {
  nlevels_map <- sapply(y_data, nlevels)
  
  edge_table$from_levels <- nlevels_map[edge_table$from]
  edge_table$to_levels   <- nlevels_map[edge_table$to]
  
  edge_table$ident_tier <- ifelse(
    edge_table$from_levels > 2 & edge_table$to_levels > 2,
    "Tier 1: Fully identifiable (both ordinal)",
    ifelse(
      edge_table$from_levels > 2 | edge_table$to_levels > 2,
      "Tier 2: Asymmetric parameterization (ordinal x binary)",
      "Tier 3: Not identifiable by ordinal machinery (both binary)"
    )
  )
  
  return(edge_table)
}

edge_table <- annotate_identifiability(edge_table, y_data)

cat("\n=== Edge Identifiability Summary ===\n")
sig_annotated <- edge_table[edge_table$prob > 0.50, ]
for (tier in unique(sig_annotated$ident_tier)) {
  cat(sprintf("\n%s:\n", tier))
  sub <- sig_annotated[sig_annotated$ident_tier == tier, ]
  for (r in seq_len(nrow(sub))) {
    cat(sprintf("  %s -> %s  (prob = %.3f)\n",
                sub$from[r], sub$to[r], sub$prob[r]))
  }
}


# ==============================================================================
# STEP 7: VISUALIZATION
# ==============================================================================

plot_dag <- function(gam, edge_probs, var_names, threshold = 0.50) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is required for plotting")
  }
  
  q <- nrow(gam)
  el <- data.frame(from = character(), to = character(),
                   weight = numeric(), stringsAsFactors = FALSE)
  
  for (i in seq_len(q)) {
    for (j in seq_len(q)) {
      if (i == j) next
      if (edge_probs[i, j] >= threshold) {
        el <- rbind(el, data.frame(
          from = var_names[j], to = var_names[i],
          weight = edge_probs[i, j], stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  if (nrow(el) == 0) {
    cat("No edges above threshold", threshold, "\n")
    return(invisible(NULL))
  }
  
  g <- igraph::graph_from_data_frame(el, directed = TRUE, vertices = var_names)
  
  # Color nodes by role
  node_colors <- ifelse(var_names %in% outcome_vars, "#E74C3C",
                        ifelse(var_names %in% need_vars[6:8], "#3498DB",
                               "#95A5A6"))
  names(node_colors) <- var_names
  
  igraph::V(g)$color <- node_colors[igraph::V(g)$name]
  igraph::E(g)$width <- el$weight * 4
  igraph::E(g)$label <- round(el$weight, 2)
  igraph::E(g)$arrow.size <- 0.5
  
  plot(g,
       vertex.size = 25,
       vertex.label.cex = 0.7,
       vertex.label.color = "black",
       edge.curved = 0.2,
       edge.label.cex = 0.6,
       main = "swa-oBN Estimated Causal Network")
  
  legend("bottomleft",
         legend = c("Outcomes (o_mh, o_wb)",
                    "Needs - Ordinal (clust)",
                    "Needs - Binary (Help)"),
         fill = c("#E74C3C", "#3498DB", "#95A5A6"),
         cex = 0.8, bty = "n")
}

cat("\n=== Plotting DAG (edges with prob > 0.50) ===\n")
plot_dag(fit_boot$gam, fit_boot$edge_probs, dag_vars, threshold = 0.50)


# ==============================================================================
# STEP 8: COMPARISON WITH NAIVE MODELS
# ==============================================================================

cat("\n=== Comparison: swa-oBN vs Naive Models ===\n\n")

# Model 1: No covariates, no weights, no constraints (fully naive)
cat("Running naive model (no Z, no weights, no constraints)...\n")
fit_naive <- swa_obn(
  y = y_data, Z = NULL, weights = NULL,
  search = "greedy", ic = "bic", link = "probit",
  blacklist = NULL, whitelist = NULL,
  nstart = 10, boot = NULL, verbose = FALSE, maxit = 100
)

# Model 2: With Z only
cat("Running Z-adjusted model (no weights, no constraints)...\n")
fit_z_only <- swa_obn(
  y = y_data, Z = Z_data, weights = NULL,
  search = "greedy", ic = "bic", link = "probit",
  blacklist = NULL, whitelist = NULL,
  nstart = 10, boot = NULL, verbose = FALSE, maxit = 100
)

# Model 3: With Z + weights
cat("Running Z + weights model (no constraints)...\n")
fit_zw <- swa_obn(
  y = y_data, Z = Z_data, weights = w_data,
  search = "greedy", ic = "bic", link = "probit",
  blacklist = NULL, whitelist = NULL,
  nstart = 10, boot = NULL, verbose = FALSE, maxit = 100
)

# Model 4: Full swa-oBN (Z + weights + constraints) — already fit above
cat("\nEdge counts by model:\n")
cat(sprintf("  Naive:              %d edges\n", sum(fit_naive$gam)))
cat(sprintf("  Z-adjusted:         %d edges\n", sum(fit_z_only$gam)))
cat(sprintf("  Z + weighted:       %d edges\n", sum(fit_zw$gam)))
cat(sprintf("  swa-oBN (full):     %d edges\n", sum(fit_point$gam)))

cat("\nBIC by model:\n")
cat(sprintf("  Naive:              %.2f\n", fit_naive$ic_best))
cat(sprintf("  Z-adjusted:         %.2f\n", fit_z_only$ic_best))
cat(sprintf("  Z + weighted:       %.2f\n", fit_zw$ic_best))
cat(sprintf("  swa-oBN (full):     %.2f\n", fit_point$ic_best))

# Count constraint violations in unconstrained models
count_violations <- function(gam, blacklist, var_names) {
  violations <- 0
  for (r in seq_len(nrow(blacklist))) {
    from <- blacklist[r, 1]; to <- blacklist[r, 2]
    if (gam[to, from] == 1) violations <- violations + 1
  }
  return(violations)
}

cat("\nBlacklist violations (outcome -> need edges) in unconstrained models:\n")
cat(sprintf("  Naive:      %d / %d\n",
            count_violations(fit_naive$gam, blacklist, dag_vars), nrow(blacklist)))
cat(sprintf("  Z-adjusted: %d / %d\n",
            count_violations(fit_z_only$gam, blacklist, dag_vars), nrow(blacklist)))
cat(sprintf("  Z+weighted: %d / %d\n",
            count_violations(fit_zw$gam, blacklist, dag_vars), nrow(blacklist)))

save.image(file.path(here::here(), "analyses", "finals", "2026-03-13-output.Rds"))
