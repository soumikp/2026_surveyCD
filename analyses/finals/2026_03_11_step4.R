# ==============================================================================
# Publication-Quality DAG Plotting
# Requires: ggraph, tidygraph, ggplot2, dplyr
# ==============================================================================

pacman::p_load(ggraph, tidygraph, ggplot2, dplyr, ggparty)


#' Plot a DAG from an adjacency matrix
#'
#' @param gam        Adjacency matrix (q x q). gam[i,j] = 1 means j -> i.
#' @param node_names Character vector of node names (length q)
#' @param edge_probs Optional matrix (q x q) of bootstrap edge inclusion probs.
#'                   If provided, edge width and label reflect probability.
#' @param node_roles Named list: list(need = c(...), outcome = c(...))
#'                   for color-coding. Unassigned nodes get a default color.
#' @param threshold  Minimum edge probability to display (only used if edge_probs given)
#' @param title      Plot title
#' @param node_size  Base node size
#' @param label_size Node label text size
#' @param show_edge_labels Show probability on edges (if edge_probs given)
#' @param layout     ggraph layout: "sugiyama" (layered), "fr", "kk", "stress", "circle"
#' @param short_names Optional named vector for shorter display labels,
#'                    e.g., c(n_HelpCaregive = "Caregiving", o_mh = "Mental\nHealth")
#'
#' @return A ggplot2 object
plot_dag <- function(gam,
                     node_names = NULL,
                     edge_probs = NULL,
                     node_roles = NULL,
                     threshold = 0.0,
                     title = NULL,
                     node_size = 18,
                     label_size = 3.2,
                     show_edge_labels = TRUE,
                     layout = "sugiyama",
                     short_names = NULL) {
  
  q <- nrow(gam)
  if (is.null(node_names)) node_names <- colnames(gam)
  if (is.null(node_names)) node_names <- paste0("X", seq_len(q))
  
  # --- Build edge list ---
  edges <- data.frame(from = character(), to = character(),
                      prob = numeric(), stringsAsFactors = FALSE)
  for (i in seq_len(q)) {
    for (j in seq_len(q)) {
      if (i == j) next
      if (!is.null(edge_probs)) {
        if (edge_probs[i, j] < threshold) next
        edges <- rbind(edges, data.frame(
          from = node_names[j], to = node_names[i],
          prob = edge_probs[i, j], stringsAsFactors = FALSE))
      } else {
        if (gam[i, j] == 0) next
        edges <- rbind(edges, data.frame(
          from = node_names[j], to = node_names[i],
          prob = 1, stringsAsFactors = FALSE))
      }
    }
  }
  
  # --- Build node table ---
  nodes <- data.frame(name = node_names, stringsAsFactors = FALSE)
  
  # Assign roles
  nodes$role <- "Other"
  if (!is.null(node_roles)) {
    if (!is.null(node_roles$need))    nodes$role[nodes$name %in% node_roles$need] <- "Need"
    if (!is.null(node_roles$outcome)) nodes$role[nodes$name %in% node_roles$outcome] <- "Outcome"
    if (!is.null(node_roles$need_ordinal))
      nodes$role[nodes$name %in% node_roles$need_ordinal] <- "Need (ordinal)"
    if (!is.null(node_roles$need_binary))
      nodes$role[nodes$name %in% node_roles$need_binary] <- "Need (binary)"
  }
  
  # Display labels
  nodes$label <- nodes$name
  if (!is.null(short_names)) {
    matched <- match(nodes$name, names(short_names))
    has_short <- !is.na(matched)
    nodes$label[has_short] <- short_names[matched[has_short]]
  }
  
  # --- Build tidygraph ---
  g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  # --- Color palette ---
  role_colors <- c(
    "Outcome"        = "#C0392B",
    "Need"           = "#2980B9",
    "Need (ordinal)" = "#2980B9",
    "Need (binary)"  = "#7F8C8D",
    "Other"          = "#BDC3C7"
  )
  
  role_fills <- c(
    "Outcome"        = "#FADBD8",
    "Need"           = "#D6EAF8",
    "Need (ordinal)" = "#D6EAF8",
    "Need (binary)"  = "#EAECEE",
    "Other"          = "#F2F3F4"
  )
  
  # --- Edge aesthetics ---
  if (!is.null(edge_probs)) {
    edge_aes <- aes(
      width = edges$prob,
      alpha = edges$prob,
      label = ifelse(show_edge_labels, sprintf("%.2f", edges$prob), "")
    )
  } else {
    edge_aes <- aes()
  }
  
  # --- Plot ---
  p <- ggraph(g, layout = layout) +
    
    # Edges
    geom_edge_link(
      edge_aes,
      arrow = arrow(length = unit(4, "mm"), type = "closed"),
      end_cap = circle(node_size * 0.55, "mm"),
      start_cap = circle(node_size * 0.55, "mm"),
      colour = "#2C3E50",
      edge_width = if (is.null(edge_probs)) 0.6 else NULL
    ) +
    
    # Edge labels (probability)
    {if (!is.null(edge_probs) && show_edge_labels)
      geom_edge_label(
        aes(label = sprintf("%.2f", edges$prob)),
        size = 2.3,
        fill = "white",
        label.padding = unit(0.12, "lines"),
        alpha = 1
      )
    } +
    
    # Node outlines
    geom_node_point(
      aes(colour = role),
      size = node_size,
      shape = 21,
      fill = "white",
      stroke = 2
    ) +
    
    # Node fills (lighter)
    geom_node_point(
      aes(fill = role),
      size = node_size - 2,
      shape = 21,
      colour = NA
    ) +
    
    # Node labels
    geom_node_text(
      aes(label = label),
      size = label_size,
      colour = "#2C3E50",
      fontface = "bold"
    ) +
    
    # Scales
    scale_colour_manual(values = role_colors, name = "Node Type") +
    scale_fill_manual(values = role_fills, name = "Node Type") +
    
    {if (!is.null(edge_probs))
      scale_edge_width(range = c(0.3, 2.5), name = "Edge\nProbability",
                       limits = c(0, 1))
    } +
    {if (!is.null(edge_probs))
      scale_edge_alpha(range = c(0.25, 1), guide = "none")
    } +
    
    # Theme
    theme_void(base_size = 12) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                                margin = margin(b = 10)),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#7F8C8D",
                                   margin = margin(b = 10)),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  if (!is.null(title)) {
    edge_count <- nrow(edges)
    p <- p + labs(
      title = title,
      subtitle = sprintf("%d nodes, %d edges", q, edge_count)
    )
  }
  
  return(p)
}


# ==============================================================================
# USAGE EXAMPLE WITH YOUR DATA
# ==============================================================================

# --- Define your variable roles and short names ---
my_roles <- list(
  need_binary  = c(),
  need_ordinal = c("n_clust_basics", "n_clust_work_edu", "n_clust_discrim_iso", "n_HelpCaregive", "n_HelpChildcare", "n_HelpTransport",
                   "n_HelpInternet", "n_HelpLegal"),
  outcome      = c("o_mh", "o_wb")
)

my_short_names <- c(
  n_HelpCaregive      = "Caregiving",
  n_HelpChildcare     = "Childcare",
  n_HelpTransport     = "Transport",
  n_HelpInternet      = "Internet",
  n_HelpLegal         = "Legal",
  n_clust_basics      = "Basic\nNeeds",
  n_clust_work_edu    = "Work/\nEducation",
  n_clust_discrim_iso = "Discrim./\nIsolation",
  o_mh                = "Mental\nHealth",
  o_wb                = "Well-\nbeing"
)


# --- Plot point estimate (no bootstrap) ---
p1 <- plot_dag(
  gam        = fit_point$gam,
  node_names = dag_vars,
  node_roles = my_roles,
  short_names = my_short_names,
  title      = "swa-oBN: Point Estimate",
  layout     = "sugiyama"
)
print(p1)
ggsave("dag_point_estimate.pdf", p1, width = 10, height = 7)


# --- Plot bootstrap version (with edge probabilities) ---
# Only show edges with inclusion prob > 0.50
p2 <- plot_dag(
  gam         = fit_boot$gam,
  node_names  = dag_vars,
  edge_probs  = fit_boot$edge_probs,
  node_roles  = my_roles,
  short_names = my_short_names,
  threshold   = 0.50,
  title       = "swa-oBN: Bootstrap (prob > 0.50)",
  layout      = "sugiyama"
)
print(p2)
ggsave("dag_bootstrap_050.pdf", p2, width = 10, height = 7)


# --- Side-by-side comparison at different thresholds ---
p3 <- plot_dag(
  gam         = fit_boot$gam,
  node_names  = dag_vars,
  edge_probs  = fit_boot$edge_probs,
  node_roles  = my_roles,
  short_names = my_short_names,
  threshold   = 0.30,
  title       = "Bootstrap (prob > 0.30)",
  layout      = "sugiyama"
)

p4 <- plot_dag(
  gam         = fit_boot$gam,
  node_names  = dag_vars,
  edge_probs  = fit_boot$edge_probs,
  node_roles  = my_roles,
  short_names = my_short_names,
  threshold   = 0.70,
  title       = "Bootstrap (prob > 0.70)",
  layout      = "sugiyama"
)

# Combine with patchwork if available
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined <- p3 + p4 + plot_layout(guides = "collect")
  print(combined)
  ggsave("dag_threshold_comparison.pdf", combined, width = 16, height = 7)
}
