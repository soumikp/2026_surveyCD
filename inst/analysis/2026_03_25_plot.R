#' # ==============================================================================
#' # Publication-Quality DAG Plotting
#' # Requires: ggraph, tidygraph, ggplot2, dplyr
#' # ==============================================================================
#' 
#' pacman::p_load(ggraph, tidygraph, ggplot2, dplyr, ggparty)
#' 
#' 
#' #' Plot a DAG from an adjacency matrix
#' #'
#' #' @param gam        Adjacency matrix (q x q). gam[i,j] = 1 means j -> i.
#' #' @param node_names Character vector of node names (length q)
#' #' @param edge_probs Optional matrix (q x q) of bootstrap edge inclusion probs.
#' #'                   If provided, edge width and label reflect probability.
#' #' @param node_roles Named list: list(need = c(...), outcome = c(...))
#' #'                   for color-coding. Unassigned nodes get a default color.
#' #' @param threshold  Minimum edge probability to display (only used if edge_probs given)
#' #' @param title      Plot title
#' #' @param node_size  Base node size
#' #' @param label_size Node label text size
#' #' @param show_edge_labels Show probability on edges (if edge_probs given)
#' #' @param layout     ggraph layout: "sugiyama" (layered), "fr", "kk", "stress", "circle"
#' #' @param short_names Optional named vector for shorter display labels,
#' #'                    e.g., c(n_HelpCaregive = "Caregiving", o_mh = "Mental\nHealth")
#' #'
#' #' @return A ggplot2 object
#' plot_dag <- function(gam,
#'                      node_names = NULL,
#'                      edge_probs = NULL,
#'                      node_roles = NULL,
#'                      threshold = 0.0,
#'                      title = NULL,
#'                      node_size = 30,
#'                      label_size = 5,
#'                      show_edge_labels = TRUE,
#'                      layout = "sugiyama",
#'                      short_names = NULL) {
#'   
#'   q <- nrow(gam)
#'   if (is.null(node_names)) node_names <- colnames(gam)
#'   if (is.null(node_names)) node_names <- paste0("X", seq_len(q))
#'   
#'   # --- Build edge list ---
#'   edges <- data.frame(from = character(), to = character(),
#'                       prob = numeric(), stringsAsFactors = FALSE)
#'   for (i in seq_len(q)) {
#'     for (j in seq_len(q)) {
#'       if (i == j) next
#'       if (!is.null(edge_probs)) {
#'         if (edge_probs[i, j] < threshold) next
#'         edges <- rbind(edges, data.frame(
#'           from = node_names[j], to = node_names[i],
#'           prob = edge_probs[i, j], stringsAsFactors = FALSE))
#'       } else {
#'         if (gam[i, j] == 0) next
#'         edges <- rbind(edges, data.frame(
#'           from = node_names[j], to = node_names[i],
#'           prob = 1, stringsAsFactors = FALSE))
#'       }
#'     }
#'   }
#'   
#'   # --- Remove isolated nodes (no edges in or out) ---
#'   connected_nodes <- unique(c(edges$from, edges$to))
#'   if (length(connected_nodes) < length(node_names)) {
#'     removed <- setdiff(node_names, connected_nodes)
#'     if (length(removed) > 0) {
#'       message(sprintf("Removed %d isolated node(s): %s",
#'                       length(removed), paste(removed, collapse = ", ")))
#'     }
#'     node_names <- node_names[node_names %in% connected_nodes]
#'   }
#'   
#'   # --- Build node table ---
#'   nodes <- data.frame(name = node_names, stringsAsFactors = FALSE)
#'   
#'   # Assign roles
#'   nodes$role <- "Other"
#'   if (!is.null(node_roles)) {
#'     if (!is.null(node_roles$need))    nodes$role[nodes$name %in% node_roles$need] <- "Need"
#'     if (!is.null(node_roles$outcome)) nodes$role[nodes$name %in% node_roles$outcome] <- "Outcome (ordinal)"
#'     if (!is.null(node_roles$need_ordinal))
#'       nodes$role[nodes$name %in% node_roles$need_ordinal] <- "Need (ordinal)"
#'     if (!is.null(node_roles$need_binary))
#'       nodes$role[nodes$name %in% node_roles$need_binary] <- "Need (binary)"
#'   }
#'   
#'   # Display labels
#'   nodes$label <- nodes$name
#'   if (!is.null(short_names)) {
#'     matched <- match(nodes$name, names(short_names))
#'     has_short <- !is.na(matched)
#'     nodes$label[has_short] <- short_names[matched[has_short]]
#'   }
#'   
#'   # --- Build tidygraph ---
#'   g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
#'   
#'   # --- Color palette ---
#'   role_colors <- c(
#'     "Outcome (ordinal)"        = "#C0392B",
#'     "Need (ordinal)"           = "#2980B9",
#'     "Need (ordinal)" = "#2980B9",
#'     "Need (binary)"  = "#7F8C8D",
#'     "Other"          = "#BDC3C7"
#'   )
#'   
#'   role_fills <- c(
#'     "Outcome (ordinal)"        = "#FADBD8",
#'     "Need (ordinal)"           = "#D6EAF8",
#'     "Need (ordinal)" = "#D6EAF8",
#'     "Need (binary)"  = "#EAECEE",
#'     "Other"          = "#F2F3F4"
#'   )
#'   
#'   # --- Plot ---
#'   p <- ggraph(g, layout = layout) +
#'     
#'     # Edges (with optional probability labels baked in)
#'     {if (!is.null(edge_probs) && show_edge_labels)
#'       geom_edge_link(
#'         aes(width = prob, alpha = prob, label = sprintf("%.2f", prob)),
#'         arrow = arrow(length = unit(4, "mm"), type = "closed"),
#'         end_cap = circle(node_size * 0.55, "mm"),
#'         start_cap = circle(node_size * 0.55, "mm"),
#'         colour = "#2C3E50",
#'         angle_calc = "along",
#'         label_dodge = unit(3, "mm"),
#'         label_size = 5
#'       )
#'       else if (!is.null(edge_probs))
#'         geom_edge_link(
#'           aes(width = prob, alpha = prob),
#'           arrow = arrow(length = unit(4, "mm"), type = "closed"),
#'           end_cap = circle(node_size * 0.55, "mm"),
#'           start_cap = circle(node_size * 0.55, "mm"),
#'           colour = "#2C3E50"
#'         )
#'       else
#'         geom_edge_link(
#'           arrow = arrow(length = unit(4, "mm"), type = "closed"),
#'           end_cap = circle(node_size * 0.55, "mm"),
#'           start_cap = circle(node_size * 0.55, "mm"),
#'           colour = "#2C3E50",
#'           edge_width = 0.6
#'         )
#'     } +
#'     
#'     # Node outlines
#'     geom_node_point(
#'       aes(colour = role),
#'       size = node_size,
#'       shape = 21,
#'       fill = "white",
#'       stroke = 2
#'     ) +
#'     
#'     # Node fills (lighter)
#'     geom_node_point(
#'       aes(fill = role),
#'       size = node_size - 2,
#'       shape = 21,
#'       colour = NA
#'     ) +
#'     
#'     # Node labels
#'     geom_node_text(
#'       aes(label = label),
#'       size = label_size,
#'       colour = "#2C3E50",
#'       fontface = "bold"
#'     ) +
#'     
#'     # Scales
#'     scale_colour_manual(values = role_colors, name = "Node Type") +
#'     scale_fill_manual(values = role_fills, name = "Node Type") +
#'     
#'     {if (!is.null(edge_probs))
#'       scale_edge_width(range = c(0.3, 2.5), name = "Edge\nProbability",
#'                        limits = c(0, 1))
#'     } +
#'     {if (!is.null(edge_probs))
#'       scale_edge_alpha(range = c(0.25, 1), guide = "none")
#'     } +
#'     
#'     # Theme
#'     theme_void(base_size = 12) +
#'     theme(
#'       legend.position = "right",
#'       legend.text = element_text(size = 10),
#'       legend.title = element_text(size = 11, face = "bold"),
#'       plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
#'                                 margin = margin(b = 10)),
#'       plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#7F8C8D",
#'                                    margin = margin(b = 10)),
#'       plot.margin = margin(15, 15, 15, 15)
#'     )
#'   
#'   if (!is.null(title)) {
#'     edge_count <- nrow(edges)
#'     p <- p + labs(
#'       title = title,
#'       subtitle = sprintf("%d nodes, %d edges", length(node_names), edge_count)
#'     )
#'   }
#'   
#'   return(p)
#' }
#' 
#' 
#' 
#' my_roles <- list(
#'   need_binary  = c(),
#'   need_ordinal = c("n_clust_basics", "n_clust_work_edu", "n_clust_discrim_iso", "n_HelpCaregive", "n_HelpChildcare", "n_HelpTransport",
#'                    "n_HelpInternet", "n_HelpLegal"),
#'   outcome      = c("o_mh", "o_wb")
#' )
#' 
#' my_short_names <- c(
#'   n_HelpCaregive      = "Caregiving",
#'   n_HelpChildcare     = "Childcare",
#'   n_HelpTransport     = "Transport",
#'   n_HelpInternet      = "Internet",
#'   n_HelpLegal         = "Legal",
#'   n_clust_basics      = "Basic\nNeeds",
#'   n_clust_work_edu    = "Work/\nEducation",
#'   n_clust_discrim_iso = "Discrim./\nIsolation",
#'   o_mh                = "Mental\nHealth",
#'   o_wb                = "Well-\nbeing"
#' )
#' # 
#' # 
#' # # --- Plot point estimate (no bootstrap) ---
#' # p1 <- plot_dag(
#' #   gam        = fit_point$gam,
#' #   node_names = dag_vars,
#' #   node_roles = my_roles,
#' #   short_names = my_short_names,
#' #   title      = "swa-oBN: Point Estimate",
#' #   layout     = "sugiyama"
#' # )
#' # print(p1)
#' # ggsave("dag_point_estimate.pdf", p1, width = 10, height = 7)
#' # 
#' # 
#' # # --- Plot bootstrap version (with edge probabilities) ---
#' # # Only show edges with inclusion prob > 0.50
#' # p2 <- plot_dag(
#' #   gam         = fit_boot$gam,
#' #   node_names  = dag_vars,
#' #   edge_probs  = fit_boot$edge_probs,
#' #   node_roles  = my_roles,
#' #   short_names = my_short_names,
#' #   threshold   = 0.50,
#' #   title       = "swa-oBN: Bootstrap (prob > 0.50)",
#' #   layout      = "sugiyama"
#' # )
#' # print(p2)
#' # ggsave("dag_bootstrap_050.pdf", p2, width = 10, height = 7)
#' # 
#' # 
#' # # --- Side-by-side comparison at different thresholds ---
#' # p3 <- plot_dag(
#' #   gam         = fit_boot$gam,
#' #   node_names  = dag_vars,
#' #   edge_probs  = fit_boot$edge_probs,
#' #   node_roles  = my_roles,
#' #   short_names = my_short_names,
#' #   threshold   = 0.30,
#' #   title       = "Bootstrap (prob > 0.30)",
#' #   layout      = "sugiyama"
#' # )
#' # 
#' # p4 <- plot_dag(
#' #   gam         = fit_boot$gam,
#' #   node_names  = dag_vars,
#' #   edge_probs  = fit_boot$edge_probs,
#' #   node_roles  = my_roles,
#' #   short_names = my_short_names,
#' #   threshold   = 0.70,
#' #   title       = "Bootstrap (prob > 0.70)",
#' #   layout      = "sugiyama"
#' # )
#' # 
#' # # Combine with patchwork if available
#' # if (requireNamespace("patchwork", quietly = TRUE)) {
#' #   library(patchwork)
#' #   combined <- p3 + p4 + plot_layout(guides = "collect")
#' #   print(combined)
#' #   ggsave("dag_threshold_comparison.pdf", combined, width = 16, height = 7)
#' # }
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 


















# ==============================================================================
# Publication-Quality DAG Plotting  (Suggestive / Robust edge support)
# Requires: ggraph, tidygraph, ggplot2, dplyr
# ==============================================================================

pacman::p_load(ggraph, tidygraph, ggplot2, dplyr, ggparty)


#' Plot a DAG from an adjacency matrix
#'
#' Edges are binned into support categories based on their bootstrap inclusion
#' probability:
#'   Weak        p in [0, prob_breaks[1])
#'   Suggestive  p in [prob_breaks[1], prob_breaks[2])
#'   Robust      p in [prob_breaks[2], 1]
#'
#' @param gam        Adjacency matrix (q x q). gam[i,j] = 1 means j -> i.
#' @param node_names Character vector of node names (length q)
#' @param edge_probs Optional matrix (q x q) of bootstrap edge inclusion probs.
#' @param node_roles Named list: list(need_ordinal = c(...), outcome = c(...))
#' @param threshold  Minimum edge probability to display (default 0.50)
#' @param title      Plot title
#' @param node_size  Base node size
#' @param label_size Node label text size
#' @param show_edge_labels Show probability on edges
#' @param layout     ggraph layout: "sugiyama", "fr", "kk", "stress", "circle"
#' @param short_names Optional named vector for shorter display labels
#' @param prob_breaks Numeric vector of length 2 giving the inner cut points
#'                   (default c(0.50, 0.75))
#' @param support_labels Character vector of length 3 naming the three bins,
#'                   in order of increasing support
#' @param show_counts_in_legend If TRUE, append "(n = k)" to each legend entry
#'
#' @return A ggplot2 object
plot_dag <- function(gam,
                     node_names = NULL,
                     edge_probs = NULL,
                     node_roles = NULL,
                     threshold = 0.50,
                     title = NULL,
                     subtitle = NULL,                 # NEW: markdown string
                     node_size = 30,
                     label_size = 5,
                     show_edge_labels = TRUE,
                     layout = "sugiyama",
                     short_names = NULL,
                     prob_breaks = c(0.50, 0.75),
                     support_labels = c("Weak (< 0.50)",
                                        "Suggestive (0.50\u20130.75)",
                                        "Robust (\u2265 0.75)"),
                     show_counts_in_legend = TRUE) {
  
  stopifnot(length(prob_breaks) == 2,
            length(support_labels) == 3,
            prob_breaks[1] < prob_breaks[2])
  
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
  
  # --- Assign support category and build legend labels ---
  # Only relevant when edge_probs is supplied; for the non-bootstrap case
  # we fall back to the original uniform styling further down.
  used_levels   <- character(0)
  legend_labels <- character(0)
  if (!is.null(edge_probs) && nrow(edges) > 0) {
    edges$support <- cut(
      edges$prob,
      breaks = c(-Inf, prob_breaks[1], prob_breaks[2], Inf),
      labels = support_labels,
      right = FALSE,
      include.lowest = TRUE
    )
    edges$support <- factor(edges$support,
                            levels = support_labels, ordered = TRUE)
    
    # Keep only levels that actually occur, in ascending order
    counts       <- table(edges$support)
    used_levels  <- support_labels[counts > 0]
    
    if (show_counts_in_legend) {
      legend_labels <- setNames(
        paste0(used_levels, "  (n = ", as.integer(counts[used_levels]), ")"),
        used_levels
      )
    } else {
      legend_labels <- setNames(used_levels, used_levels)
    }
  }
  
  # --- Remove isolated nodes ---
  connected_nodes <- unique(c(edges$from, edges$to))
  if (length(connected_nodes) < length(node_names)) {
    removed <- setdiff(node_names, connected_nodes)
    if (length(removed) > 0) {
      message(sprintf("Removed %d isolated node(s): %s",
                      length(removed), paste(removed, collapse = ", ")))
    }
    node_names <- node_names[node_names %in% connected_nodes]
  }
  
  # --- Build node table ---
  nodes <- data.frame(name = node_names, stringsAsFactors = FALSE)
  nodes$role <- "Other"
  if (!is.null(node_roles)) {
    if (!is.null(node_roles$need))    nodes$role[nodes$name %in% node_roles$need] <- "Need"
    if (!is.null(node_roles$outcome)) nodes$role[nodes$name %in% node_roles$outcome] <- "Outcome (ordinal)"
    if (!is.null(node_roles$need_ordinal))
      nodes$role[nodes$name %in% node_roles$need_ordinal] <- "Need (ordinal)"
    if (!is.null(node_roles$need_binary))
      nodes$role[nodes$name %in% node_roles$need_binary] <- "Need (binary)"
  }
  
  nodes$label <- nodes$name
  if (!is.null(short_names)) {
    matched <- match(nodes$name, names(short_names))
    has_short <- !is.na(matched)
    nodes$label[has_short] <- short_names[matched[has_short]]
  }
  
  # --- Build tidygraph ---
  g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  # --- Colour palettes ---
  role_colors <- c(
    "Outcome (ordinal)" = "#C0392B",
    "Need (ordinal)"    = "#2980B9",
    "Need (binary)"     = "#7F8C8D",
    "Other"             = "#BDC3C7"
  )
  role_fills <- c(
    "Outcome (ordinal)" = "#FADBD8",
    "Need (ordinal)"    = "#D6EAF8",
    "Need (binary)"     = "#EAECEE",
    "Other"             = "#F2F3F4"
  )
  
  # --- Support -> (linetype, width) mappings, indexed by the FULL label set
  #     so that changes to `support_labels` flow through automatically.
  support_linetype_values <- setNames(c("dotted", "dashed", "solid"),
                                      support_labels)
  support_width_values    <- setNames(c(0.4, 0.9, 2.0),
                                      support_labels)
  
  # --- Plot ---
  p <- ggraph(g, layout = layout) +
    
    # Edges
    {if (!is.null(edge_probs) && show_edge_labels)
      geom_edge_link(
        aes(edge_width = support, edge_linetype = support,
            label = sprintf("%.2f", prob)),
        arrow = arrow(length = unit(4, "mm"), type = "closed"),
        end_cap = circle(node_size * 0.35, "mm"),
        start_cap = circle(node_size * 0.35, "mm"),
        colour = "#2C3E50",
        alpha = 0.75,
        angle_calc = "along",
        label_dodge = unit(3, "mm"),
        label_size = 4.2,
        label_colour     = "#7B0000",      # dark red; swap "#1A237E" for navy
        label_fontface   = "bold",
        label_fill       = "white"         # white box behind each label
      )
      else if (!is.null(edge_probs))
        geom_edge_link(
          aes(edge_width = support, edge_linetype = support),
          arrow = arrow(length = unit(4, "mm"), type = "closed"),
          end_cap = circle(node_size * 0.35, "mm"),
          start_cap = circle(node_size * 0.35, "mm"),
          colour = "#2C3E50",
          alpha = 0.9,
          label_colour     = "#7B0000",      # dark red; swap "#1A237E" for navy
          label_fontface   = "bold",
          label_fill       = "white"         # white box behind each label
        )
      else
        geom_edge_link(
          arrow = arrow(length = unit(4, "mm"), type = "closed"),
          end_cap = circle(node_size * 0.35, "mm"),
          start_cap = circle(node_size * 0.35, "mm"),
          colour = "#2C3E50",
          edge_width = 0.25
        )
    } +
    
    # Node outlines
    geom_node_point(
      aes(colour = role),
      size = node_size, shape = 21, fill = "white", stroke = 2
    ) +
    
    # Node fills
    geom_node_point(
      aes(fill = role),
      size = node_size - 2, shape = 21, colour = NA
    ) +
    
    # Node labels
    geom_node_text(
      aes(label = label),
      size = label_size, colour = "#2C3E50", fontface = "bold"
    ) +
    
    # Node scales
    scale_colour_manual(values = role_colors, name = "Node type", guide = "none") +
    scale_fill_manual  (values = role_fills,  name = "Node type", guide = "none") +
    
    # Edge-support scales (merged by sharing name/breaks/labels)
    {if (!is.null(edge_probs))
      scale_edge_width_manual(
        values = support_width_values,
        breaks = used_levels,
        labels = legend_labels,
        name   = "Edge support",
        drop   = FALSE
      )
    } +
    {if (!is.null(edge_probs))
      scale_edge_linetype_manual(
        values = support_linetype_values,
        breaks = used_levels,
        labels = legend_labels,
        name   = "Edge support",
        drop   = FALSE
      )
    } +
    
    # Theme
    theme_void(base_size = 15) +
    theme(
      legend.position = "bottom",
      legend.text     = element_text(size = 15),
      legend.title    = element_text(size = 15, face = "bold"),
      legend.key.width = unit(1.2, "cm"),
      plot.title      = element_text(size = 15, face = "bold", hjust = 0.5,
                                     margin = margin(b = 10)),
      plot.subtitle   = element_text(size = 15, hjust = 0.5, color = "#2C3E50",
                                     margin = margin(b = 10)),
      plot.margin     = margin(10, 10, 10, 10)
    )
  
  # --- Title / subtitle ---
  if (!is.null(title) || !is.null(subtitle)) {
    labs_args <- list()
    if (!is.null(title))    labs_args$title    <- title
    if (!is.null(subtitle)) {
      labs_args$subtitle <- subtitle
    } else if (!is.null(title)) {
      # Fall back to the old edge-count subtitle when caller passes no subtitle
      labs_args$subtitle <- sprintf("%d nodes, %d edges",
                                    length(node_names), nrow(edges))
    }
    p <- p + do.call(labs, labs_args)
  }
  
  return(p)
}



# =============================================================================
# Patch to plot_dag: accept a pre-computed layout, keep isolated nodes
# =============================================================================
#
# Two changes to your existing plot_dag:
#   1. New argument `layout_df`: a data.frame with columns (name, x, y).
#      When supplied, it overrides `layout` and preserves all nodes.
#   2. When layout_df is supplied, the "remove isolated nodes" block is skipped
#      so that node positions stay identical across panels.
#
# Easiest integration: add these lines inside plot_dag.
#
# (A) After the function signature, before the "Build edge list" loop, add:
#
#     use_manual <- !is.null(layout_df)
#     if (use_manual) {
#       stopifnot(all(c("name", "x", "y") %in% names(layout_df)))
#       stopifnot(all(node_names %in% layout_df$name))
#     }
#
# (B) Wrap the "Remove isolated nodes" block so it's skipped for manual layouts:
#
#     if (!use_manual) {
#       connected_nodes <- unique(c(edges$from, edges$to))
#       if (length(connected_nodes) < length(node_names)) {
#         removed <- setdiff(node_names, connected_nodes)
#         if (length(removed) > 0) {
#           message(sprintf("Removed %d isolated node(s): %s",
#                           length(removed), paste(removed, collapse = ", ")))
#         }
#         node_names <- node_names[node_names %in% connected_nodes]
#       }
#     }
#
# (C) Change the ggraph() call to use the manual layout when supplied:
#
#     p <- if (use_manual) {
#       lay <- layout_df[match(nodes$name, layout_df$name), c("x", "y")]
#       ggraph(g, layout = "manual", x = lay$x, y = lay$y)
#     } else {
#       ggraph(g, layout = layout)
#     } +
#       ... (rest of plot unchanged)
#
# Also add `layout_df = NULL` to the function signature.


# =============================================================================
# plot_dag_aligned: variant of plot_dag that takes an explicit (x, y) layout
# and keeps isolated nodes. For side-by-side comparison panels.
# =============================================================================

pacman::p_load(ggraph, tidygraph, ggplot2, dplyr, patchwork)

plot_dag_aligned <- function(edge_probs,
                             node_names,
                             layout_df,                 # data.frame(name, x, y)
                             node_roles   = NULL,
                             short_names  = NULL,
                             threshold    = 0.50,
                             title        = NULL,
                             node_size    = 24,
                             label_size   = 4,
                             show_edge_labels = TRUE,
                             prob_breaks  = c(0.50, 0.75),
                             support_labels = c("Weak (< 0.50)",
                                                "Suggestive (0.50\u20130.75)",
                                                "Robust (\u2265 0.75)")) {
  
  stopifnot(nrow(edge_probs) == length(node_names),
            all(node_names %in% layout_df$name))
  
  q <- length(node_names)
  
  # --- Build edge list (filter at threshold) ---
  # Convention: edge_probs[i, j] = P(j -> i)
  ed_idx <- which(edge_probs >= threshold, arr.ind = TRUE)
  ed_idx <- ed_idx[ed_idx[, "row"] != ed_idx[, "col"], , drop = FALSE]
  edges <- data.frame(
    from = node_names[ed_idx[, "col"]],    # parent
    to   = node_names[ed_idx[, "row"]],    # child
    prob = edge_probs[ed_idx],
    stringsAsFactors = FALSE
  )
  
  # --- Support bins ---
  edges$support <- cut(
    edges$prob,
    breaks = c(-Inf, prob_breaks[1], prob_breaks[2], Inf),
    labels = support_labels, right = FALSE, include.lowest = TRUE
  )
  edges$support <- factor(edges$support, levels = support_labels, ordered = TRUE)
  
  # --- Nodes in a FIXED order matching layout_df ---
  # Critical: tbl_graph respects the order of the `nodes` data frame,
  # and ggraph's manual layout indexes x/y by that order.
  lay <- layout_df[match(node_names, layout_df$name), ]
  nodes <- data.frame(name = lay$name, x = lay$x, y = lay$y,
                      stringsAsFactors = FALSE)
  
  # Assign roles
  nodes$role <- "Other"
  if (!is.null(node_roles)) {
    if (!is.null(node_roles$need_ordinal))
      nodes$role[nodes$name %in% node_roles$need_ordinal] <- "Need (ordinal)"
    if (!is.null(node_roles$need_binary))
      nodes$role[nodes$name %in% node_roles$need_binary]  <- "Need (binary)"
    if (!is.null(node_roles$outcome))
      nodes$role[nodes$name %in% node_roles$outcome]      <- "Outcome (ordinal)"
  }
  
  # Display labels
  nodes$label <- nodes$name
  if (!is.null(short_names)) {
    m <- match(nodes$name, names(short_names))
    nodes$label[!is.na(m)] <- short_names[m[!is.na(m)]]
  }
  
  g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  # --- Color palette ---
  role_colors <- c("Outcome (ordinal)" = "#C0392B",
                   "Need (ordinal)"    = "#2980B9",
                   "Need (binary)"     = "#7F8C8D",
                   "Other"             = "#BDC3C7")
  role_fills  <- c("Outcome (ordinal)" = "#FADBD8",
                   "Need (ordinal)"    = "#D6EAF8",
                   "Need (binary)"     = "#EAECEE",
                   "Other"             = "#F2F3F4")
  
  support_linetype <- setNames(c("dotted", "dashed", "solid"), support_labels)
  support_width    <- setNames(c(0.4, 0.9, 2.0),              support_labels)
  
  # --- Build plot with MANUAL layout using nodes$x, nodes$y ---
  # Passing x and y explicitly here is the bit that was failing before.
  p <- ggraph(g, layout = "manual", x = nodes$x, y = nodes$y)
  
  if (nrow(edges) > 0) {
    if (show_edge_labels) {
      p <- p + geom_edge_arc(
        aes(edge_width = support, edge_linetype = support,
            label = sprintf("%.2f", prob)),
        strength    = 0.12,                       # <- NEW: curvature; 0 = straight
        arrow       = arrow(length = unit(3.5, "mm"), type = "closed"),
        end_cap     = circle(node_size * 0.40, "mm"),
        start_cap   = circle(node_size * 0.40, "mm"),
        colour      = "#2C3E50", alpha = 0.9,
        angle_calc  = "along",
        label_pos   = 0.70,
        label_dodge = unit(3.5, "mm"),
        label_size  = 4,
        label_colour     = "#7B0000",      # dark red; swap "#1A237E" for navy
      )
    } else {
      p <- p + geom_edge_arc(
        aes(edge_width = support, edge_linetype = support),
        strength    = 0.12,
        arrow       = arrow(length = unit(3.5, "mm"), type = "closed"),
        end_cap     = circle(node_size * 0.40, "mm"),
        start_cap   = circle(node_size * 0.40, "mm"),
        colour      = "#2C3E50", alpha = 0.9,
        label_colour     = "#7B0000"  
      )
    }
  }
  
  
  
  p <- p +
    geom_node_point(aes(colour = role),
                    size = node_size, shape = 21, fill = "white", stroke = 2) +
    geom_node_point(aes(fill = role),
                    size = node_size - 2, shape = 21, colour = NA) +
    geom_node_text(aes(label = label),
                   size = label_size, colour = "#2C3E50", fontface = "bold") +
    scale_colour_manual(values = role_colors, guide = "none") +
    scale_fill_manual  (values = role_fills,  guide = "none") +
    scale_edge_width_manual   (values = support_width,    name = "Edge support",
                               drop = FALSE) +
    scale_edge_linetype_manual(values = support_linetype, name = "Edge support",
                               drop = FALSE) +
    theme_void(base_size = 13) +
    theme(legend.position  = "bottom",
          legend.text      = element_text(size = 12),
          legend.title     = element_text(size = 12, face = "bold"),
          legend.key.width = unit(1.0, "cm"),
          plot.title       = element_text(size = 13, face = "bold",
                                          hjust = 0.5, margin = margin(b = 8)),
          plot.margin      = margin(8, 8, 8, 8))
  
  if (!is.null(title)) p <- p + labs(title = title)
  p
}
