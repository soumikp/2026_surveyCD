library(igraph)

plot_cpdag <- function(result, layout_fn = layout_with_sugiyama,
                       vertex_size = 30, vertex_color = "white",
                       edge_color_directed = "black",
                       edge_color_undirected = "grey40",
                       label_cex = 1.0, edge_width = 1.5,
                       main = "CPDAG", ...) {
  # result: output from pc_stable_weighted (needs $amat, $varnames)
  
  amat <- result$amat
  p <- nrow(amat)
  vn <- result$varnames
  
  # Build edge list with type info
  edges_directed <- list()
  edges_undirected <- list()
  
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i == j) next
      if (amat[i, j] == 1 && amat[j, i] == 0) {
        # Directed: i -> j
        edges_directed[[length(edges_directed) + 1]] <- c(i, j)
      } else if (amat[i, j] == 1 && amat[j, i] == 1 && i < j) {
        # Undirected: i -- j (only once)
        edges_undirected[[length(edges_undirected) + 1]] <- c(i, j)
      }
    }
  }
  
  # Create a directed igraph object
  el <- do.call(rbind, c(edges_directed, edges_undirected))
  
  if (is.null(el) || nrow(el) == 0) {
    # No edges — plot isolated nodes
    g <- make_empty_graph(n = p, directed = TRUE)
    V(g)$name <- vn
    plot(g, vertex.size = vertex_size, vertex.color = vertex_color,
         vertex.label.cex = label_cex, main = main, ...)
    return(invisible(g))
  }
  
  n_dir <- length(edges_directed)
  n_undir <- length(edges_undirected)
  
  # For undirected edges, add both directions so layout algorithms see them,
  # but we'll draw them with no arrow
  if (n_undir > 0) {
    reverse_undir <- do.call(rbind, lapply(edges_undirected, rev))
    el_full <- rbind(el, reverse_undir)
  } else {
    el_full <- el
  }
  
  # Build graph from named edge list so vertex names are correct
  el_named <- matrix(vn[el_full], ncol = 2)
  g <- graph_from_edgelist(el_named, directed = TRUE)
  
  # Ensure all vertices present (isolated nodes)
  missing <- setdiff(vn, V(g)$name)
  if (length(missing) > 0) g <- g + vertices(missing)
  
  # Edge attributes
  n_total <- ecount(g)
  arrow_mode <- rep(2, n_total)       # 2 = forward arrow (directed)
  e_color <- rep(edge_color_directed, n_total)
  e_lty <- rep(1, n_total)
  
  # Mark undirected edges (both directions) as no-arrow, dashed
  el_g <- as_edgelist(g)
  for (k in seq_len(n_total)) {
    from <- el_g[k, 1]
    to <- el_g[k, 2]
    is_undir <- any(sapply(edges_undirected, function(e) {
      (vn[e[1]] == from && vn[e[2]] == to) || (vn[e[1]] == to && vn[e[2]] == from)
    }))
    if (is_undir) {
      arrow_mode[k] <- 0  # no arrow
      e_color[k] <- edge_color_undirected
    }
  }
  
  # Remove duplicate undirected edges (keep only i<j direction for drawing)
  keep <- rep(TRUE, n_total)
  for (k in seq_len(n_total)) {
    if (arrow_mode[k] == 0) {
      from <- el_g[k, 1]; to <- el_g[k, 2]
      if (from > to) keep[k] <- FALSE
    }
  }
  g_plot <- delete_edges(g, which(!keep))
  arrow_mode <- arrow_mode[keep]
  e_color <- e_color[keep]
  e_lty <- e_lty[keep]
  
  # Layout
  if (identical(layout_fn, layout_with_sugiyama)) {
    lay <- layout_fn(g_plot)$layout
  } else {
    lay <- layout_fn(g_plot)
  }
  
  plot(g_plot,
       layout = lay,
       vertex.size = vertex_size,
       vertex.color = vertex_color,
       vertex.frame.color = "black",
       vertex.label = vn,
       vertex.label.cex = label_cex,
       vertex.label.color = "black",
       edge.arrow.size = 0.5,
       edge.arrow.mode = arrow_mode,
       edge.color = e_color,
       edge.width = edge_width,
       edge.lty = e_lty,
       main = main,
       ...)
  
  invisible(g_plot)
}

# ==============================================================================
# Usage
# ==============================================================================
# plot_cpdag(result)                          # Sugiyama (layered DAG) layout
# plot_cpdag(result, layout_fn = layout_nicely)
# plot_cpdag(result, layout_fn = layout_in_circle)
# plot_cpdag(result, vertex_color = "lightblue", main = "My CPDAG")