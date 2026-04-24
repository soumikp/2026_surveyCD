rm(list = ls())
pacman::p_load(here, tidyverse, pheatmap, reshape2, stringr, dplyr, bnlearn)
source(file.path(here(), "analyses", "clusterSHEP.R"))

# bnlearn expects a 2-column data frame: "from" and "to"
create_bnlearn_blacklist <- function(var_names) {
  
  get_tier <- function(name) {
    if (startsWith(name, "r_")) return(1) # Tier 1: Risks
    if (startsWith(name, "n_")) return(2) # Tier 2: Needs
    if (startsWith(name, "o_")) return(3) # Tier 3: Outcomes
    return(0)
  }
  
  # Initialize empty data frame
  blacklist <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
  
  # Loop through all pairs
  for (i in var_names) {
    for (j in var_names) {
      tier_from <- get_tier(i)
      tier_to <- get_tier(j)
      
      # If 'from' is a HIGHER tier than 'to', it's a backward edge. Forbid it!
      if (tier_from > tier_to) {
        blacklist <- rbind(blacklist, data.frame(from = i, to = j))
      }
    }
  }
  
  return(blacklist)
}

# Install these if you don't have them yet:
# install.packages(c("igraph", "ggraph", "ggplot2", "patchwork"))

pacman::p_load(igraph, ggraph, ggplot2, patchwork, dplyr)

# ==============================================================================
# 1. THE BEAUTIFUL PLOT GENERATOR
# ==============================================================================
plot_beautiful_dag <- function(bn_dag, title) {
  
  # 1. Convert bnlearn object to an igraph object
  g <- bnlearn::as.igraph(bn_dag)
  
  # 2. Extract node names to assign tiers for coloring
  nodes <- V(g)$name
  
  # Assign categories based on your prefixes
  tiers <- case_when(
    startsWith(nodes, "r_") ~ "1. Social Risks",
    startsWith(nodes, "n_") ~ "2. Social Needs",
    startsWith(nodes, "o_") ~ "3. Health Outcomes",
    TRUE ~ "Other"
  )
  V(g)$tier <- tiers
  
  # 3. Build the ggraph plot
  p <- ggraph(g, layout = "sugiyama") + 
    # Draw the arrows with a gap at the end so they don't pierce the labels
    geom_edge_link(arrow = arrow(length = unit(3, 'mm'), type = "closed"),
                   end_cap = circle(8, 'mm'), 
                   start_cap = circle(8, 'mm'),
                   color = "gray40", 
                   width = 1.2, 
                   alpha = 0.8) +
    # Draw the nodes as clean, colored text boxes
    geom_node_label(aes(label = name, fill = tier), 
                    color = "black", 
                    fontface = "bold", 
                    size = 3.5,
                    label.padding = unit(0.4, "lines"),
                    label.r = unit(0.3, "lines")) +
    # Define your exact color palette
    scale_fill_manual(values = c("1. Social Risks" = "#A9D0F5",      # Soft Blue
                                 "2. Social Needs" = "#A9F5A9",      # Soft Green
                                 "3. Health Outcomes" = "#F5A9A9")) + # Soft Red
    # Remove all background grids and axes
    theme_void() +
    # Style the title and legend
    labs(title = title, fill = "Domain") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 15)),
          legend.position = "bottom",
          legend.title = element_text(face = "bold"))
  
  return(p)
}


# Ensure your data is stored as factors (bnlearn requires this for discrete data)
df_clean <- shep_clustered %>% mutate(across(everything(), as.factor))

cat("\nRunning DAG 1 (Risks -> Outcomes)...\n")
df_dag1 <- df_clean %>% select(starts_with("r_"), starts_with("o_"))
bl_dag1 <- create_bnlearn_blacklist(colnames(df_dag1))
res_dag1 <- hc(data.frame(df_dag1), blacklist = bl_dag1, score = "bic") 
ggsave(file.path(here(), "analyses", "2026_02_25_clusteredRisksToOutcomes.pdf"), 
       plot_beautiful_dag(res_dag1, "DAG 1 (Risks -> Outcomes)"),
       height = 10, 
       width = 15, 
       units = 'in', 
       device = "pdf")
# ==============================================================================
# 3. DAG 2: NEEDS AND OUTCOMES ONLY
# ==============================================================================
cat("Running DAG 2 (Needs -> Outcomes)...\n")
df_dag2 <- df_clean %>% select(starts_with("n_"), starts_with("o_"))
bl_dag2 <- create_bnlearn_blacklist(colnames(df_dag2))
res_dag2 <-  hc(data.frame(df_dag2), blacklist = bl_dag2, score = "bic") 
ggsave(file.path(here(), "analyses", "2026_02_25_clusteredNeedsToOutcomes.pdf"), 
       plot_beautiful_dag(res_dag2, "DAG 2 (Needs -> Outcomes)"),
       height = 10, 
       width = 15, 
       units = 'in', 
       device = "pdf")


# ==============================================================================
# 4. DAG 3: RISKS, NEEDS, AND OUTCOMES
# ==============================================================================
cat("Running DAG 3 (Risks -> Needs -> Outcomes)...\n")
df_dag3 <- df_clean %>% select(starts_with("r_"), starts_with("n_"), starts_with("o_"))
bl_dag3 <- create_bnlearn_blacklist(colnames(df_dag3))
res_dag3 <-  hc(data.frame(df_dag3), blacklist = bl_dag3, score = "bic") 
ggsave(file.path(here(), "analyses", "2026_02_25_clusteredRisksToClusteredNeedsToOutcomes.pdf"), 
       plot_beautiful_dag(res_dag3, "DAG 3 (Risks -> Needs -> Outcomes)"),
       height = 10, 
       width = 15, 
       units = 'in', 
       device = "pdf")

