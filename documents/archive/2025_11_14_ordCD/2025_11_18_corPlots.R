rm(list=ls())
source("2025_11_17_cleanData.R") 
#
pacman::p_load(here, VGAM, lavaan, 
               pheatmap, qgraph, patchwork, ggpubr, stringr)

dat_fin <- dat_fin |> 
  rowwise() |> 
  mutate(wellbeing_outcome = sum(10*(as.numeric(WBS_FullySat_T1) - 1),
                                 10*(as.numeric(WBS_Involved_T1) - 1), 
                                 10*(as.numeric(WBS_Function_T1) - 1), na.rm = TRUE )/3) |> 
  ungroup() |> 
  mutate(mentalHealth_outcome = RateMentalHealth_T1)

dat_fin <- dat_fin |> select(-contains("WBS")) ## dropping uncoded wellbeing outcomes
dat_fin <- dat_fin |> select(-contains("Rate")) ## dropping uncoded mental health outcomes

dat_fin$wellbeing_outcome <- cut(dat_fin$wellbeing_outcome, 
                                  breaks = c(0, 25, 50, 75, 100), 
                                  labels = c("Q1 (worst)", "Q2", "Q3", "Q4 (best)"))


needs <- colnames(dat_fin)[grepl("help", colnames(dat_fin), ignore.case = TRUE)]
needs_names <- c("basics", "caregiving", "childcare", "work", "food", "housing", 
                 "transport", "internet", "isolation", "loneliness", "discrimination", 
                 "legal", "education")

risks <- colnames(dat_fin)[grepl("_T1", colnames(dat_fin), ignore.case = TRUE)]
risks <- risks[!(risks %in% needs)]
risks_names <- c("basics", "caregiving", "childcare", "work", "food", "housing", 
                 "transport", "internet", "isolation", "loneliness", 
                 "discrimination_rec", "discrimination_gender", "discrimination_age",
                 "discrimination_other", "legal", "education")
outcomes <- colnames(dat_fin)[grepl("outcome", colnames(dat_fin), ignore.case = TRUE)]
outcomes_names <- c("wellbeing", "mental health")



#### outcomes with risks ####
outcomes_risks <- dat_fin |> 
  select(all_of(c(risks, outcomes))) |> 
  drop_na()

outcomes_risks$HaveHousing_T1 <- 5 - as.numeric(outcomes_risks$HaveHousing_T1)
outcomes_risks$AccessInternet_T1 <- 5 - as.numeric(outcomes_risks$AccessInternet_T1)
outcomes_risks$LegalNeed_T1 <- factor(3 - as.numeric(outcomes_risks$LegalNeed_T1), levels = c("1", "2"))

outcomes_risks <- outcomes_risks |> 
  mutate(across(all_of(risks[-8]), 
                ~case_match(as.character(.x), c("1", "8") ~ "1", .default = as.character(.x)))) |> 
  mutate(across(all_of(risks),~factor(.x, levels = c("1", "2", "3", "4"))))


outcomes_risks$mentalHealth_outcome <- factor(outcomes_risks$mentalHealth_outcome, 
                                              levels = c(0, 1, 2, 3, 4))
colnames(outcomes_risks) <- c(risks_names, outcomes_names)

cor_matrix <- lavaan::lavCor(outcomes_risks, ordered = colnames(outcomes_risks))

cor_map_risks <- pheatmap(
  cor_matrix,
  main = ("Polychoric Correlation Heatmap: Outcomes vs Risks\n(UPDATED 12.04.25 - ordering fixed for housing, legal, internet)"), 
  fontsize = 12,
  display_numbers = TRUE,  # Show the correlation coefficient in each cell
  cluster_rows = TRUE,   # Cluster rows to group similar variables
  cluster_cols = TRUE,#,   # Cluster columns
  #color = colorRampPalette(c("#67001F", "white", "#053061"))(50) # Red-White-Blue palette
)

factor <- 1.5
ggsave("C:/Users/VHAPTHPURKAS/OneDrive - University of Pittsburgh/CHERP/2024_shep_hausmann/2025_11_14_ordCD/2025_11_18_corPlot_risks.pdf", 
       cor_map_risks, 
       height = factor*8.5,
       width = factor*10, 
       units = "in", 
       device = "pdf")


cor_mat <- matrix(as.numeric(cor_matrix), ncol = ncol(cor_matrix), nrow =  ncol(cor_matrix))
colnames(cor_mat) <- rownames(cor_mat) <- colnames(outcomes_risks)

undir_graph_all_risks <- qgraph(
  cor_mat,
  graph = "cor",            # Use correlation-specific settings
  layout = "spring",        # Use a force-directed layout
  labels = colnames(outcomes_risks), # Use variable names as labels
  theme = "colorblind",
  #minimum = quantile(cor_mat, 2/3),            # Hide weak correlations (edges) < |0.1|
  details = TRUE,           # Show the legend
  title = "Correlation Network (no trimming)", 
  label.scale = FALSE
)

undir_graph_trim_risks <- qgraph(
  cor_mat,
  graph = "cor",            # Use correlation-specific settings
  layout = "spring",        # Use a force-directed layout
  labels = colnames(outcomes_risks), # Use variable names as labels
  theme = "colorblind",
  minimum = quantile(cor_mat, 0.5),            # Hide weak correlations (edges) < |0.1|
  details = TRUE,           # Show the legend
  title = "Correlation Network (trimming 1/2 of edges)", 
  label.scale = FALSE
)





#### outcomes with needs ####
outcomes_needs <- dat_fin |> 
  select(all_of(c(needs, outcomes))) |> 
  drop_na()

outcomes_needs$mentalHealth_outcome <- factor(outcomes_needs$mentalHealth_outcome, 
                                              levels = c(0, 1, 2, 3, 4))
colnames(outcomes_needs) <- c(needs_names, outcomes_names)

cor_matrix <- lavaan::lavCor(outcomes_needs, ordered = colnames(outcomes_needs))

cor_map_needs <- pheatmap(
  cor_matrix,
  main = "Polychoric Correlation Heatmap",
  fontsize = 12,
  display_numbers = TRUE,  # Show the correlation coefficient in each cell
  cluster_rows = TRUE,   # Cluster rows to group similar variables
  cluster_cols = TRUE#,   # Cluster columns
  #color = colorRampPalette(c("#67001F", "white", "#053061"))(50) # Red-White-Blue palette
)

factor <- 1.5
ggsave("C:/Users/VHAPTHPURKAS/OneDrive - University of Pittsburgh/CHERP/2024_shep_hausmann/2025_11_14_ordCD/2025_11_18_corPlot_needs.pdf", 
       cor_map_needs, 
       height = factor*8.5,
       width = factor*10, 
       units = "in", 
       device = "pdf")


cor_mat <- matrix(as.numeric(cor_matrix), ncol = ncol(cor_matrix), nrow =  ncol(cor_matrix))
colnames(cor_mat) <- rownames(cor_mat) <- colnames(outcomes_needs)

undir_graph_all_needs <- qgraph(
  cor_mat,
  graph = "cor",            # Use correlation-specific settings
  layout = "spring",        # Use a force-directed layout
  labels = colnames(outcomes_needs), # Use variable names as labels
  theme = "colorblind",
  #minimum = quantile(cor_mat, 2/3),            # Hide weak correlations (edges) < |0.1|
  details = TRUE,           # Show the legend
  title = "Correlation Network of outcomes with needs (no trimming)", 
  label.scale = FALSE
)

undir_graph_trim_needs <- qgraph(
  cor_mat,
  graph = "cor",            # Use correlation-specific settings
  layout = "spring",        # Use a force-directed layout
  labels = colnames(outcomes_needs), # Use variable names as labels
  theme = "colorblind",
  minimum = quantile(cor_mat, 0.5),            # Hide weak correlations (edges) < |0.1|
  details = TRUE,           # Show the legend
  title = "Correlation Network of outcomes with needs (trimming 1/2 of edges)", 
  label.scale = FALSE
)

