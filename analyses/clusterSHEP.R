pacman::p_load(here, tidyverse, pheatmap, reshape2, stringr)
source(file.path(here(), "analyses", "correlationsSHEP.R"))

shep_clustered <- data %>%
  rowwise() %>%
  mutate(
    # ==========================================
    r_clust_discrimination = round(mean(c(as.numeric(r_DiscrimREC), as.numeric(r_DiscrimGend), 
                                          as.numeric(r_DiscrimAge), as.numeric(r_DiscrimOther)), na.rm = TRUE)),
    
    r_clust_isolation = round(mean(c(as.numeric(r_FeelLonely), as.numeric(r_LackSocial)), na.rm = TRUE)),
    
    r_clust_basics = round(mean(c(as.numeric(r_PayBasics), as.numeric(r_LackFood), 
                                  as.numeric(r_HaveHousing)), na.rm = TRUE)),
    
    r_clust_work_edu = round(mean(c(as.numeric(r_FindWork), as.numeric(r_LackEduc)), na.rm = TRUE)),
    
    n_clust_basics = round(mean(c(as.numeric(n_HelpBasics), as.numeric(n_HelpFood), 
                                  as.numeric(n_HelpHousing)), na.rm = TRUE)),
    
    n_clust_work_edu = round(mean(c(as.numeric(n_HelpWork), as.numeric(n_HelpEduc)), na.rm = TRUE)),
    
    n_clust_discrim_iso = round(mean(c(as.numeric(n_HelpDiscrim), as.numeric(n_HelpIsolation), 
                                       as.numeric(n_HelpLoneliness)), na.rm = TRUE))
  ) %>%
  ungroup() %>%
mutate(across(
  c(r_clust_discrimination, r_clust_isolation, r_clust_basics, r_clust_work_edu,
    n_clust_basics, n_clust_work_edu, n_clust_discrim_iso),
  ~ factor(.x, ordered = TRUE)
)) %>%
# ==========================================
select(
  -r_DiscrimREC, -r_DiscrimGend, -r_DiscrimAge, -r_DiscrimOther,
  -r_FeelLonely, -r_LackSocial, -r_PayBasics, -r_LackFood, -r_HaveHousing,
  -r_FindWork, -r_LackEduc,
  
  -n_HelpBasics, -n_HelpFood, -n_HelpHousing,
  -n_HelpWork, -n_HelpEduc,
  -n_HelpDiscrim, -n_HelpIsolation, -n_HelpLoneliness
)




clust_cor <- lavaan::lavCor(data.frame(shep_clustered[-c(1)]), 
                           ordered = colnames(shep_clustered)[-c(1, 2)], 
                           sampling.weights = "WEIGHT")

order <- sort(colnames(shep_clustered)[-c(1, 2)])

melted_cor <- melt(clust_cor)

melted_cor$Var1 <- factor(melted_cor$Var1, levels = order)
melted_cor$Var2 <- factor(melted_cor$Var2, levels = rev(order))

corplot <- melted_cor |> 
  ggplot(aes(x=Var1, y = Var2, fill = value)) + 
  geom_tile(color = "black") +
  geom_text(aes(label = sprintf("%0.2f", value))) + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16), 
        axis.text.y = element_text(hjust = 1, size = 16), 
        plot.title = element_text(hjust = 0, size = 20, face = "bold"), 
        plot.subtitle = element_text(hjust = 0, size = 20, face = "bold"), 
        legend.position = "none") + 
  scale_fill_gsea(limits = c(-1, 1)) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  labs(x = "", y = "", 
       title = str_wrap("Survey-weighted polychoric correlations between clustered needs (n), clustered risks (r), and outcomes (o).", 90), 
       subtitle = "n = 5396 complete cases")

factor <- 1.5
ggsave(file.path(here(), "analyses", "2026_02_24_cluster_correlations.pdf"), 
       corplot, 
       height = factor*10, 
       width = factor*10, 
       units = "in", 
       device = "pdf")
