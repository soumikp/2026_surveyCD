rm(list=ls())
pacman::p_load(here, VGAM, lavaan, 
               pheatmap, qgraph, patchwork, ggpubr)
source(file.path(here(), "analysis", "2025_06_04_numericWBS_summary.R"))

outcome_need<-c('HelpBasics_T1'
                ,'HelpCaregive_T1'
                ,'HelpChildcare_T1'
                ,'HelpWork_T1'
                ,'HelpFood_T1'
                ,'HelpHousing_T1'
                ,'HelpTransport_T1'
                ,'HelpInternet_T1'
                ,'HelpIsolation_T1'
                ,'HelpLoneliness_T1'
                ,'HelpDiscrim_T1'
                ,'HelpLegal_T1'
                ,'HelpEduc_T1'
)

need.binary<-dat_fin%>%
  select(all_of(outcome_need)
         ,SURVID)

#create new grid variables to indicate support received or not received 

#create binary indicator for support needed

need.binary <- need.binary %>%  
  mutate(across(                                      
    .cols = all_of(c(outcome_need)),  ## for each column listed and "outcome"
    .fns = ~case_when(                              
      . %in% c(2, 3)   ~ 2,           ## needs support =2
      . %in% c(1) ~ 1,           ##  no support needed=1
      TRUE                            ~ NA_real_)    ## otherwise set to missing
  ))%>%
  mutate_at(outcome_need,factor)%>%
  mutate(across((!SURVID), ~fct_recode(.,"No support needed"= "1",
                                       "Needs support" = "2")))%>%
  mutate(across((!SURVID),~fct_relevel(., "No support needed", "Needs support")))

outcomes <- dat_fin |> 
  rowwise() |> 
  mutate(wellbeing_outcome = sum(10*(as.numeric(WBS_FullySat_T1) - 1),
                                 10*(as.numeric(WBS_Involved_T1) - 1), 
                                 10*(as.numeric(WBS_Function_T1) - 1), na.rm = TRUE )/3) |> 
  ungroup() |> 
  mutate(mentalHealth_outcome = RateMentalHealth_T1) |> 
  rowwise() |> 
  mutate(mentalHealth_outcome_v2 = case_when(mentalHealth_outcome == "Poor" ~ "Negative", 
                                             (mentalHealth_outcome == "Excellent" | 
                                                mentalHealth_outcome == "Very good" | 
                                                mentalHealth_outcome == "Good" | 
                                                mentalHealth_outcome == "Fair" ) ~ "Positive", 
                                             TRUE ~ NA_character_)) |> 
  ungroup() |> 
  select(wellbeing_outcome, mentalHealth_outcome_v2, SURVID)


#set reference level
dat_fin$age_v4<- relevel(factor(dat_fin$age_v4, 
                                levels = c(1,2,3,4,5), 
                                labels = c("18 to 44", "45 to 54", "55 to 64", "65 to 74", "75 or more"))
                         , ref = "65 to 74")
dat_fin$sexor_v1<- relevel(as.factor(dat_fin$sexor_v1), ref = "Straight")
dat_fin$race_sex_v2 <- relevel(as.factor(dat_fin$race_sex_v2), ref = "White men") 

dat_ccdf <- dat_fin |> 
  select(SURVID, age_v4, sexor_v1, race_sex_v2, WEIGHT) |> 
  left_join(outcomes) |> 
  left_join(need.binary)

dat_ccdf$mentalHealth_outcome_v2 <- relevel(as.factor(dat_ccdf$mentalHealth_outcome_v2), ref = "Positive")
rm(list = ls()[!ls() %in% c("dat_ccdf", "outcome_need")])

dat_ccdf$wellbeing_outcome <- cut(dat_ccdf$wellbeing_outcome, 
                                  breaks = c(0, 25, 50, 75, 100), 
                                  labels = c("Q1 (worst)", "Q2", "Q3", "Q4 (best)"))


nodes <- colnames(dat_ccdf)[-c(1:5)]
new_names <- c("wellbeing", "mentalHealth", "basics", 
               "caregiving", "childcare", "work", "food", "housing", 
               "transport", "internet", "isolation", "loneliness", "discrimination", 
               "legal", "education")
data <- dat_ccdf |> dplyr::select(all_of(nodes)) |> drop_na()
colnames(data) <- new_names
data$mentalHealth <- factor(data$mentalHealth, levels = c("Negative", "Positive"))

cor_matrix <- lavaan::lavCor(data, ordered = colnames(data))

cor_map <- pheatmap(
  cor_matrix,
  main = "Polychoric Correlation Heatmap",
  fontsize = 12,
  display_numbers = TRUE,  # Show the correlation coefficient in each cell
  cluster_rows = TRUE,   # Cluster rows to group similar variables
  cluster_cols = TRUE#,   # Cluster columns
  #color = colorRampPalette(c("#67001F", "white", "#053061"))(50) # Red-White-Blue palette
)

factor <- 1.25
ggsave("C:/Users/VHAPTHPURKAS/OneDrive - University of Pittsburgh/CHERP/2024_shep_hausmann/2025_11_14_ordCD/2025_11_14_corPlot.pdf", 
       cor_map, 
       height = factor*8.5,
       width = factor*10, 
       units = "in", 
       device = "pdf")

cor_mat <- matrix(as.numeric(cor_matrix), ncol = 15, nrow = 15)
colnames(cor_mat) <- rownames(cor_mat) <- colnames(data)


undir_graph_all <- qgraph(
  cor_mat,
  graph = "cor",            # Use correlation-specific settings
  layout = "spring",        # Use a force-directed layout
  labels = colnames(data), # Use variable names as labels
  theme = "colorblind",
  #minimum = quantile(cor_mat, 2/3),            # Hide weak correlations (edges) < |0.1|
  details = TRUE,           # Show the legend
  title = "Correlation Network (no trimming)", 
  label.scale = FALSE
)



undir_graph_trim <- qgraph(
  cor_mat,
  graph = "cor",            # Use correlation-specific settings
  layout = "spring",        # Use a force-directed layout
  labels = colnames(data), # Use variable names as labels
  theme = "colorblind",
  minimum = quantile(cor_mat, 0.5),            # Hide weak correlations (edges) < |0.1|
  details = TRUE,           # Show the legend
  title = "Correlation Network (trimming 1/2 of edges)", 
  label.scale = FALSE
)

