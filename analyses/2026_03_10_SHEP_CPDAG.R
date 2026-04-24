rm(list = ls())
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

source(file.path(here(), "analyses", "2026_03_10_weightedOrdinalPC.R"))

cpdag <- pc_stable_weighted(data.frame(shep_clustered[,str_detect(colnames(shep_clustered), "o_|n_")]), 
                            shep_clustered |> pull("WEIGHT"))
                        

source(file.path(here(), "analyses", "2026_03_10_cpdagPlot.R"))
plot_cpdag(cpdag, layout_fn = layout_in_circle)

edges <- which(cpdag$amat == 1, arr.ind = TRUE)

cbind(cpdag$varnames[edges[,1]], cpdag$varnames[edges[,2]])
