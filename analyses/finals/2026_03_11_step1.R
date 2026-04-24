rm(list = ls())
pacman::p_load(here, tidyverse, pheatmap, reshape2, stringr)
source("C:/Users/VHAPTHPURKAS/OneDrive - University of Pittsburgh/CHERP/2025_gordon_ccdf/analysis/2025_05_28_02_multinomial.R")
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
dplyr::select(
  -r_DiscrimREC, -r_DiscrimGend, -r_DiscrimAge, -r_DiscrimOther,
  -r_FeelLonely, -r_LackSocial, -r_PayBasics, -r_LackFood, -r_HaveHousing,
  -r_FindWork, -r_LackEduc,
  
  -n_HelpBasics, -n_HelpFood, -n_HelpHousing,
  -n_HelpWork, -n_HelpEduc,
  -n_HelpDiscrim, -n_HelpIsolation, -n_HelpLoneliness
)

rm(list = ls()[!ls() %in% c("dat_fin", "shep_clustered")])
data <- shep_clustered |> mutate(SURVID = as.numeric(SURVID)) |> left_join(dat_fin |> dplyr::select(SURVID, age_v4, race_sex_v2))


