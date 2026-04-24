
rm(list=ls())
source("2025_11_17_cleanData.R") 

pacman::p_load(here, VGAM, lavaan, bnlearn,
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


pacman::p_load(bnlearn, psych) # For polychoric correlation

# OUTCOME Cluster (To be kept separate from each other)
outcomes <- c("wellbeing", "mentalHealth")

# PREDICTOR Clusters (Internal edges allowed, Cross-cluster edges forbidden)
c_material <- c("housing", "basics", "food")
c_ses      <- c("work", "education")
c_psycho   <- c("discrimination_rec", "discrimination_gender", "discrimination_age", "discrimination_other", "isolation", "loneliness")

# Independent variables (Treating them as clusters of size 1)
c_transport <- c("transport")
c_internet  <- c("internet")
c_childcare <- c("childcare")
c_caregiving<- c("caregiving")
c_legal     <- c("legal")

predictor_clusters_list <- list(
  c_material,
  c_ses,
  c_psycho,
  c_transport,
  c_internet,
  c_childcare,
  c_caregiving,
  c_legal
)

all_predictors <- unlist(predictor_clusters_list)


# No arcs between any need clusters.
bl_silos <- data.frame(from = character(), to = character())

n_clusters <- length(predictor_clusters_list)

for (i in 1:(n_clusters - 1)) {
  for (j in (i + 1):n_clusters) {
    # Get nodes in cluster i and cluster j
    nodes_i <- predictor_clusters_list[[i]]
    nodes_j <- predictor_clusters_list[[j]]
    
    # Block i -> j AND j -> i
    # expand.grid creates all combinations of nodes between the two sets
    bl_silos <- rbind(bl_silos, expand.grid(from = nodes_i, to = nodes_j))
    bl_silos <- rbind(bl_silos, expand.grid(from = nodes_j, to = nodes_i))
  }
}

# This allows Predictor -> Outcome, but forbids Outcome -> Predictor
tiers <- list(all_predictors, outcomes)
bl_reverse <- tiers2blacklist(tiers)

# I dont want any arrow between the two outcomes
bl_outcomes <- data.frame(
  from = c("wellbeing", "mentalHealth"),
  to   = c("mentalHealth", "wellbeing")
)

bl_final <- rbind(bl_silos, 
                  bl_reverse, bl_outcomes)
bl_final <- unique(bl_final) # Remove duplicates

# ==============================================================================
# 3. VERIFICATION & EXECUTION
# ==============================================================================

# Print summary of constraints
cat("Total constraints enforcing modular structure:", nrow(bl_final), "\n")

# Run the Structure Learning (Example using Hill-Climbing)
# Assuming 'shep_data' is your dataset
outcomes_risks <- data.frame(outcomes_risks)
colnames(outcomes_risks)[colnames(outcomes_risks) == "mental.health"] <- "mentalHealth"
outcomes_risks <- droplevels(outcomes_risks)
bn_structure <- hc(outcomes_risks, blacklist = bl_final)

# VISUALIZATION
plot(bn_structure)
