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

# dat_fin$wellbeing_outcome <- cut(dat_fin$wellbeing_outcome, 
#                                  breaks = c(0, 25, 50, 75, 100), 
#                                  labels = c("Q1 (worst)", "Q2", "Q3", "Q4 (best)"))


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


#### outcomes with needs ####
outcomes_needs <- dat_fin |> 
  select(all_of(c(needs, outcomes))) |> 
  drop_na()

outcomes_needs$mentalHealth_outcome <- factor(outcomes_needs$mentalHealth_outcome, 
                                              levels = c(0, 1, 2, 3, 4))
colnames(outcomes_needs) <- c(needs_names, outcomes_names)


outcomes_needs_2 <- outcomes_needs

outcomes_needs_2$cluster1 <- as.factor((ifelse(as.numeric(outcomes_needs$housing)!=1, 2, 1) + 
                                          ifelse(as.numeric(outcomes_needs$basics)!=1, 2, 1) + 
                                          ifelse(as.numeric(outcomes_needs$food)!=1, 2, 1)))

outcomes_needs_2$cluster2 <- as.factor((ifelse(as.numeric(outcomes_needs$work)!=1, 2, 1) + 
                                          ifelse(as.numeric(outcomes_needs$education)!=1, 2, 1)))

outcomes_needs_2$cluster3 <- as.factor((ifelse(as.numeric(outcomes_needs$discrimination)!=1, 2, 1) + 
                                          ifelse(as.numeric(outcomes_needs$isolation)!=1, 2, 1) + 
                                          ifelse(as.numeric(outcomes_needs$loneliness)!=1, 2, 1)))

outcomes_needs_2$childcare <- as.factor((ifelse(as.numeric(outcomes_needs$childcare)!=1, 2, 1)))

outcomes_needs_2$caregiving <- as.factor((ifelse(as.numeric(outcomes_needs$caregiving)!=1, 2, 1)))

outcomes_needs_2$legal <- as.factor((ifelse(as.numeric(outcomes_needs$legal)!=1, 2, 1)))

outcomes_needs_2$transport <- as.factor((ifelse(as.numeric(outcomes_needs$transport)!=1, 2, 1)))

outcomes_needs_2$internet <- as.factor((ifelse(as.numeric(outcomes_needs$internet)!=1, 2, 1)))


outcomes_needs_2 <- outcomes_needs_2 |> select(c("wellbeing", "mental health", "cluster1", "cluster2", "cluster3", "childcare", "caregiving", "legal", "transport", "internet"))

bl_outcomes <- data.frame(
  from = c("wellbeing", "mental health"),
  to   = c("mental health", "wellbeing")
)

bl_needs = data.frame(
  from = c(rep("mental health", times = 8), rep("wellbeing", times = 8)), 
  to = c(paste0("cluster", 1:3), c("childcare", "caregiving", "legal", "transport", "internet"), paste0("cluster", 1:3), c("childcare", "caregiving", "legal", "transport", "internet"))
)

overall_dag <- plot(hc(outcomes_needs_2, blacklist = rbind(bl_outcomes, bl_needs)))

cluster1_dag <- plot(hc(data.frame(outcomes_needs |> select(c("housing", "basics", "food")))))
cluster2_dag <- plot(hc(data.frame(outcomes_needs |> select(c("work", "education")))))
cluster3_dag <- plot(hc(data.frame(outcomes_needs |> select(c("discrimination", "isolation", "loneliness")))))
