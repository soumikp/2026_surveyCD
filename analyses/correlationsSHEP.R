pacman::p_load(here, tidyverse, pheatmap, reshape2)
source(file.path(here(), "analyses", "pullSHEP.R"))

data <- dat_fin |> 
  dplyr::select(need_vars, risk_vars, WEIGHT, WBS_total, RateHealth_T1, RateMentalHealth_T1, SURVID) |> 
  drop_na()

needs <- tibble(name = need_vars, y = "n") |> 
  slice_head(n=-1) |> 
  mutate(newname = str_remove(paste0(y, "_", name), "_T1"))

risks <- tibble(name = risk_vars, y = "r") |> 
  slice_head(n=-1) |> 
  mutate(newname = str_remove(paste0(y, "_", name), "_T1"))

needs <- data |> 
  dplyr::select(need_vars, SURVID, WEIGHT) |> 
  pivot_longer(cols = -c(SURVID, WEIGHT)) |>
  left_join(needs, by = "name") |> 
  dplyr::select(-c(name, y)) |> 
  pivot_wider(names_from = newname, values_from = value)

risks <- data |> 
  dplyr::select(risk_vars, SURVID, WEIGHT) |> 
  pivot_longer(cols = -c(SURVID, WEIGHT)) |>
  left_join(risks, by = "name") |> 
  dplyr::select(-c(name, y)) |> 
  pivot_wider(names_from = newname, values_from = value)

outcomes <- data |> 
  dplyr::select(SURVID, WEIGHT, RateMentalHealth_T1, WBS_total)

data <- inner_join(needs, outcomes) |> inner_join(risks)

#### some manual cleaning for risks
data$r_GetCaregiving <- (fct_collapse(data$r_GetCaregiving, "1" = c("1", "8")))
data$r_GetChildcare <- (fct_collapse(data$r_GetChildcare, "1" = c("1", "8")))
data$r_FindWork <- (fct_collapse(data$r_FindWork, "1" = c("1", "8")))
for(i in 18:ncol(data)){
  data[,i] <- droplevels(data[,i])
}
data$r_HaveHousing <- factor(5 - as.numeric(data$r_HaveHousing), levels = c(1, 2, 3, 4))
data$r_AccessInternet <- factor(5 - as.numeric(data$r_AccessInternet), levels = c(1, 2, 3, 4))
data$r_LegalNeed <- factor(3 - as.numeric(data$r_LegalNeed), levels = c(1, 2))
#### some manual cleaning for neeeds
data <- data |> mutate(across(starts_with("n"), 
                              ~ ordered(as.factor(.), levels = c("1", "2", "3"))))
#### some manual cleaning for outcomes
data$o_mh <- as.factor(6 - as.numeric(data$RateMentalHealth_T1))
data$o_wb <- cut(data$WBS_total, breaks = c(0, 25, 50, 75, 100), include.lowest = TRUE, 
                      labels = c("1", "2", "3", "4"))
data <- data |> dplyr::select(-c(RateMentalHealth_T1, WBS_total))
data <- data |> dplyr::select(c("SURVID", "WEIGHT", colnames(needs)[3:15], "o_mh", "o_wb", colnames(risks)[3:18]))

data_cor <- lavaan::lavCor(data.frame(data[-c(1)]), 
                           ordered = colnames(data)[-c(1, 2)], 
                           sampling.weights = "WEIGHT")

order <- colnames(data)[-c(1, 2)]

melted_cor <- melt(data_cor)

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
       title = "Survey-weighted polychoric correlations between needs (n), risks (r), and outcomes (o).", 
       subtitle = "n = 5396 complete cases")

factor <- 2
ggsave(file.path(here(), "analyses", "2026_02_24_correlations.pdf"), 
       corplot, 
       height = factor*10, 
       width = factor*10, 
       units = "in", 
       device = "pdf")




data_cor_on <- lavaan::lavCor(data.frame(data[c(2:17)]), 
                           ordered = colnames(data)[c(3:17)], 
                           sampling.weights = "WEIGHT")
cor_map_on <- pheatmap(
  data_cor_on,
  main = "Survey-weighted polychoric correlation between outcomes and needs",
  fontsize = 20,
  display_numbers = TRUE,  # Show the correlation coefficient in each cell
  cluster_rows = FALSE,   # Cluster rows to group similar variables
  cluster_cols = TRUE   # Cluster columns
)

ggsave(file.path(here(), "analyses", "2026_02_24_clustering_outcomeNeed.pdf"), 
       cor_map_on, 
       height = factor*10, 
       width = factor*10, 
       units = "in", 
       device = "pdf")

  