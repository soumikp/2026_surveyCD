pacman::p_load(here, tidyverse, pheatmap, reshape2)
source(file.path(here(), "analyses", "pullSHEP.R"))

data <- dat_fin |> 
  select(need_vars, risk_vars, WEIGHT, WBS_total, RateHealth_T1, RateMentalHealth_T1, SURVID) |> 
  drop_na()

needs <- tibble(name = need_vars, y = "n") |> 
  slice_head(n=-1) |> 
  mutate(newname = str_remove(paste0(y, "_", name), "_T1"))

risks <- tibble(name = risk_vars, y = "r") |> 
  slice_head(n=-1) |> 
  mutate(newname = str_remove(paste0(y, "_", name), "_T1"))

needs <- data |> 
  select(need_vars, SURVID, WEIGHT) |> 
  pivot_longer(cols = -c(SURVID, WEIGHT)) |>
  left_join(needs, by = "name") |> 
  select(-c(name, y)) |> 
  pivot_wider(names_from = newname, values_from = value)

risks <- data |> 
  select(risk_vars, SURVID, WEIGHT) |> 
  pivot_longer(cols = -c(SURVID, WEIGHT)) |>
  left_join(risks, by = "name") |> 
  select(-c(name, y)) |> 
  pivot_wider(names_from = newname, values_from = value)

outcomes <- data |> 
  select(SURVID, WEIGHT, RateMentalHealth_T1, WBS_total)

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
                              ~ fct_collapse(as.factor(.), "2" = c("2", "3"))))
#### some manual cleaning for outcomes
data$o_mh <- as.factor(6 - as.numeric(data$RateMentalHealth_T1))
data$o_wb <- cut(data$WBS_total, breaks = c(0, 25, 50, 75, 100), include.lowest = TRUE, 
                 labels = c("1", "2", "3", "4"))
data <- data |> select(-c(RateMentalHealth_T1, WBS_total))
data <- data |> select(c("SURVID", "WEIGHT", colnames(needs)[3:15], "o_mh", "o_wb", colnames(risks)[3:18]))

op <- OrdCD(y = data.frame(data[str_detect(colnames(data), "n_|o_")]), 
      w = data$WEIGHT, 
      blacklist = expand.grid(from = colnames(data)[str_detect(colnames(data), "o_")], 
                              to = colnames(data)[str_detect(colnames(data), "n_")]))
