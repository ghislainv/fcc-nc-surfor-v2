# Libraries
library(readr)
library(glue)
library(dplyr)
library(here)
library(caret)  # Confusion matrices

# =============
# Load datasets
# =============

# Change on period P1 2000--2008
file_change_p1 <- "raster-intersect_ceo-44650-samples_GP-2000-2010.csv"
file_start_p1 <- "raster-intersect_ceo-44650-samples_GP-2000-2010_v2_FnF.csv"
df_change_p1 <- read_csv(here("data_raw", file_change_p1), show_col_types=FALSE)
df_start_p1 <- read_csv(here("data_raw", file_start_p1), show_col_types=FALSE)

# Correcting variable names
names(df_change_p1) <- gsub("_ ", "_", names(df_change_p1))

# Join tables
df_p1 <- df_change_p1 |>
  dplyr::left_join(df_start_p1, by="sampleid", keep=FALSE) |>
  rename(plotid=plotid.x) |>
  select(-plotid.y) |>
  mutate(fcc_obs=Classification) |>
  mutate(fcc_obs=ifelse(Classification=="Stable" & Classification_MAJ=="Forêt",
                        "StableF", fcc_obs)) |>
  mutate(fcc_obs=ifelse(Classification=="Stable" & Classification_MAJ=="Non Forêt",
                        "StableNF", fcc_obs)) |>
  mutate(fcc_obs=ifelse(Classification=="Stable" & Classification_MAJ=="Non interprétable",
                        "NotInt", fcc_obs)) |>
  mutate(fcc_obs=ifelse(Classification=="Non interprétable",
                        "NotInt", fcc_obs)) |>
  mutate(fcc_obs=ifelse(fcc_obs=="Perte",
                        "Loss", fcc_obs)) |>
  mutate(year_start=2000, year_end=2008) |>
  write_csv(here("data", "df_p1.csv"))

# Change on period P2 2008--2021
file_change_p2 <- "raster-intersect_ceo-45022-samples_GP-2010-2020.csv"
file_start_p2 <- "raster-intersect_ceo-45022-samples_GP-2010-2020_v2_FnF.csv"
df_change_p2 <- read_csv(here("data_raw", file_change_p2), show_col_types=FALSE)
df_start_p2 <- read_csv(here("data_raw", file_start_p2), show_col_types=FALSE)

# Correcting variable names
names(df_change_p2) <- gsub("_ ", "_", names(df_change_p2))

# Join tables
df_p2 <- df_change_p2 |>
  dplyr::left_join(df_start_p2, by="sampleid", keep=FALSE) |>
  rename(plotid=plotid.x) |>
  select(-plotid.y) |>
  mutate(fcc_obs=Classification) |>
  mutate(fcc_obs=ifelse(Classification=="Stable" & Classification_MAJ=="Forêt",
                        "StableF", fcc_obs)) |>
  mutate(fcc_obs=ifelse(Classification=="Stable" & Classification_MAJ=="Non Forêt",
                        "StableNF", fcc_obs)) |>
  mutate(fcc_obs=ifelse(Classification=="Stable" & Classification_MAJ=="Non interprétable",
                        "NotInt", fcc_obs)) |>
  mutate(fcc_obs=ifelse(Classification=="Non inerprétable",
                        "NotInt", fcc_obs)) |>
  mutate(fcc_obs=ifelse(fcc_obs=="Perte",
                        "Loss", fcc_obs)) |>
  mutate(year_start=2008, year_end=2021) |>
  write_csv(here("data", "df_p2.csv"))

# Data on forest cover in 2021
file_for2021 <- "raster-intersect_ceo-38675-samples_F-NF-2021.csv"
df_for2021 <- read_csv(here("data_raw", file_for2021), show_col_types=FALSE)
names(df_for2021) <- gsub("_ ", "_", names(df_for2021))
df_for2021 <- df_for2021 |>
  mutate(for2021=Classification) |>
  mutate(for2021=ifelse(Classification=="Non interprétable", "NotInt", for2021)) |>
  mutate(for2021=ifelse(Classification=="Forêt", "Forest", for2021)) |>
  mutate(for2021=ifelse(Classification=="Non Forêt", "NonForest", for2021)) |>
  write_csv(here("data", "df_for2021.csv"))

# ========================================================
# Best percentage of tree cover for GFC based on year 2021
# ========================================================

# ========================================================
# Confusion matrix for forest cover in 2000, 2008 and 2021
# ========================================================

# Observed forest p1
df_p1 <- df_p1 |>
  mutate(for2000_obs=ifelse(fcc_obs=="StableF" | fcc_obs=="Loss", "Forest", fcc_obs)) |>
  mutate(for2000_obs=ifelse(fcc_obs=="StableNF" | fcc_obs=="Gain", "NonForest", for2000_obs)) |>
  mutate(for2000_obs=ifelse(fcc_obs=="NotInt", NA, for2000_obs))
df_p1 <- df_p1 |>
  mutate(for2008_obs=ifelse(fcc_obs=="StableF" | fcc_obs=="Gain", "Forest", fcc_obs)) |>
  mutate(for2008_obs=ifelse(fcc_obs=="StableNF" | fcc_obs=="Loss", "NonForest", for2008_obs)) |>
  mutate(for2008_obs=ifelse(fcc_obs=="NotInt", NA, for2008_obs))

# Observed forest p2
df_p2 <- df_p2 |>
  mutate(for2008_obs=ifelse(fcc_obs=="StableF" | fcc_obs=="Loss", "Forest", fcc_obs)) |>
  mutate(for2008_obs=ifelse(fcc_obs=="StableNF" | fcc_obs=="Gain", "NonForest", for2008_obs)) |>
  mutate(for2008_obs=ifelse(fcc_obs=="NotInt", NA, for2008_obs))
df_p2 <- df_p2 |>
  mutate(for2021_obs=ifelse(fcc_obs=="StableF" | fcc_obs=="Gain", "Forest", fcc_obs)) |>
  mutate(for2021_obs=ifelse(fcc_obs=="StableNF" | fcc_obs=="Loss", "NonForest", for2021_obs)) |>
  mutate(for2021_obs=ifelse(fcc_obs=="NotInt", NA, for2021_obs))
  
# Predicted forest
df_p1 <- df_p1 |>
  mutate(for2000_pred_gfc_60=ifelse((gfc_treecover2000 >= 60), "Forest", "NonForest")) |>
  mutate(for2008_pred_gfc_60=ifelse((gfc_treecover2000 >= 60) & !(gfc_lossyear %in% c(1:8)), "Forest", "NonForest"))
df_p2 <- df_p2 |>
  mutate(for2008_pred_gfc_60=ifelse((gfc_treecover2000 >= 60) & !(gfc_lossyear %in% c(1:8)), "Forest", "NonForest")) |> # 8 for forest at end of year 2008
  mutate(for2021_pred_gfc_60=ifelse((gfc_treecover2000 >= 60) & !(gfc_lossyear %in% c(1:21)), "Forest", "NonForest"))  # 21 for forest at end of year 2021

# Observed and predicted forest for2021
df_for2021 <- df_for2021 |>
  mutate(for2021_obs=ifelse(for2021=="NotInt", NA, for2021)) |>
  mutate(for2021_pred_gfc_60=ifelse((gfc_treecover2000 >= 60) & !(gfc_lossyear %in% c(1:21)), "Forest", "NonForest")) |> # 21 for forest at end of year 2021
  mutate(for2021_pred_tmf=ifelse(tmf_ac_Dec2021 %in% c(1, 2, 4), "Forest", "NonForest"))

# 2000
cmat_2000 <- df_p1 |>
  select(1:6, for2000_obs, for2000_pred_gfc_60) |>
  mutate(for2000_obs=as.factor(for2000_obs)) |>
  mutate(for2000_pred_gfc_60=as.factor(for2000_pred_gfc_60)) |>
  na.omit()
conf_mat_2000 <- caret::confusionMatrix(
  data=cmat_2000$for2000_pred_gfc_60,
  reference=cmat_2000$for2000_obs)
print(conf_mat_2000)

# 2008
cmat_2008 <- data.frame(
  obs=c(df_p1$for2008_obs, df_p2$for2008_obs),
  pred_gfc_60=c(df_p1$for2008_pred_gfc_60, df_p2$for2008_pred_gfc_60))
conf_mat_2008 <- caret::confusionMatrix(
  data=as.factor(cmat_2008$pred_gfc_60),
  reference=as.factor(cmat_2008$obs))
print(conf_mat_2008)

# 2021 at the pixel level
# gfc
conf_mat_2021_gfc <- caret::confusionMatrix(
  data=as.factor(df_for2021$for2021_obs),
  reference=as.factor(df_for2021$for2021_pred_gfc_60))
# tmf
conf_mat_2021_tmf <- caret::confusionMatrix(
  data=as.factor(df_for2021$for2021_obs),
  reference=as.factor(df_for2021$for2021_pred_tmf))
