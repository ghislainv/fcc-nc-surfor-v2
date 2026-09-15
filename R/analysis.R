# Libraries
library(readr)
library(glue)
library(dplyr)
library(here)
library(caret)  # Confusion matrices

# Load data
df <- read_csv(here("data", "df_surfor_gee.csv"), show_col_types=FALSE)

# ========================================================
# Best percentage of tree cover for GFC based on year 2021
# ========================================================

# ========================================================
# Confusion matrix for forest cover in 2000, 2008 and 2021
# ========================================================

# Predicted forest from GFC (considering tree cover >= 60%)
df_cm <- df |>
  mutate(for2000_pred_gfc_60=ifelse((treecover2000 >= 60), "Forest", "NonForest")) |>
  mutate(for2008_pred_gfc_60=ifelse((treecover2000 >= 60) & !(lossyear %in% c(1:7)), "Forest", "NonForest")) |> # 7 for forest at end of year 2007
  mutate(for2021_pred_gfc_60=ifelse((treecover2000 >= 60) & !(lossyear %in% c(1:20)), "Forest", "NonForest"))  # 20 for forest at end of year 2020

# Predicted forest from TMF
df_cm <- df_cm |>
  mutate(for2000_pred_tmf=ifelse(Dec1999 %in% c(1, 2, 4), "Forest", "NonForest")) |>
  mutate(for2008_pred_tmf=ifelse(Dec2007 %in% c(1, 2, 4), "Forest", "NonForest")) |>
  mutate(for2021_pred_tmf=ifelse(Dec2020 %in% c(1, 2, 4), "Forest", "NonForest"))

# 2000
cmat_2000 <- df_cm |>
  mutate(for2000_obs=as.factor(for2000)) |>
  mutate(for2000_pred_gfc_60=as.factor(for2000_pred_gfc_60)) |>
  select(for2000_obs, for2000_pred_gfc_60) |>
  na.omit()
conf_mat_2000 <- caret::confusionMatrix(
  data=cmat_2000$for2000_pred_gfc_60,
  reference=cmat_2000$for2000_obs)
print(conf_mat_2000)

# 2008
cmat_2008 <- df_cm |>
  mutate(for2008_obs=as.factor(for2008)) |>
  mutate(for2008_pred_gfc_60=as.factor(for2008_pred_gfc_60)) |>
  select(for2008_obs, for2008_pred_gfc_60) |>
  na.omit()
conf_mat_2008 <- caret::confusionMatrix(
  data=cmat_2008$for2008_pred_gfc_60,
  reference=cmat_2008$for2008_obs)
print(conf_mat_2008)

# 2021
cmat_2021 <- df_cm |>
  mutate(for2021_obs=as.factor(for2021)) |>
  mutate(for2021_pred_gfc_60=as.factor(for2021_pred_gfc_60)) |>
  select(for2021_obs, for2021_pred_gfc_60) |>
  na.omit()
conf_mat_2021 <- caret::confusionMatrix(
  data=cmat_2021$for2021_pred_gfc_60,
  reference=cmat_2021$for2021_obs)
print(conf_mat_2021)

# fnf2021
cmat_2021_fnf <- df_cm |>
  mutate(for2021_obs=as.factor(for2021_fnf)) |> # ! for2021fnf here
  mutate(for2021_pred_gfc_60=as.factor(for2021_pred_gfc_60)) |>
  select(for2021_obs, for2021_pred_gfc_60) |>
  na.omit()
conf_mat_2021_fnf <- caret::confusionMatrix(
  data=cmat_2021_fnf$for2021_pred_gfc_60,
  reference=cmat_2021_fnf$for2021_obs)
print(conf_mat_2021_fnf)
