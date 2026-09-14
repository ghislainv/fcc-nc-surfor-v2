# Libraries
library(readr)
library(glue)
library(dplyr)
library(here)
library(caret)  # Confusion matrices

# =============
# Load datasets
# =============

file_2000 <- "38354-samples_2000.csv"
file_2008 <- "38356-samples_2010.csv"
file_2021 <- "38355-samples_2021.csv"
# 5 decimals for degrees give a precision of about 1 m
df_2000 <- read_csv(here("data_raw", file_2000), show_col_types=FALSE) |>
  mutate(latitude=round(lat, 5), longitude=round(lon, 5))
df_2008 <- read_csv(here("data_raw", file_2008), show_col_types=FALSE) |>
  mutate(latitude=round(lat, 5), longitude=round(lon, 5))
df_2021 <- read_csv(here("data_raw", file_2021), show_col_types=FALSE) |>
  mutate(latitude=round(lat, 5), longitude=round(lon, 5))

# Majority of the observations per pixel
df_2000 <- df_2000 |>
  group_by(longitude, latitude) |>
  filter(Classification == names(which.max(table(Classification)))) |>
  slice(1) |>
  ungroup()

df_2008 <- df_2008 |>
  group_by(longitude, latitude) |>
  filter(Classification == names(which.max(table(Classification)))) |>
  slice(1) |>
  ungroup()

df_2021 <- df_2021 |>
  group_by(longitude, latitude) |>
  filter(Classification == names(which.max(table(Classification)))) |>
  slice(1) |>
  ungroup()

# Combine datasets
df_year <- df_2008 |>
  left_join(df_2000, by=c("longitude", "latitude"), suffix=c("", "_t1")) |>
  rename(Classification_2008=Classification) |>
  rename(Classification_2000=Classification_t1) |>
  left_join(df_2021, by=c("longitude", "latitude"), suffix=c("", "_t3")) |>
  rename(Classification_2021=Classification) |>
  select(!ends_with("_t1")) |>
  select(!ends_with("_t3")) |>
  write_csv(here("data", "df_year.csv"))

# Get change
df_year <- df_year |>
  # fcc_p1
  mutate(fcc_p1=ifelse(Classification_2000=="Forêt" & Classification_2008=="Forêt", "stableF", NA)) |>
  mutate(fcc_p1=ifelse(Classification_2000=="Non Forêt" & Classification_2008=="Non Forêt", "stableNF", fcc_p1)) |>
  mutate(fcc_p1=ifelse(Classification_2000=="Forêt" & Classification_2008=="Non Forêt", "loss", fcc_p1)) |>
  mutate(fcc_p1=ifelse(Classification_2000=="Non Forêt" & Classification_2008=="Forêt", "gain", fcc_p1)) |>
  # fcc_p2
  mutate(fcc_p2=ifelse(Classification_2008=="Forêt" & Classification_2021=="Forêt", "stableF", NA)) |>
  mutate(fcc_p2=ifelse(Classification_2008=="Non Forêt" & Classification_2021=="Non Forêt", "stableNF", fcc_p2)) |>
  mutate(fcc_p2=ifelse(Classification_2008=="Forêt" & Classification_2021=="Non Forêt", "loss", fcc_p2)) |>
  mutate(fcc_p2=ifelse(Classification_2008=="Non Forêt" & Classification_2021=="Forêt", "gain", fcc_p2)) |>
  write_csv(here("data", "df_year.csv"))
