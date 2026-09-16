# Libraries
library(readr)
library(glue)
library(dplyr)
library(here)
library(caret)  # Confusion matrices

# NOTES:
# Between 2000 and 2008, plotid and sampleid correspond between dates.
# But between 2008 and 2010, it is not the case:
# - plotid in 2008 correspond to plotid in 2021 + 1000.
# - sampleid in 2008 correspond to sampleid in 2021 + 16000.

# Load datasets for change
file_2000 <- "38354-samples_2000.csv"
file_2008 <- "38356-samples_2010.csv"
file_2021 <- "38355-samples_2021.csv"
df_2000 <- read_csv(here("data_raw", file_2000), show_col_types=FALSE)
df_2008 <- read_csv(here("data_raw", file_2008), show_col_types=FALSE)
df_2021 <- read_csv(here("data_raw", file_2021), show_col_types=FALSE)

# Majority of the observations per pixel
df_2000 <- df_2000 |>
  group_by(plotid, sampleid) |>
  filter(Classification == names(which.max(table(Classification)))) |>
  slice(1) |>
  ungroup()

df_2008 <- df_2008 |>
  group_by(plotid, sampleid) |>
  filter(Classification == names(which.max(table(Classification)))) |>
  slice(1) |>
  ungroup()

df_2021 <- df_2021 |>
  group_by(plotid, sampleid) |>
  filter(Classification == names(which.max(table(Classification)))) |>
  slice(1) |>
  ungroup() |>
  mutate(plotid=plotid + 1000) |>
  mutate(sampleid=sampleid + 16000)

# Combining datasets
# Start with 2008 which is the central year and includes
# all plots used for interpreting change.
df_change <- df_2008 |>
  full_join(df_2000, by=c("plotid", "sampleid"), suffix=c("", "_t1")) |>
  full_join(df_2021, by=c("plotid", "sampleid"), suffix=c("", "_t3")) |>
  rename(for2008=Classification) |>
  rename(for2000=Classification_t1) |>
  rename(for2021=Classification_t3) |>
  select(!ends_with("_t1")) |>
  select(!ends_with("_t3"))

# Get data on forest cover change
df_change <- df_change |>
  # fcc_p1
  mutate(fcc_p1=ifelse(for2000=="Forêt" & for2008=="Forêt", "stableF", NA)) |>
  mutate(fcc_p1=ifelse(for2000=="Non Forêt" & for2008=="Non Forêt", "stableNF", fcc_p1)) |>
  mutate(fcc_p1=ifelse(for2000=="Forêt" & for2008=="Non Forêt", "loss", fcc_p1)) |>
  mutate(fcc_p1=ifelse(for2000=="Non Forêt" & for2008=="Forêt", "gain", fcc_p1)) |>
  # fcc_p2
  mutate(fcc_p2=ifelse(for2008=="Forêt" & for2021=="Forêt", "stableF", NA)) |>
  mutate(fcc_p2=ifelse(for2008=="Non Forêt" & for2021=="Non Forêt", "stableNF", fcc_p2)) |>
  mutate(fcc_p2=ifelse(for2008=="Forêt" & for2021=="Non Forêt", "loss", fcc_p2)) |>
  mutate(fcc_p2=ifelse(for2008=="Non Forêt" & for2021=="Forêt", "gain", fcc_p2))

# Get data on forest cover in 2021
file_for2021 <- "raster-intersect_ceo-38675-samples_F-NF-2021.csv"
df_for2021 <- read_csv(here("data_raw", file_for2021), show_col_types=FALSE) |>
  select(-c(17:22)) |>
  mutate(dataset="for2021") |>
  rename(for2021_fnf=Classification)

# Combine the two types of data (change and cover)
df_surfor <- df_change |>
  mutate(dataset="change") |>
  mutate(collection_time=as.character(collection_time)) |>
  bind_rows(df_for2021) |>
  write_csv(here("data", "df_surfor.csv"))
  
# Run the python script to get GEE data for GFC and TFM
# system("python ../python/get-gee-data.py")

# Combining data-sets
df_gee <- read_csv(here("data", "extracting_gfc_tmf.csv"), show_col_types=FALSE) |>
  select(-lat, -lon)
df_surfor_gee <- read_csv(here("data", "df_surfor.csv"), show_col_types=FALSE) |>
  bind_cols(df_gee) |>
  mutate(for2000 = factor(case_match(
    for2000,
    "Forêt" ~ "Forest",
    "Non Forêt" ~ "NonForest",
    .default = NA
  ))) |>
  mutate(for2008 = factor(case_match(
    for2008,
    "Forêt" ~ "Forest",
    "Non Forêt" ~ "NonForest",
    .default = NA
  ))) |>
  mutate(for2021 = factor(case_match(
    for2021,
    "Forêt" ~ "Forest",
    "Non Forêt" ~ "NonForest",
    .default = NA
  ))) |>
  mutate(for2021_fnf = factor(case_match(
    for2021_fnf,
    "Forêt" ~ "Forest",
    "Non Forêt" ~ "NonForest",
    .default = NA
  ))) |>
  select(
    plotid, sampleid, sample_internal_id, lon, lat, imagery_title,
    for2000, for2008, for2021,
    fcc_p1, fcc_p2, for2021_fnf,
    Dec1999, Dec2000, Dec2007, Dec2008, Dec2020, Dec2021,
    gain, lossyear, treecover2000) |>
  write_csv(here("data", "df_surfor_gee.csv"))

# End
