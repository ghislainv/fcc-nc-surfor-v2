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
# Only pixels in 2008 en 2021_fnf have three interpreters (A, O, and S).
# Pixels in 2008 and 2021_fnf have either 1 or 3 interpreters (not 2).

# Load datasets for change
file_2000 <- "38354-samples_2000.csv"
file_2008 <- "38356-samples_2010.csv"
file_2021 <- "38355-samples_2021.csv"
df_2000 <- read_csv(here("data_raw", file_2000), show_col_types=FALSE)
df_2008 <- read_csv(here("data_raw", file_2008), show_col_types=FALSE)
df_2021 <- read_csv(here("data_raw", file_2021), show_col_types=FALSE) |>
    mutate(plotid=plotid + 1000) |>
    mutate(sampleid=sampleid + 16000)

# Remove stagiaire observations
# and select only pixels were both interpreters agree for 2008
select_obs <- FALSE
if (select_obs == TRUE) {
  suffix <- "selobs"
  df_2000 <- df_2000 |>
    filter(email != "stagiaire@oeil.nc")
  
  df_2008 <- df_2008 |>
    filter(email != "stagiaire@oeil.nc") |>
    group_by(plotid, sampleid) |>
    # Select only pixels were both interpreters agree
    filter(n_distinct(Classification) == 1) |>
    slice(1) |>
    ungroup()
  
  df_2021 <- df_2021 |>
    filter(email != "stagiaire@oeil.nc")
} else {
  suffix <- "allobs"
  set.seed(1234)
  df_2000 <- df_2000 |>
    group_by(plotid, sampleid) |>
    # Select one observation at random
    slice_sample(n=1) |>
    ungroup()

  df_2008 <- df_2008 |>
    group_by(plotid, sampleid) |>
    # Majority of the observations per pixel
    filter(Classification == names(which.max(table(Classification)))) |>
    slice(1) |>
    ungroup()

  df_2021 <- df_2021 |>
    group_by(plotid, sampleid) |>
    # Select one observation at random
    slice_sample(n=1) |>
    ungroup()
}

# Combining datasets
# Start with 2008 which is the central year and includes
# all plots used for interpreting change.
df_change <- df_2008 |>
  left_join(df_2000, by=c("plotid", "sampleid"), suffix=c("", "_t1")) |>
  left_join(df_2021, by=c("plotid", "sampleid"), suffix=c("", "_t3")) |>
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
file_for2021_fnf <- "raster-intersect_ceo-38675-samples_F-NF-2021.csv"
df_for2021_fnf <- read_csv(here("data_raw", file_for2021_fnf), show_col_types=FALSE) |>
  select(-c(17:22)) |>
  mutate(dataset="for2021") |>
  group_by(plotid, sampleid) |>
  # Majority of the observations per pixel
  filter(Classification == names(which.max(table(Classification)))) |>
  slice(1) |>
  ungroup() |>
  rename(for2021_fnf=Classification)

# Combine the two types of data (change and cover)
df_surfor <- df_change |>
  mutate(dataset="change") |>
  mutate(collection_time=as.character(collection_time)) |>
  bind_rows(df_for2021_fnf) |>
  write_csv(here("data", glue("df_surfor_{suffix}.csv")))
  
# Run the python script to get GEE data for GFC and TFM
# system("python ../python/get-gee-data.py")

# Combining data-sets
df_gee <- read_csv(here("data", glue("extract_gfc_tmf_{suffix}.csv")), show_col_types=FALSE) |>
  select(-lat, -lon)
df_surfor_gee <- read_csv(here("data", glue("df_surfor_{suffix}.csv")), show_col_types=FALSE) |>
  bind_cols(df_gee) |>
  mutate(for2000 = factor(recode_values(
    for2000,
    "Forêt" ~ "Forest",
    "Non Forêt" ~ "NonForest",
    default = NA
  ))) |>
  mutate(for2008 = factor(recode_values(
    for2008,
    "Forêt" ~ "Forest",
    "Non Forêt" ~ "NonForest",
    default = NA
  ))) |>
  mutate(for2021 = factor(recode_values(
    for2021,
    "Forêt" ~ "Forest",
    "Non Forêt" ~ "NonForest",
    default = NA
  ))) |>
  mutate(for2021_fnf = factor(recode_values(
    for2021_fnf,
    "Forêt" ~ "Forest",
    "Non Forêt" ~ "NonForest",
    default = NA
  )))

# Select variables and write
df_surfor_gee <- df_surfor_gee |>
  select(
    plotid, sampleid, sample_internal_id, lon, lat, imagery_title,
    email, dataset,
    for2000, for2008, for2021,
    fcc_p1, fcc_p2, for2021_fnf,
    Dec1999, Dec2000, Dec2007, Dec2008, Dec2020, Dec2021,
    gain, lossyear, treecover2000) |>
  write_csv(here("data", glue("df_surfor_gee_{suffix}.csv")))

# End
