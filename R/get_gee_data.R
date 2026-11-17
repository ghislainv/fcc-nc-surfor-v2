# Libraries
library(readr)
library(dplyr)
library(here)

# NOTES:
# Between 2004 and 2008, plotid and sampleid correspond between dates.
# But between 2008 and 2021, it is not the case:
# - plotid in 2008 correspond to plotid in 2021 + 1000.
# - sampleid in 2008 correspond to sampleid in 2021 + 16000.
# Only pixels in 2008 en 2021_fnf have three interpreters (A, O, and S).
# Pixels in 2008 and 2021_fnf have either 1 or 3 interpreters (not 2).

# Load datasets for change
file_2004 <- "38354-samples_2000.csv"
file_2008 <- "38356-samples_2010.csv"
file_2021 <- "38355-samples_2021.csv"
df_2004 <- read_csv(here("data_raw", file_2004), show_col_types=FALSE) |>
  mutate(year=2004, dataset="change_0408")
df_2008 <- read_csv(here("data_raw", file_2008), show_col_types=FALSE) |>
  mutate(year=2008)
df_2021 <- read_csv(here("data_raw", file_2021), show_col_types=FALSE) |>
  mutate(plotid=plotid + 1000) |>
  mutate(sampleid=sampleid + 16000) |>
  mutate(year=2021, dataset="change_0821")

# Add the period to observations in 2008 (the pivotal year including all plots)
id_2004 <- paste(df_2004$plotid, df_2004$sampleid, sep="_")
id_2008 <- paste(df_2008$plotid, df_2008$sampleid, sep="_")
df_2008 <- df_2008 |>
  mutate(dataset=ifelse(id_2008 %in% id_2004, "change_0408", "change_0821"))

df_2004 <- df_2004 |>
  group_by(plotid, sampleid) |>
  # Select one observation
  slice(1) |>
  ungroup()

df_2008 <- df_2008 |>
  group_by(plotid, sampleid) |>
  slice(1) |>
  ungroup()

df_2021 <- df_2021 |>
  group_by(plotid, sampleid) |>
  slice(1) |>
  ungroup()

# Combining datasets
# Start with 2008 which is the central year and includes
# all plots used for interpreting change.
df_change <- df_2008 |>
  left_join(df_2004, by=c("plotid", "sampleid"), suffix=c("", "_t1")) |>
  left_join(df_2021, by=c("plotid", "sampleid"), suffix=c("", "_t3")) |>
  mutate(collection_time=as.character(collection_time))

# Get data on forest cover in 2021
file_for2021_fnf <- "raster-intersect_ceo-38675-samples_F-NF-2021.csv"
df_for2021_fnf <- read_csv(here("data_raw", file_for2021_fnf), show_col_types=FALSE) |>
  select(-c(17:22)) |>
  mutate(dataset="for2021") |>
  group_by(plotid, sampleid) |>
  slice(1) |>
  ungroup()

# Combine the two types of data (change and cover)
df_surfor <- df_change |>
  bind_rows(df_for2021_fnf) |>
  select(plotid, sampleid, sample_internal_id, lon, lat, dataset) |>
  write_csv(here("data", "df_points_for_gee.csv"))
  
# Run the python script to get GEE data for GFC and TFM
# system("python ../python/get-gee-data.py")

# Combining data-sets
extract_gee <- read_csv(here("data", "extract_gee_gfc_tmf.csv"), show_col_types=FALSE) |>
  select(-lat, -lon)
df_gee <- read_csv(here("data", "df_points_for_gee.csv"), show_col_types=FALSE) |>
  bind_cols(extract_gee)

# Select variables and write
df_gee <- df_gee |>
  select(
    plotid, sampleid, sample_internal_id, lon, lat, dataset,
    Dec2003, Dec2004, Dec2007, Dec2008, Dec2020, Dec2021,
    gain, lossyear, treecover2000) |>
  write_csv(here("data", "df_gee.csv"))

# End
