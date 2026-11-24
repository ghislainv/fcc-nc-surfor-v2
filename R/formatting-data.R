# Libraries
library(readr)
library(glue)
library(dplyr)
library(here)
library(caret)  # Confusion matrices

# NOTES:
# Between 2004 and 2008, plotid and sampleid correspond between dates.
# But between 2008 and 2021, it is not the case:
# - plotid in 2008 correspond to plotid in 2021 + 1000.
# - sampleid in 2008 correspond to sampleid in 2021 + 16000.
# Only pixels in 2008 en 2021_fnf have three interpreters (A, O, and S).
# Pixels in 2008 and 2021_fnf have either 1 or 3 interpreters (not 2).
  
# Function "format_data"
format_data <- function() {

  # -------------------
  # Message
  # -------------------

  cat("Formatting data\n")

  # ------------------------
  # Load datasets for change
  # ------------------------
  
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
  
  # ---------------------------------------------
  # Define consensus and Not Interpretable pixels
  # ---------------------------------------------

  set.seed(1234) # Set seed for reproducibility of the results
  df_2004 <- df_2004 |>
    group_by(plotid, sampleid) |>
    # Set "Non interprétable" if two different interpretations
    mutate(Classification = ifelse(length(unique(Classification)) == 2,
                                   "Non interprétable", Classification)) |>
    # Majority of the observations per pixel
    filter(Classification == names(which.max(table(Classification)))) |>
    slice(1) |>
    ungroup()

  df_2008 <- df_2008 |>
    group_by(plotid, sampleid) |>
    # Set "Non interprétable" if three different interpretations (in fact it is never the case)
    mutate(Commentaire = paste("agreement:", max(table(Classification)),
                               "/", length(Classification))) |>
    mutate(Classification = ifelse(length(unique(Classification)) == 3,
                                   "Non interprétable", Classification)) |> 
    # Majority of the observations per pixel
    filter(Classification == names(which.max(table(Classification)))) |>
    slice(1) |>
    ungroup()

  df_2021 <- df_2021 |>
    group_by(plotid, sampleid) |>
    # Set "Non interprétable" if two different interpretations
    mutate(Classification = ifelse(length(unique(Classification)) == 2,
                                   "Non interprétable", Classification)) |>
    # Majority of the observations per pixel
    filter(Classification == names(which.max(table(Classification)))) |>
    slice(1) |>
    ungroup()

  # ------------------
  # Combining datasets
  # ------------------
  # Start with 2008 which is the central year and includes
  # all plots used for interpreting change.
  df_change <- df_2008 |>
    full_join(df_2004, by=c("plotid", "sampleid", "dataset"), suffix=c("", "_t1")) |>
    full_join(df_2021, by=c("plotid", "sampleid", "dataset"), suffix=c("", "_t3")) |>
    rename(for2008=Classification) |>
    rename(for2004=Classification_t1) |>
    rename(for2021=Classification_t3) |>
    rename(email2008=email) |>
    rename(email2004=email_t1) |>
    rename(email2021=email_t3) |>
    select(!ends_with("_t1")) |>
    select(!ends_with("_t3"))

  # Estimate forest cover change from forest cover at two dates
  df_change <- df_change |>
    # fcc_p1
    mutate(fcc_p1=ifelse(for2004=="Forêt" & for2008=="Forêt", "stableF", NA)) |>
    mutate(fcc_p1=ifelse(for2004=="Non Forêt" & for2008=="Non Forêt", "stableNF", fcc_p1)) |>
    mutate(fcc_p1=ifelse(for2004=="Forêt" & for2008=="Non Forêt", "loss", fcc_p1)) |>
    mutate(fcc_p1=ifelse(for2004=="Non Forêt" & for2008=="Forêt", "gain", fcc_p1)) |>
    # fcc_p2
    mutate(fcc_p2=ifelse(for2008=="Forêt" & for2021=="Forêt", "stableF", NA)) |>
    mutate(fcc_p2=ifelse(for2008=="Non Forêt" & for2021=="Non Forêt", "stableNF", fcc_p2)) |>
    mutate(fcc_p2=ifelse(for2008=="Forêt" & for2021=="Non Forêt", "loss", fcc_p2)) |>
    mutate(fcc_p2=ifelse(for2008=="Non Forêt" & for2021=="Forêt", "gain", fcc_p2)) |>
    # Change time as character
    mutate(collection_time=as.character(collection_time))

  # Get data on forest cover in 2021
  file_for2021_fnf <- "raster-intersect_ceo-38675-samples_F-NF-2021.csv"
  df_for2021_fnf <- read_csv(here("data_raw", file_for2021_fnf), show_col_types=FALSE) |>
    select(-c(17:22)) |>
    mutate(dataset="for2021") |>
    rename(email2021_fnf=email) |>
    group_by(plotid, sampleid) |>
    # Set "Non interprétable" if three different interpretations (in fact it is never the case)
    mutate(Commentaire = paste("agreement:", max(table(Classification)),
                               "/", length(Classification))) |>
    mutate(Classification = ifelse(length(unique(Classification)) == 3,
                                   "Non interprétable", Classification)) |> 
    # Majority of the observations per pixel
    filter(Classification == names(which.max(table(Classification)))) |>
    slice(1) |>
    ungroup() |>
    rename(for2021_fnf=Classification)

  # Combine the two types of data (change and cover)
  df_surfor <- df_change |>
    bind_rows(df_for2021_fnf) |>
    select(-Commentaire) |>
    mutate(for2004 = factor(recode_values(
      for2004,
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
    ))) |>
    write_csv(here("data", "df_interp_consensus.csv"))
  
  # Combine with GEE extract
  df_gee <- read_csv(here("data", glue("df_gee.csv")), show_col_types=FALSE) |>
    select(-lat, -lon)
  df_surfor_gee <- read_csv(here("data", "df_interp_consensus.csv"), show_col_types=FALSE) |>
    left_join(df_gee, by=c("plotid", "sampleid", "dataset"))

  # Select variables and write
  df_surfor_gee <- df_surfor_gee |>
    select(
      plotid, sampleid, dataset, lon, lat, imagery_title,
      email2004, email2008, email2021, email2021_fnf,
      for2004, for2008, for2021,
      fcc_p1, fcc_p2, for2021_fnf,
      Dec2003, Dec2004, Dec2007, Dec2008, Dec2020, Dec2021,
      gain, lossyear, treecover2000) |>
    write_csv(here("data", "df_surfor_consensus.csv"))
}

# Use format_data function
format_data()

# End
