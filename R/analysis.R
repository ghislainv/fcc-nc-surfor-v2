# Libraries
library(readr)
library(glue)
library(dplyr)
library(here)
library(caret)  # Confusion matrices
library(ggplot2)
library(sf)
library(tidyr)

# Output directory
dir.create(here("outputs"), showWarnings=FALSE)

# Load data
df <- read_csv(here("data", "df_surfor_gee_rmobs.csv"), show_col_types=FALSE)

# ========================================================
# Disagreement between interpreters
# ========================================================

# Formatting data
# ---------------

# Variables to keep
keep_var <- c("plotid", "sampleid", "email", "dataset", "Classification")

# Load datasets for change
file_2000 <- "38354-samples_2000.csv"
file_2008 <- "38356-samples_2010.csv"
file_2021 <- "38355-samples_2021.csv"
df_2000 <- read_csv(here("data_raw", file_2000), show_col_types=FALSE) |>
  mutate(dataset="for2000") |>
  select(keep_var)
df_2008 <- read_csv(here("data_raw", file_2008), show_col_types=FALSE) |>
  mutate(dataset="for2008") |>
  select(keep_var)
df_2021 <- read_csv(here("data_raw", file_2021), show_col_types=FALSE) |>
  mutate(dataset="for2021") |>
  select(keep_var)

# Load dataset for forest cover 2021
file_for2021 <- "raster-intersect_ceo-38675-samples_F-NF-2021.csv"
df_for2021 <- read_csv(here("data_raw", file_for2021), show_col_types=FALSE) |>
  mutate(dataset="for2021_fnf") |>
  select(keep_var)

# Combine data-sets
df_disag <- df_2000 |>
  bind_rows(df_2008) |>
  bind_rows(df_2021) |>
  bind_rows(df_for2021) |>
  mutate(Classification=ifelse(Classification %in% c("Indéterminé", "Non interprétable"), "NotInt", Classification))

# Number of observations per date and interpreter
df_nobs_inter <- df_disag |>
  summarize(n=n(), .by=c(dataset, email)) |>
  pivot_wider(id_cols=email, names_from=dataset, values_from=n) |>
  write_csv(here("outputs", "df_nobs_interpreter.csv"))
  
# Wide data
data_wide <- df_disag |>
  mutate(Classification = case_when(
    Classification == "Forêt"     ~ 1,
    Classification == "Non Forêt" ~ 2,
    Classification == "NotInt"    ~ 3
  )) |>
  pivot_wider(
    id_cols     = c(plotid, sampleid, dataset),
    names_from  = email,
    values_from = Classification
  )

# Raw disagreement per dataset
# ----------------------------
datasets <- sort(unique(data_wide$dataset))
n_datasets <- length(datasets)
raw_disag <- data.frame(
  dataset=datasets)
ff1 <- function(row){ifelse(sum(is.na(row)) < 2, 1, 0)}
ff2 <- function(row){length(unique(row[!is.na(row)])) != 1}
for (i in 1:n_datasets) {
  df <- data_wide[data_wide$dataset==datasets[i], 4:6]
  raw_disag[i, "n_pix"] <- nrow(df)
  raw_disag[i, "n_pix_rep"] <- sum(apply(df, 1, ff1))
  raw_disag[i, "n_disag"] <- sum(apply(df, 1, ff2))
}
raw_disag[n_datasets + 1, 2:4] <- apply(raw_disag[, 2:4], 2, sum)
raw_disag[n_datasets + 1, 1] <- "All combined"
raw_disag$perc_rep <- round(100 * (raw_disag$n_pix_rep / raw_disag$n_pix), 1)
raw_disag$perc_disag <- round(100 * (raw_disag$n_disag / raw_disag$n_pix_rep), 1)
raw_disag$n_plots <- raw_disag$n_pix / 9
raw_disag <- raw_disag |>
  relocate(n_plots, .before=n_pix) |>
  relocate(perc_rep, .before=n_disag)

# Krippendorff's Alpha
# --------------------
get_alpha <- function(data_wide, dataset_names) {
  # --- Prepare the matrix of interpreters only ---
  mat <- data_wide |>
    filter(dataset %in% dataset_names) |>
    select(-plotid, -sampleid, -dataset) |>
    as.matrix()

  # --- Observed disagreement (D_o) ---
  # For each pixel, count the disagreeing pairs
  do_pixel <- apply(mat, 1, function(row) {
    vals <- na.omit(row)
    n <- length(vals)
    if (n < 2) return(NA)
    # Number of disagreeing pairs over total pairs
    paires_desaccord <- sum(outer(vals, vals, "!=")) / 2
    paires_total     <- n * (n - 1) / 2
    paires_desaccord / paires_total
  })
  D_o <- mean(do_pixel, na.rm = TRUE)

  # --- Expected disagreement by chance (D_e) ---
  # Frequency of each class across all annotations
  toutes_vals <- na.omit(as.vector(mat))
  n_total     <- length(toutes_vals)
  freq        <- table(toutes_vals) / n_total

  # D_e = 1 - sum(freq_k^2)
  # [probability that 2 randomly drawn annotations differ]
  D_e <- 1 - sum(freq^2)

  # --- Alpha ---
  alpha <- 1 - (D_o / D_e)
  return(list(D_o=round(100 * D_o,  1),
              D_e=round(100 * D_e,  1),
              alpha=round(100 * alpha, 1)))
}

# Loop
df_alpha <- data.frame(dataset=c(datasets, "All combined"))
for (i in 1:n_datasets) {
  K_alpha <- get_alpha(data_wide, datasets[i])
  df_alpha[i, "D_o"] <- K_alpha$D_o
  df_alpha[i, "D_e"] <- K_alpha$D_e
  df_alpha[i, "K_alpha"] <- K_alpha$alpha
}
K_alpha <- get_alpha(data_wide, datasets)
df_alpha[n_datasets + 1, 2:4] <- K_alpha

# Combine raw disagreement and alpha
df_disag_interpreters <- raw_disag |>
  left_join(df_alpha, by="dataset") |>
  write_csv(here("outputs", "df_disag_interp.csv"))

# Bad interpreter
# ---------------
# Removing interpreter one by one to see how K_alpha increases

K_alpha_full <- df_alpha$K_alpha[df_alpha$dataset=="All combined"]
df_bad <- data.frame(removed=c("S", "A", "O"))
for (i in 1:3) {
  dw <- data_wide |>
    select(-c(3 + i))
  K_alpha <- get_alpha(dw, datasets)
  df_bad[i, "alpha"] <- K_alpha$alpha
  K_change <- K_alpha$alpha - K_alpha_full
  df_bad[i, "alpha_change"] <- K_change
}
df_bad |>
  write_csv(here("outputs", "df_bad_interpreter.csv"))

# ========================================================
# Best percentage of tree cover for GFC based on year 2021
# ========================================================

# Loop on tree cover
perc <- seq(10, 90, by=10)
mat_acc <- data.frame(perc, acc=NA)
for (i in 1:length(perc)) {
  for2021_pred_gfc <- ifelse((df$treecover2000 >= perc[i]) & !(df$lossyear %in% c(1:20)),
                                "Forest", "NonForest")
  df_cm <- data.frame(
    obs=as.factor(df$for2021_fnf),
    pred=as.factor(for2021_pred_gfc))
  df_cm <- na.omit(df_cm)
  conf_mat <- caret::confusionMatrix(
    data=df_cm$pred,
    reference=df_cm$obs)
  mat_acc$acc[i] <- conf_mat$overall["Accuracy"] 
}

# Save results
mat_acc |> write_csv(here("outputs", "mat_accuracy.csv"))

# Which perc
gfc_acc <- mat_acc$acc[which.max(mat_acc$acc)]
gfc_perc <- mat_acc$perc[which.max(mat_acc$acc)]

# Plot
p <- mat_acc |>
  ggplot(aes(perc, acc)) +
  geom_point() +
  geom_line(col=grey(0.5)) +
  xlab("Tree cover (%)") +
  ylab("Overall accuracy") +
  scale_x_continuous(limits=c(0, 100), breaks=seq(from=0, to=100, by=20)) +
  scale_y_continuous(limits=c(0.65, 0.85), breaks=seq(from=0.60, to=0.85, by=0.05)) +
  theme_bw(base_size=18) +
  theme(
    panel.border = element_blank(),  # Supprime le cadre complet de la zone de tracé
    axis.line = element_line(color = "black") # Ajoute explicitement les lignes X et Y
  )
ggsave(here("outputs", "plot_accuracy_gfc.pdf"))

# Confusion matrix for best percentage
for2021_pred_gfc <- ifelse((df$treecover2000 >= gfc_perc) & !(df$lossyear %in% c(1:20)),
                           "Forest", "NonForest")
df_cm <- data.frame(
  obs=as.factor(df$for2021_fnf),
  pred=as.factor(for2021_pred_gfc))
df_cm <- na.omit(df_cm)
conf_mat <- caret::confusionMatrix(
  data=df_cm$pred,
  reference=df_cm$obs)
# Save results
sink(here("outputs", "conf_mat_2021_fnf_gfc60.txt"))
conf_mat
sink()

# ========================================================
# TMF fnf 2021
# ========================================================

for2021_pred_tmf <- ifelse(df$Dec2020 %in% c(1, 2, 4), "Forest", "NonForest")
df_cm <- data.frame(
  obs=as.factor(df$for2021_fnf),
  pred=as.factor(for2021_pred_tmf))
df_cm <- na.omit(df_cm)
conf_mat <- caret::confusionMatrix(
  data=df_cm$pred,
  reference=df_cm$obs)
tmf_acc <- as.numeric(conf_mat$overall["Accuracy"])

# Save results
sink(here("outputs", "conf_mat_2021_fnf_tmf.txt"))
conf_mat
sink()

# Comparison with GFC
comp_gfc_tmf_fnf2021 <- data.frame(
  source=c("GFC", "TMF"),
  perc=c(gfc_perc, NA),
  year=rep(2021, 2),
  acc=sapply(c(gfc_acc, tmf_acc), round, digits=2)) |>
  write_csv(here("outputs", "comp_gfc_tmf_fnf2021.csv"))

# ========================================================
# Comparison with AMAP map (~2015)
# ========================================================

# See python script comparing_amap_gfc_tmf_maps.py

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

# Variables
years <- c(2000, 2008, 2021)
nyears <- length(years)
maps <- c("gfc_60", "tmf")
nmaps <- length(maps)

# Data-frame to store results
df_res <- data.frame(
  years=rep(years, nmaps),
  maps=rep(maps, each=nyears),
  ff=NA, fnf=NA, nff=NA, nfnf=NA, n_f=NA, n_nf=NA,
  oa=NA, sen=NA, spe=NA, kappa=NA)

# Loops
for (i in 1:nmaps) {
  for (j in 1:nyears) {  
    cmat <- df_cm |>
      mutate(obs=as.factor(.data[[paste0("for", years[j])]])) |>
      mutate(pred=as.factor(.data[[paste0("for", years[j], "_pred_", maps[i])]])) |>
      select(obs, pred) |>
      na.omit()
    conf_mat <- caret::confusionMatrix(
      data=cmat$pred,
      reference=cmat$obs)
    w <- which(df_res$years==years[j] & df_res$maps==maps[i])
    df_res$oa[w] <- conf_mat$overall["Accuracy"]
    df_res$kappa[w] <- conf_mat$overall["Kappa"]
    df_res$sen[w] <- conf_mat$byClass["Sensitivity"]
    df_res$spe[w] <- conf_mat$byClass["Specificity"]
    df_res$ff[w] <- conf_mat$table[1, 1]
    df_res$fnf[w] <- conf_mat$table[2, 1]
    df_res$nff[w] <- conf_mat$table[1, 2]
    df_res$nfnf[w] <- conf_mat$table[2, 2]
  }
}

# Save
df_res <- df_res |>
  mutate(n_f = ff + fnf, n_nf = nfnf + nff) |>
  mutate(across(c(oa, sen, spe, kappa), ~ round(.x, 3))) |>
  arrange(years) |>
  write_csv(here("outputs", "conf_mat_forest_cover_2000_2008_2021.csv"))

# Much lower OA for forest/non-forest classification (~0.70). Much of the forest is classified as non-forest by the two products.

# =====================================================================
# Confusion matrix for forest cover change in 2000--2008 and 2008--2021
# =====================================================================

# Predicted change from GFC (considering tree cover >= 60%)
df_fcc <- df_cm |>
  # p1
  mutate(fcc_p1_gfc_60=ifelse(for2000_pred_gfc_60=="Forest" &
                                for2008_pred_gfc_60=="Forest", "stableF", NA)) |>
  mutate(fcc_p1_gfc_60=ifelse(for2000_pred_gfc_60=="Forest" &
                                for2008_pred_gfc_60=="NonForest", "loss", fcc_p1_gfc_60)) |>
  mutate(fcc_p1_gfc_60=ifelse(for2000_pred_gfc_60=="NonForest" &
                                for2008_pred_gfc_60=="Forest", "gain", fcc_p1_gfc_60)) |>
  mutate(fcc_p1_gfc_60=ifelse(for2000_pred_gfc_60=="NonForest" &
                                for2008_pred_gfc_60=="NonForest", "stableNF", fcc_p1_gfc_60)) |>
  # p2
  mutate(fcc_p2_gfc_60=ifelse(for2008_pred_gfc_60=="Forest" &
                                for2021_pred_gfc_60=="Forest", "stableF", NA)) |>
  mutate(fcc_p2_gfc_60=ifelse(for2008_pred_gfc_60=="Forest" &
                                for2021_pred_gfc_60=="NonForest", "loss", fcc_p2_gfc_60)) |>
  mutate(fcc_p2_gfc_60=ifelse(for2008_pred_gfc_60=="NonForest" &
                                for2021_pred_gfc_60=="Forest", "gain", fcc_p2_gfc_60)) |>
  mutate(fcc_p2_gfc_60=ifelse(for2008_pred_gfc_60=="NonForest" &
                                for2021_pred_gfc_60=="NonForest", "stableNF", fcc_p2_gfc_60))

# Predicted change from TMF
df_fcc <- df_fcc |>
  # p1
  mutate(fcc_p1_tmf=ifelse(for2000_pred_tmf=="Forest" &
                                for2008_pred_tmf=="Forest", "stableF", NA)) |>
  mutate(fcc_p1_tmf=ifelse(for2000_pred_tmf=="Forest" &
                                for2008_pred_tmf=="NonForest", "loss", fcc_p1_tmf)) |>
  mutate(fcc_p1_tmf=ifelse(for2000_pred_tmf=="NonForest" &
                                for2008_pred_tmf=="Forest", "gain", fcc_p1_tmf)) |>
  mutate(fcc_p1_tmf=ifelse(for2000_pred_tmf=="NonForest" &
                                for2008_pred_tmf=="NonForest", "stableNF", fcc_p1_tmf)) |>
  # p2
  mutate(fcc_p2_tmf=ifelse(for2008_pred_tmf=="Forest" &
                                for2021_pred_tmf=="Forest", "stableF", NA)) |>
  mutate(fcc_p2_tmf=ifelse(for2008_pred_tmf=="Forest" &
                                for2021_pred_tmf=="NonForest", "loss", fcc_p2_tmf)) |>
  mutate(fcc_p2_tmf=ifelse(for2008_pred_tmf=="NonForest" &
                                for2021_pred_tmf=="Forest", "gain", fcc_p2_tmf)) |>
  mutate(fcc_p2_tmf=ifelse(for2008_pred_tmf=="NonForest" &
                                for2021_pred_tmf=="NonForest", "stableNF", fcc_p2_tmf))

# Variables
periods <- c("p1", "p2")
n_periods <- length(periods)
maps <- c("gfc_60", "tmf")
n_maps <- length(maps)
fcc_classes <- c("gain", "loss", "stableF", "stableNF")

# Data-frame to store results
df_freq <- data.frame()
df_acc <- data.frame(
  periods=rep(periods, n_maps),
  maps=rep(maps, each=n_periods),
  sen_gain=NA, spe_gain=NA,
  sen_loss=NA, spe_loss=NA)

# Loops
for (i in 1:n_maps) {
  for (j in 1:n_periods) {
    # Confusion matrix
    cmat <- df_fcc |>
      mutate(obs=factor(.data[[paste0("fcc_", periods[j])]], levels=fcc_classes)) |>
      mutate(pred=factor(.data[[paste0("fcc_", periods[j], "_", maps[i])]],
                         levels=fcc_classes)) |>
      select(obs, pred) |>
      na.omit()
    conf_mat <- caret::confusionMatrix(
      data=cmat$pred,
      reference=cmat$obs)
    # Frequencies
    tab <- data.frame(as.matrix(conf_mat))
    tab$map <- maps[i]
    tab$period <- periods[j]
    df_freq <- df_freq |> bind_rows(tab)
    # Accuracy indices
    w <- which(df_acc$periods==periods[j] & df_acc$maps==maps[i])
    df_acc$sen_gain[w] <- conf_mat$byClass["Class: gain", "Sensitivity"]
    df_acc$spe_gain[w] <- conf_mat$byClass["Class: gain", "Specificity"]
    df_acc$sen_loss[w] <- conf_mat$byClass["Class: loss", "Sensitivity"]
    df_acc$spe_loss[w] <- conf_mat$byClass["Class: loss", "Specificity"]
  }
}

# Save
df_freq |> write.csv(here("outputs", "conf_mat_fcc.csv"), row.names=TRUE)
df_acc <- df_acc |>
  mutate(across(c("sen_gain", "spe_gain", "sen_loss", "spe_loss"), ~ round(.x, 3))) |>
  write_csv(here("outputs", "accurracy_fcc.csv"))


