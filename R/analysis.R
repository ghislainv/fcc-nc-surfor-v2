# Libraries
library(readr)
library(glue)
library(dplyr)
library(here)
library(caret)  # Confusion matrices
library(ggplot2)
library(sf)

# Output directory
dir.create(here("outputs"))

# Load data
df <- read_csv(here("data", "df_surfor_gee.csv"), show_col_types=FALSE)

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

# Much lower OA for forest/non-forest classification (~0.70)A lot of forest is classified as non-forest by the two products

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
status <- c("loss", "gain", "stableF", "stableNF")
n_status <- length(status)
periods <- c("p1", "p2")
n_periods <- length(periods)
maps <- c("gfc_60", "tmf")
n_maps <- length(maps)

# Data-frame to store results
df_res <- data.frame(
  periods=rep(periods, n_maps),
  maps=rep(maps, each=n_periods),
  00=NA, 01=NA, nfnf=NA, nff=NA, n_f=NA, n_nf=NA,
  oa=NA, sen=NA, spe=NA, kappa=NA)

# Loops
for (i in 1:n_maps) {
  for (j in 1:n_periods) {
    for (k in 1:n_status) {
    cmat <- df_fcc |>
      mutate(obs=as.factor(.data[[paste0("fcc_", periods[j])]])) |>
      mutate(pred=as.factor(.data[[paste0("fcc_", periods[j], "_", maps[i])]])) |>
      select(obs, pred) |>
      na.omit()
    conf_mat <- caret::confusionMatrix(
      data=cmat$pred,
      reference=cmat$obs)
    w <- which(df_res$periods==periods[j] & df_res$maps==maps[i])
    df_res$oa[w] <- conf_mat$overall["Accuracy"]
    df_res$kappa[w] <- conf_mat$overall["Kappa"]
    df_res$sen[w] <- conf_mat$byClass["Sensitivity"]
    df_res$spe[w] <- conf_mat$byClass["Specificity"]
    ## df_res$ff[w] <- conf_mat$table[1, 1]
    ## df_res$fnf[w] <- conf_mat$table[2, 1]
    ## df_res$nff[w] <- conf_mat$table[1, 2]
    ## df_res$nfnf[w] <- conf_mat$table[2, 2]
  }
}

# Save
df_res <- df_res |>
  mutate(n_f = ff + fnf, n_nf = nfnf + nff) |>
  mutate(across(c(oa, sen, spe, kappa), ~ round(.x, 3))) |>
  arrange(years) |>
  write_csv(here("outputs", "conf_mat_forest_cover_2000_2008_2021.csv"))


