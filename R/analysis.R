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
# Comparison with Birnbaum map (~2015)
# ========================================================

# Data download
url <- "https://zenodo.org/records/12731044/files/amap_carto_3k_20240715.zip?download=1"
destfile <- here("data_raw", "amap_carto_3k_20240715.zip")
download.file(url, destfile, mode="wb")

# Uncompress
out_dir <- here("data_raw", "amap_carto_3k_20240715")
unzip(destfile, exdir=out_dir)

# Get shapefile
shp_file <- list.files(out_dir, pattern="\\.shp$", full.names=TRUE, recursive=TRUE)

# Get CRS
st_crs(st_read(shp_file, quiet=TRUE))

# Shapefile -> GeoPackage
gpkg_file <- file.path(out_dir, "amap_carto_3k_20240715.gpkg")
sf::gdal_utils(
  util="vectortranslate",
  source=shp_file,
  destination=gpkg_file,
  options=c(
    "-f", "GPKG",
    "-nln", "forest_nc",
    "-t_srs", "EPSG:32758"
  )
)

# Rasterize with extent and resolution to define the grid
xmin <- 344000
ymin <- 7488000
xmax <- 765000
ymax <- 7839000
res  <- 30
raster_file <- file.path(out_dir, "forest_birnbaum.tif")
sf::gdal_utils(
  util="rasterize",
  source=gpkg_file,
  destination=raster_file,
  options=c(
    "-tap",
    "-l", "forest_nc",
    "-burn", "1",
    "-a_nodata", "0",
    "-te", as.character(c(xmin, ymin, xmax, ymax)),
    "-tr", as.character(c(res, res)),               
    "-ot", "Byte",
    "-co", "COMPRESS=DEFLATE"
  )
)

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
