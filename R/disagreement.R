# Table Of Content
# 1. Disagreement (disag)

# Libraries
library(readr)
library(glue)
library(dplyr)
library(here)
library(caret)  # Confusion matrices
library(ggplot2)
library(sf)
library(tidyr)
library(tibble)

# Output directories
out_dir <- "outputs"
dir.create(here(out_dir), showWarnings=FALSE)

# ========================================================
# 1. Disagreement between interpreters
# ========================================================

subdir <- "1_disag"
dir.create(here(out_dir, subdir), showWarnings=FALSE)

# Formatting data
# ---------------

# Variables to keep
keep_var <- c("plotid", "sampleid", "email", "dataset", "Classification")

# Load datasets for change
file_2004 <- "38354-samples_2000.csv"
file_2008 <- "38356-samples_2010.csv"
file_2021 <- "38355-samples_2021.csv"
df_2004 <- read_csv(here("data_raw", file_2004), show_col_types=FALSE) |>
  mutate(dataset="for2004") |>
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
df_disag <- df_2004 |>
  bind_rows(df_2008) |>
  bind_rows(df_2021) |>
  bind_rows(df_for2021) |>
  mutate(Classification=ifelse(Classification %in% c("Indéterminé", "Non interprétable"), "NotInt", Classification))

# Number of observations per date and interpreter
df_nobs_inter <- df_disag |>
  summarize(n=n(), .by=c(dataset, email)) |>
  pivot_wider(id_cols=email, names_from=dataset, values_from=n) |>
  write_csv(here(out_dir, subdir, "df_disag_nobs_by_interpreter.csv"))
  
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
  ) |>
  rename(S = `stagiaire@oeil.nc`, A = `adrien.bertaud@oeil.nc`, O = `oriane.bruyere@oeil.nc`)

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
get_alpha <- function(data_wide, dataset_names, interp = c("S","A","O"), ni_as_na = FALSE) {
  mat <- data_wide |> filter(dataset %in% dataset_names) |>
    select(all_of(interp)) |> as.matrix()
  if (ni_as_na) mat[mat == 3] <- NA
  m_u <- rowSums(!is.na(mat))
  keep <- m_u >= 2 # keep only paired units
  mat <- mat[keep, , drop = FALSE]; m_u <- m_u[keep]

  # D_o : ordered disagreeing pairs, weighted by 1/(m_u - 1)
  mis <- apply(mat, 1, function(r) { v <- r[!is.na(r)]; sum(outer(v, v, "!=")) })
  D_o <- sum(mis / (m_u - 1)) / sum(m_u)

  # D_e with n(n-1) correction
  vals <- mat[!is.na(mat)]; n <- length(vals); nc <- table(vals)
  D_e <- (n^2 - sum(nc^2)) / (n * (n - 1))

  list(D_o = round(100 * D_o, 1), D_e = round(100 * D_e, 1),
       alpha = round(100 * (1 - D_o / D_e), 1))
}

d3  <- c("for2021_fnf")
ints <- c("S", "A", "O")

# (a) Reference alpha with 3 interpreters, then removing each one, same perimeter
for (ni in c(FALSE, TRUE)) {
  ref <- get_alpha(data_wide, d3, ints, ni)$alpha
  print(sapply(ints, function(i)
    get_alpha(data_wide, d3, setdiff(ints, i), ni)$alpha - ref))
}
##   S    A    O 
## 1.4 -2.6  0.9 
##   S    A    O 
## 1.2 -3.2  1.4

# (b) Bootstrap by plot
dw3 <- data_wide |> filter(dataset %in% d3) |>
  filter(rowSums(!is.na(across(all_of(ints)))) >= 2) |>
  mutate(key = paste(dataset, plotid))
by_plot <- split(dw3, dw3$key)
set.seed(1)
bt <- replicate(1000, {
  b <- bind_rows(by_plot[sample(names(by_plot), replace = TRUE)])
  ref <- get_alpha(b, d3)$alpha
  sapply(ints, function(i) get_alpha(b, d3, setdiff(ints, i))$alpha - ref)
})
apply(bt, 1, quantile, c(.025, .5, .975))     # 95% CI of gain per interpreter
##             S       A       O
## 2.5%  -3.0025 -7.5000 -3.5000
## 50%    1.3000 -2.7000  0.9000
## 97.5%  6.3000  2.3025  5.5025

# (c) Who is in the minority in 2-1 splits?
m <- as.matrix(dw3[, ints])
m <- m[complete.cases(m) & apply(m, 1, function(r) length(unique(r)) == 2), ]
minor <- apply(m, 1, function(r) {
  tb <- table(r); ints[r != as.numeric(names(tb)[which.max(tb)])]
})
table(minor); chisq.test(table(minor))        # H0: 1/3 each
## minor
##   A   O   S 
##  92 115 115 

## 	Chi-squared test for given probabilities

## data:  table(minor)
## X-squared = 3.2857, df = 2, p-value = 0.1934

# -------------------------
library(dplyr); library(tidyr)
ints <- c("S", "A", "O")

# 1. Classe de référence et statut de chaque pixel, par date
ref_one <- function(v) {
  v <- v[!is.na(v)]; n <- length(v); tb <- table(v)
  if (n == 1) {
    if (v == 3) return(c(NA, "excl_NI_single"))
    return(c(v, "ok"))
  }
  if (n == 2) {
    if (v[1] != v[2])  return(c(NA, "excl_disagree_2"))
    if (v[1] == 3)     return(c(NA, "excl_NI_agreed_2"))
    return(c(v[1], "ok"))
  }
  # n == 3
  if (max(tb) == 1)    return(c(NA, "excl_no_majority"))
  maj <- as.numeric(names(tb)[which.max(tb)])
  if (maj == 3)        return(c(NA, "excl_NI_majority"))
  c(maj, "ok")
}

res <- t(apply(data_wide[, ints], 1, ref_one))
ref <- data_wide |>
  select(dataset, plotid, sampleid) |>
  mutate(n_int = rowSums(!is.na(data_wide[, ints])),
         ref_class = as.integer(res[, 1]),      # 1 = forest, 2 = non-forest
         status = res[, 2])

# 2. Décompte par jeu de données et nombre d'interprètes
by_dataset <- ref |>
  count(dataset, n_int, status) |>
  pivot_wider(names_from = status, values_from = n, values_fill = 0) |>
  mutate(n_total = rowSums(across(-c(dataset, n_int))))
by_dataset

# Chiffres du texte (ex. NI à un seul interprète : n sur y)
ref |> filter(n_int == 1) |>
  summarize(n_NI = sum(status == "excl_NI_single"), y = n(), .by = dataset)

# 3. Exclusions pour les changements
# Hypothèse : un même pixel porte le même (plotid, sampleid) à toutes les dates.
count_change <- function(d1, d2) {
  a <- ref |> filter(dataset == d1) |> select(plotid, sampleid, st1 = status)
  b <- ref |> filter(dataset == d2) |> select(plotid, sampleid, st2 = status)
  j <- inner_join(a, b, by = c("plotid", "sampleid"))
  tibble(period = paste(d1, d2, sep = " - "),
         n_pixels = nrow(j),
         excl_date1_only = sum(j$st1 != "ok" & j$st2 == "ok"),
         excl_date2_only = sum(j$st1 == "ok" & j$st2 != "ok"),
         excl_both       = sum(j$st1 != "ok" & j$st2 != "ok"),
         n_excluded      = sum(j$st1 != "ok" | j$st2 != "ok"),   # sans double comptage
         n_kept          = sum(j$st1 == "ok" & j$st2 == "ok"))
}
bind_rows(count_change("for2004", "for2008"),
          count_change("for2008", "for2021"))
