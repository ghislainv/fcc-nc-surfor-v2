# Table Of Content
# 1. Disagreement (disag)

# Libraries
library(readr)
library(dplyr)
library(here)
library(tidyr)
library(tibble)

# Output directories
out_dir <- "outputs"
dir.create(here(out_dir), showWarnings=FALSE)
subdir <- "1_disag"
dir.create(here(out_dir, subdir), showWarnings=FALSE)

# ---------------
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
  write_csv(here(out_dir, subdir, "out_disag_nobs_by_interpreter.csv"))
  
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

# ---------------------
# Excluded observations
# ---------------------
ints <- c("S", "A", "O")

# Reference class and pixel status per date
ref_one <- function(v) {
  v <- v[!is.na(v)]; n <- length(v); tb <- table(v)
  if (n == 1) {
    if (v == 3) return(c(NA, "excl_NI_single"))
    return(c(v, "ok"))
  }
  if (n == 2) {
    if (v[1] != v[2])  return(c(NA, "excl_disag_2"))
    if (v[1] == 3)     return(c(NA, "excl_NI_agreed_2"))
    return(c(v[1], "ok"))
  }
  # n == 3
  if (max(tb) == 1)    return(c(NA, "excl_no_maj"))
  maj <- as.numeric(names(tb)[which.max(tb)])
  if (maj == 3)        return(c(NA, "excl_NI_maj"))
  c(maj, "ok")
}

res <- t(apply(data_wide[, ints], 1, ref_one))
ref <- data_wide |>
  select(dataset, plotid, sampleid) |>
  mutate(n_int = rowSums(!is.na(data_wide[, ints])),
         ref_class = as.integer(res[, 1]),      # 1 = forest, 2 = non‑forest
         status = res[, 2])

# Count by dataset and number of interpreters
by_dataset <- ref |>
  count(dataset, n_int, status) |>
  pivot_wider(names_from = status, values_from = n, values_fill = 0)
if (!("excl_no_maj" %in% names(by_dataset))) by_dataset["excl_no_maj"] <- 0
by_dataset <- by_dataset |>
  mutate(n_total = rowSums(across(-c(dataset, n_int))))
write_csv(by_dataset, here(out_dir, subdir, "out_disag_excluded_obs_dataset.csv"))

# Valid obs per dataset
valid_obs <- by_dataset |>
  group_by(dataset) |>
  summarise(n_ok = sum(ok), n_init = sum(n_total)) |>
  ungroup() |>
  mutate(perc_valid = round(100 * n_ok / n_init, 1))
write_csv(valid_obs, here(out_dir, subdir, "out_disag_valid_obs_dataset.csv"))

# Numbers for the text (e.g., NI with a single interpreter: n out of y)
out_single_int <- ref |> 
  filter(n_int == 1) |>
  summarize(n_NI = sum(status == "excl_NI_single"), y = n(), .by = dataset)
write_csv(out_single_int, here(out_dir, subdir, "out_disag_out_single_int.csv"))

# Exclusions for changes
# The same pixel has the same (plotid, sampleid) across dates.
count_change <- function(d1, d2) {
  a <- ref |> filter(dataset == d1) |> select(plotid, sampleid, st1 = status)
  b <- ref |> filter(dataset == d2) |> select(plotid, sampleid, st2 = status)
  j <- inner_join(a, b, by = c("plotid", "sampleid"))
  tibble(period = paste(d1, d2, sep = " - "),
         n_pixels = nrow(j),
         excl_date1_only = sum(j$st1 != "ok" & j$st2 == "ok"),
         excl_date2_only = sum(j$st1 == "ok" & j$st2 != "ok"),
         excl_both       = sum(j$st1 != "ok" & j$st2 != "ok"),
         n_excluded      = sum(j$st1 != "ok" | j$st2 != "ok"),   # without double counting
         n_kept          = sum(j$st1 == "ok" & j$st2 == "ok"))
}
res <- bind_rows(count_change("for2004", "for2008"),
                 count_change("for2008", "for2021")) |>
  mutate(perc_kept = round(100 * n_kept / n_pixels, 1))
write_csv(res, here(out_dir, subdir, "out_disag_excluded_obs_change.csv"))

# ----------------------------
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

# -----------------------------
# Krippendorff's Alpha function
# -----------------------------
# NI (code 3) is treated as a class in its own right.
get_alpha <- function(data_wide, dataset_names, interp = c("S", "A", "O")) {
  mat <- data_wide |>
    filter(dataset %in% dataset_names) |>
    select(all_of(interp)) |>
    as.matrix()
  m_u <- rowSums(!is.na(mat))
  keep <- m_u >= 2 # keep only paired units
  mat <- mat[keep, , drop = FALSE]; m_u <- m_u[keep]

  # D_o : ordered disagreeing pairs, weighted by 1/(m_u - 1)
  mis <- apply(mat, 1, function(r) { v <- r[!is.na(r)]; sum(outer(v, v, "!=")) })
  D_o <- sum(mis / (m_u - 1)) / sum(m_u)

  # D_e with n(n-1) correction
  vals <- mat[!is.na(mat)]; n <- length(vals); nc <- table(vals)
  D_e <- (n^2 - sum(nc^2)) / (n * (n - 1))

  list(D_o = round(100 * D_o, 2),
       D_e = round(100 * D_e, 2),
       alpha = round(100 * (1 - D_o / D_e), 2))
}

# Loop on datasets
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
  write_csv(here(out_dir, subdir, "out_disag_by_dataset.csv"))

#--------------------------
# Leave-one-out interpreter
# -------------------------
d3  <- c("for2008", "for2021_fnf")
ints <- c("S", "A", "O")

# Reference alpha with 3 interpreters, then removing each one, same perimeter
ref_alpha <- get_alpha(data_wide, d3, ints)$alpha
loo <- sapply(ints, function(i)
  get_alpha(data_wide, d3, setdiff(ints, i))$alpha - ref_alpha)
##   S    A    O 
## 1.4 -2.6  0.9 

# Bootstrap by plot
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
loo_bt <- apply(bt, 1, quantile, c(.025, .5, .975))     # 95% CI of gain per interpreter
##              S        A     O
## 2.5%  -0.03000 -0.07000 -0.03
## 50%    0.01000 -0.03000  0.01
## 97.5%  0.06025  0.02025  0.06

# Combining leave one out results
loo_tbl <- tibble(
  interpreter = names(loo),
  gain        = unname(loo),
  ci_low      = loo_bt["2.5%",  names(loo)],
  median_bt   = loo_bt["50%",   names(loo)],
  ci_high     = loo_bt["97.5%", names(loo)]
)
write_csv(loo_tbl, here(out_dir, subdir, "out_disag_loo_interpreter_alpha_gain.csv"))

# -------------------------------------
# Who is in the minority in 2-1 splits?
# -------------------------------------
m <- as.matrix(dw3[, ints])
m <- m[complete.cases(m) & apply(m, 1, function(r) length(unique(r)) == 2), ]
minor <- apply(m, 1, function(r) {
  tb <- table(r); ints[r != as.numeric(names(tb)[which.max(tb)])]
})
minor_tbl <- as.data.frame(table(minor), stringsAsFactors = FALSE) |>
  dplyr::rename(interpreter = minor, n = Freq) |>
  dplyr::mutate(prop = n / sum(n),
                expected = sum(n) / 3)
ct <- chisq.test(table(minor))
minor_tbl <- minor_tbl |>
  dplyr::mutate(chisq = unname(ct$statistic),
                df = unname(ct$parameter),
                p_value = ct$p.value)
write_csv(minor_tbl, here(out_dir, subdir, "out_disag_minority_2_1_splits.csv"))

# End
