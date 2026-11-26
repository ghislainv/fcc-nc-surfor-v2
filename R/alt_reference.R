# Effect of interpreter variability on the user's accuracy (UA) of forest
# cover change maps (loss, gain).
#
# Outputs (outputs/9_alt_ref/):
# - alt_ref_p_single.csv: P(other interpreter says forest | label of one interpreter)
# - alt_ref_ceiling.csv: expected UA of loss/gain for a map identical to the
#   consensus, evaluated against simulated reference data (ceiling)
# - alt_ref_loss_gain.csv: UA and composition of pixels mapped as loss/gain
#   (reference classes from the consensus), + expected UA with simulated reference

# Libraries
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(here)

# Output directory
out_dir <- "outputs"
subdir <- "9_alt_ref"
dir.create(here(out_dir, subdir), showWarnings=FALSE, recursive=TRUE)

# Optimal GFC threshold (from analysis.R, section 2)
mat_acc <- read_csv(here(out_dir, "2_gfc_perc", "gfc_perc_mat_accuracy.csv"),
                    show_col_types=FALSE)
gfc_perc <- mat_acc$perc[which.max(mat_acc$acc)]

# Consensus dataset (formatting-data.R), change datasets only
df <- read_csv(here("data", "df_surfor_consensus.csv"), show_col_types=FALSE) |>
  filter(dataset != "for2021")
stopifnot("gain" %in% names(df))

# ===================================
# 1. Individual annotations (raw data)
# ===================================
# Same id shift as in formatting-data.R for the 2021 change dataset.
read_ann <- function(file, date, shift=FALSE) {
  read_csv(here("data_raw", file), show_col_types=FALSE) |>
    transmute(
      date = date,
      plotid = plotid + ifelse(shift, 1000, 0),
      sampleid = sampleid + ifelse(shift, 16000, 0),
      email,
      cls = case_when(
        Classification == "Forêt" ~ "F",
        Classification == "Non Forêt" ~ "NF",
        .default = "NI"))
}
ann <- bind_rows(
  read_ann("38354-samples_2000.csv", "for2004"),
  read_ann("38356-samples_2010.csv", "for2008"),
  read_ann("38355-samples_2021.csv", "for2021", shift=TRUE))

# ==========================================
# 2. Probability of forest for each pixel
# ==========================================
# P(F | label of one interpreter), per date, from pairs of interpreters
# (NI excluded). It includes the error of both interpreters (conservative).
p_single <- ann |>
  inner_join(ann, by=c("date", "plotid", "sampleid"),
             relationship="many-to-many") |>
  filter(email.x != email.y, cls.x != "NI", cls.y != "NI") |>
  summarise(pF_single=mean(cls.y == "F"), .by=c(date, cls.x)) |>
  rename(cls1=cls.x) |>
  arrange(date, cls1)
write_csv(p_single, here(out_dir, subdir, "alt_ref_p_single.csv"))

# Pixels kept in the consensus are selected below (consensus not NA).
# - one F/NF label (single interpreter): pF = P(F | label) from p_single
# - several labels: pF = proportion of F among the F/NF labels
pix <- ann |>
  summarise(nF=sum(cls == "F"), nNF=sum(cls == "NF"),
            .by=c(date, plotid, sampleid)) |>
  mutate(cls1=ifelse(nF == 1, "F", "NF")) |>
  left_join(p_single, by=c("date", "cls1")) |>
  mutate(pF=ifelse(nF + nNF == 1, pF_single, nF / (nF + nNF)))
write_csv(pix, here(out_dir, subdir, "alt_ref_pF_by_pix.csv"))

#' @brief Get the probability of forest for a set of pixels at a given date.
#'
#' @details
#' Looks up, in the global table `pix`, the probability of forest (`pF`)
#' of each pixel of `d` for the date `dt`. Pixels are matched on the
#' (`plotid`, `sampleid`) pair. `pF` is either P(F | label) derived from
#' `p_single` (one F/NF label) or the proportion of F among the F/NF labels
#' (several labels).
#'
#' @param d  data.frame (or tibble) with columns `plotid` and `sampleid`
#'           identifying the pixels.
#' @param dt character. Date/dataset identifier (e.g. "for2004", "for2008",
#'           "for2021"), matched against `pix$date`.
#'
#' @return numeric vector of length `nrow(d)`, in the same order as the rows
#'         of `d`: probability of forest for each pixel. `NA` if a pixel is
#'         not found in `pix` for this date.
#'
#' @note Depends on the global object `pix`.
get_p <- function(d, dt) {
  # Keep annotations for the requested date only
  pt <- filter(pix, date == dt)
  # Match pixels of d to pt on the plotid/sampleid key and return pF
  pt$pF[match(paste(d$plotid, d$sampleid), paste(pt$plotid, pt$sampleid))]
}

# ===========================
# 3. Forest predicted by maps
# ===========================
# GFC forest in year y: tree cover in 2000 >= theta and no loss up to year y-1.
# For 2008 in the period 2004-2008 only (use_gain=TRUE), pixels of GFC gain
# (2000-2012) are also considered as forest, as described in the Methods.
# Gain is not used for 2008 in the period 2008-2021, nor for 2021.
forest <- function(d, map, year, use_gain=FALSE) {
  y <- as.character(year)
  if (map == "gfc") {
    lossmax <- c(`2004`=3, `2008`=7, `2021`=20)[[y]]
    f <- ((d$treecover2000 >= gfc_perc) %in% TRUE) & !(d$lossyear %in% 1:lossmax)
    if (use_gain) f <- f | (d$gain %in% 1)
    f
  } else {
    dec <- c(`2004`="Dec2003", `2008`="Dec2007", `2021`="Dec2020")[[y]]
    d[[dec]] %in% c(1, 2, 4)
  }
}

# ===========================
# 4. Loss and gain per period
# ===========================
# Simulated reference: forest at each date ~ Bernoulli(pF), independently.
# P(loss) = p1 * (1 - p2), P(gain) = (1 - p1) * p2.
# The expected UA of a set of pixels is the mean of this probability
# (no random draw needed).

#' @brief Compute the user's accuracy (UA) of the loss/gain classes for one
#' period, for the GFC and TMF maps and for a "ceiling" map.
#'
#' @details
#' Selects from `df` the pixels of dataset `ds` having a consensus class
#' (`fcc_col` not NA), then retrieves the probability of forest at the start
#' (`y1`) and the end (`y2`) of the period with `get_p()`.
#'
#' The reference is simulated: forest at each date ~ Bernoulli(pF),
#' independently between dates, so that
#' - P(loss) = p1 * (1 - p2)
#' - P(gain) = (1 - p1) * p2
#'
#' The expected UA is the mean of these probabilities over the pixels
#' predicted in a class (no random draw).
#'
#' Two kinds of results are computed:
#' 1. **Ceiling**: UA of a map identical to the consensus (pixels with
#'    consensus = loss or gain), evaluated against the simulated reference.
#' 2. **GFC and TMF maps**: loss = forest at y1 and non-forest at y2;
#'    gain = non-forest at y1 and forest at y2 (via `forest()`). For each
#'    predicted class, computes
#'    - the UA against the consensus (`ua_cons`),
#'    - the expected UA against the simulated reference (`ua_sim`),
#'    - the consensus composition of the predicted pixels (`prop_*`).
#'
#' @param ds        character. Change dataset to analyse (value of
#'                  `df$dataset`).
#' @param y1        numeric. Start year of the period (e.g. 2004, 2008);
#'                  the annotation date used is `paste0("for", y1)`.
#' @param y2        numeric. End year of the period (e.g. 2008, 2021);
#'                  the annotation date used is `paste0("for", y2)`.
#' @param fcc_col   character. Name of the column of `df` holding the
#'                  consensus class ("loss", "gain", "stableF", "stableNF").
#'                  Pixels with NA in this column are excluded.
#' @param period    character. Period label copied into the output tables
#'                  (e.g. "2004--2008").
#' @param gfc_gain  logical. If TRUE, GFC gain (2000-2012) is counted as
#'                  forest at `y2` (see `forest()`). Default FALSE.
#'
#' @return A list of two tibbles:
#' \itemize{
#'   \item `ceiling`: columns `period`, `class` (loss/gain), `n` (number of
#'         consensus pixels in the class), `ua_ceiling` (expected UA).
#'   \item `maps`: one row per map (`map` = "gfc"/"tmf") and predicted class
#'         (`pred_class`), with `period`, `n_pred`, `ua_cons`, `ua_sim`,
#'         `prop_loss`, `prop_gain`, `prop_stableF`, `prop_stableNF`.
#'         Combinations with no predicted pixel (e.g. no GFC gain in
#'         2008--2021) are omitted.
#' }
#'
#' @note Depends on the global objects `df`, `pix` (via `get_p()`) and
#'       `gfc_perc` (via `forest()`). Stops if a pixel has no pF (NA) at
#'       either date.
#'
#' @seealso get_p, forest
run_period <- function(ds, y1, y2, fcc_col, period, gfc_gain=FALSE) {
  d <- df |> filter(dataset == ds, !is.na(.data[[fcc_col]]))
  p1 <- get_p(d, paste0("for", y1))
  p2 <- get_p(d, paste0("for", y2))
  stopifnot(!anyNA(p1), !anyNA(p2))
  cons <- d[[fcc_col]]
  p_ref <- list(loss=p1 * (1 - p2), gain=(1 - p1) * p2)

  # Ceiling: map identical to the consensus
  ceiling <- bind_rows(lapply(c("loss", "gain"), function(cl) {
    tibble(period, class=cl, n=sum(cons == cl),
           ua_ceiling=mean(p_ref[[cl]][cons == cl]))
  }))

  # GFC and TMF maps
  maps <- bind_rows(lapply(c("gfc", "tmf"), function(m) {
    f1 <- forest(d, m, y1)
    f2 <- forest(d, m, y2, use_gain=gfc_gain)
    pred <- list(loss=f1 & !f2, gain=!f1 & f2)
    bind_rows(lapply(c("loss", "gain"), function(cl) {
      if (!any(pred[[cl]])) return(NULL)  # e.g. no GFC gain in 2008-2021
      ref <- cons[pred[[cl]]]
      tibble(period, map=m, pred_class=cl, n_pred=sum(pred[[cl]]),
             ua_cons=mean(ref == cl),
             ua_sim=mean(p_ref[[cl]][pred[[cl]]]),
             prop_loss=mean(ref == "loss"), prop_gain=mean(ref == "gain"),
             prop_stableF=mean(ref == "stableF"),
             prop_stableNF=mean(ref == "stableNF"))
    }))
  }))
  list(ceiling=ceiling, maps=maps)
}

# GFC gain (2000-2012) is used for the period 2004-2008 only
res <- list(
  run_period("change_0408", 2004, 2008, "fcc_p1", "2004--2008", gfc_gain=TRUE),
  run_period("change_0821", 2008, 2021, "fcc_p2", "2008--2021", gfc_gain=FALSE))

bind_rows(lapply(res, `[[`, "ceiling")) |>
  mutate(ua_ceiling=round(ua_ceiling, 3)) |>
  write_csv(here(out_dir, subdir, "alt_ref_ceiling.csv"))

bind_rows(lapply(res, `[[`, "maps")) |>
  mutate(across(starts_with(c("ua_", "prop_")), ~ round(.x, 3))) |>
  write_csv(here(out_dir, subdir, "alt_ref_loss_gain.csv"))

# End
