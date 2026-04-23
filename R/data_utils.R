# ============================================================
# R/data_utils.R
# Shared data loading and transformation utilities
# ============================================================
#
# USAGE:
#   source(here("R", "data_utils.R"))
#
# REQUIRES:
#   nowcast_raw — parsed JSON object, loaded by the calling
#   script before sourcing this file:
#
#     json_text   <- paste(readLines(here("data",
#                     "swiss_nowcast_data.json"),
#                     encoding = "UTF-8"), collapse = "")
#     nowcast_raw <- fromJSON(json_text, flatten = TRUE)
#
# ============================================================


# ------------------------------------------------------------
# pull_series()
# Pull a single series from the nowcast JSON by key name.
#
# Arguments:
#   key  — full JSON key including suffix, e.g. "swobs085q_lvl"
#
# Returns:
#   data.frame with columns: date (Date), value (numeric)
#
# Note: dates in the JSON are strings in "DD.MM.YYYY" format,
# one date vector per series (not a shared index).
# ------------------------------------------------------------

pull_series <- function(key) {
  if (!exists("nowcast_raw", envir = parent.env(environment()))) {
    stop("nowcast_raw not found. Load the JSON before calling pull_series().")
  }
  if (!key %in% names(nowcast_raw)) {
    stop(paste("Key not found in JSON:", key))
  }
  dates  <- as.Date(nowcast_raw$dates[[key]], format = "%d.%m.%Y")
  values <- as.numeric(nowcast_raw[[key]])
  data.frame(date = dates, value = values, stringsAsFactors = FALSE)
}


# ------------------------------------------------------------
# to_quarterly()
# Aggregate a monthly series to quarterly frequency by
# averaging the three monthly observations within each quarter.
#
# Only complete quarters (all 3 months present) are returned.
# Incomplete quarters at the ragged edge of the sample are
# dropped here and handled separately in the nowcast step.
#
# Quarter convention: labelled by the first month of the
# quarter — Q1 = January (01), Q2 = April (04), etc.
# This matches the date convention used in the GDP series.
#
# Arguments:
#   df       — data.frame with columns date (Date), value (numeric)
#   varname  — character, name to give the aggregated column
#
# Returns:
#   data.frame with columns: qdate (Date), <varname> (numeric)
# ------------------------------------------------------------

to_quarterly <- function(df, varname) {
  df %>%
    dplyr::mutate(
      year    = as.integer(format(date, "%Y")),
      month   = as.integer(format(date, "%m")),
      quarter = ceiling(month / 3),
      qdate   = as.Date(paste(year,
                              sprintf("%02d", (quarter - 1) * 3 + 1),
                              "01", sep = "-"))
    ) %>%
    dplyr::group_by(qdate) %>%
    dplyr::summarise(
      n_months        = dplyr::n(),
      !!varname       := mean(value, na.rm = TRUE),
      .groups         = "drop"
    ) %>%
    dplyr::filter(n_months == 3) %>%
    dplyr::select(qdate, dplyr::all_of(varname))
}


# ------------------------------------------------------------
# get_partial_avg()
# Compute predictor averages using only the first n_months
# monthly observations of a given quarter. Used in the bridge
# model to condition on partial within-quarter information.
#
# Arguments:
#   monthly_panel — long data.frame with columns:
#                   date, value, variable, qdate, month_in_q
#                   (built once in each analysis script)
#   pred_names    — character vector of predictor names
#   tgt_date      — Date, the quarter to compute averages for
#   n_months      — integer 1, 2, or 3
#
# Returns:
#   1-row data.frame with one column per predictor, or
#   empty data.frame if data not available
# ------------------------------------------------------------

get_partial_avg <- function(monthly_panel, pred_names,
                            tgt_date, n_months) {
  monthly_panel %>%
    dplyr::filter(qdate == tgt_date, month_in_q <= n_months) %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      avg = mean(value, na.rm = TRUE),
      n   = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(n == n_months) %>%
    tidyr::pivot_wider(names_from  = variable,
                       values_from = avg,
                       id_cols     = NULL) %>%
    dplyr::select(dplyr::any_of(pred_names))
}


# ------------------------------------------------------------
# build_monthly_panel()
# Build the long-format monthly panel used by get_partial_avg().
# Call this once per script after loading the predictors list.
#
# Arguments:
#   predictors — named list of data.frames (date, value),
#                one per predictor
#
# Returns:
#   long data.frame with columns:
#   date, value, variable, year, month, quarter,
#   month_in_q, qdate
# ------------------------------------------------------------

build_monthly_panel <- function(predictors) {
  dplyr::bind_rows(lapply(names(predictors), function(nm) {
    predictors[[nm]] %>%
      dplyr::mutate(
        variable   = nm,
        year       = as.integer(format(date, "%Y")),
        month      = as.integer(format(date, "%m")),
        quarter    = ceiling(month / 3),
        month_in_q = month - (quarter - 1) * 3,
        qdate      = as.Date(paste(year,
                                   sprintf("%02d",
                                           (quarter - 1) * 3 + 1),
                                   "01", sep = "-"))
      )
  }))
}


# ------------------------------------------------------------
# build_quarterly_panel()
# Build the full quarterly panel: GDP + all predictors merged,
# with COVID dummy added.
#
# Arguments:
#   gdp        — data.frame from pull_series(), renamed gdp
#   predictors — named list of monthly data.frames
#
# Returns:
#   data.frame with columns:
#   qdate, gdp, <pred1>, ..., <predN>, covid
# ------------------------------------------------------------

build_quarterly_panel <- function(gdp, predictors) {
  quarterly <- gdp %>% dplyr::rename(qdate = date)
  for (nm in names(predictors)) {
    quarterly <- dplyr::left_join(
      quarterly,
      to_quarterly(predictors[[nm]], nm),
      by = "qdate"
    )
  }
  quarterly %>%
    dplyr::mutate(
      covid = as.integer(
        qdate %in% as.Date(c("2020-04-01", "2020-07-01")))
    )
}
