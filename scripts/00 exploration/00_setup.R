# ============================================================================
# scripts/exploration/00_setup.R
# Setup: packages, global constants, data loading, variable extraction
# Run this script first — all downstream scripts assume these objects exist
# ============================================================================

# ── Packages ─────────────────────────────────────────────────────────────────
library(here)
library(readxl)
library(jsonlite)
library(lubridate)
library(zoo)
library(ggplot2)
library(patchwork)
library(vars)
library(tseries)
library(lmtest)
library(sandwich)
library(tibble)
library(tidyr)
library(dplyr)      # load last so select/filter/mutate are not masked

# ── Plot theme & colours ──────────────────────────────────────────────────────
COL_GDP    <- "#2166ac"   # blue          — KOF sports-corrected
COL_GDP_SA <- "#92c5de"   # light blue    — SECO SA overlay
COL_CPI    <- "#d6604d"   # red-orange    — headline CPI
COL_CPI_SA <- "#f4a582"   # light salmon  — core CPI overlay
COL_BOND   <- "#1b7837"   # dark green    — bond yield level / changes
COL_GREY   <- "gray50"

THEME_BASE <- theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(size = 12, face = "bold", margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 10, color = "gray40", margin = margin(b = 8)),
    axis.title.y     = element_text(size = 9,  color = "gray30"),
    axis.text        = element_text(size = 8,  color = "gray30"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    legend.key.width = unit(1.2, "cm"),
    strip.text       = element_text(size = 10, face = "bold")
  )

theme_set(THEME_BASE)

# ── Paths ─────────────────────────────────────────────────────────────────────
PATH_DATA_JSON <- here("data", "swiss_nowcast_data.json")
PATH_DATA_META <- here("data", "swiss_nowcast_metadata.xlsx")
PATH_FIG       <- here("output", "figures")
PATH_TAB       <- here("output", "tables")
dir.create(PATH_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(PATH_TAB, recursive = TRUE, showWarnings = FALSE)

# ── Load raw data ─────────────────────────────────────────────────────────────
message("Loading metadata...")
metadata <- read_excel(PATH_DATA_META)

message("Loading JSON (this may take a moment)...")
json_text   <- paste(readLines(PATH_DATA_JSON, encoding = "UTF-8"), collapse = "")
nowcast_raw <- fromJSON(json_text, flatten = TRUE)
json_keys   <- names(nowcast_raw)
message("  JSON loaded: ", length(json_keys) - 1, " series")

# ── Helper: extract one series from the JSON ──────────────────────────────────
# key_vals      : name of the numeric vector in nowcast_raw
# key_dates     : name in nowcast_raw$dates (top-level; used for transformed series)
# key_dates_idx : name in nowcast_raw$dates$index (used for raw / level series)
extract_series <- function(key_vals,
                           key_dates     = key_vals,
                           key_dates_idx = NULL) {
  vals <- nowcast_raw[[key_vals]]
  if (is.null(vals)) stop("Key not found in JSON values: ", key_vals)
  
  if (!is.null(key_dates_idx)) {
    dates_chr <- nowcast_raw$dates$index[[key_dates_idx]]
    if (is.null(dates_chr)) stop("Key not found in dates$index: ", key_dates_idx)
  } else {
    dates_chr <- nowcast_raw$dates[[key_dates]]
    if (is.null(dates_chr)) stop("Key not found in dates: ", key_dates)
  }
  
  tibble(date  = dmy(dates_chr),
         value = as.numeric(vals))
}

# ── GDP ───────────────────────────────────────────────────────────────────────
# No level series available in the dataset — percentage growth rates only.
# dif_3m / dif_1y are in absolute CHF (billions) — not useful, excluded.
# Two flavours per transformation:
#   pct_*    = KOF sports-event corrected  (our VAR target)
#   sa_pct_* = SECO seasonally adjusted    (standard official series, overlay)
message("Extracting GDP series...")
gdp <- list(
  pct_3m    = extract_series("ch_kof_modelinput_gdpos_pct_3m"),
  pct_1y    = extract_series("ch_kof_modelinput_gdpos_pct_1y"),
  sa_pct_3m = extract_series("ch_seco_gdp_real_gdp_ssa_pct_3m"),
  sa_pct_1y = extract_series("ch_seco_gdp_real_gdp_ssa_pct_1y")
)

# ── CPI ───────────────────────────────────────────────────────────────────────
# Headline CPI vs Core CPI (excl. food/energy/seasonal ≈ SA proxy)
# Note: raw level values not stored in JSON — growth rates only
message("Extracting CPI series...")
cpi <- list(
  pct_1m      = extract_series("swconprce_pct_1m"),
  dif_1m      = extract_series("swconprce_dif_1m"),
  pct_3m      = extract_series("swconprce_pct_3m"),
  dif_3m      = extract_series("swconprce_dif_3m"),
  pct_1y      = extract_series("swconprce_pct_1y"),
  dif_1y      = extract_series("swconprce_dif_1y"),
  core_pct_1m = extract_series("swcpcoref_pct_1m"),
  core_pct_3m = extract_series("swcpcoref_pct_3m"),
  core_pct_1y = extract_series("swcpcoref_pct_1y"),
  # KOF quarterly series (VAR companion to sports-corrected GDP)
  q_pct_3m    = extract_series("ch_kof_modelinput_cpi_pct_3m"),
  q_pct_1y    = extract_series("ch_kof_modelinput_cpi_pct_1y")
)

# ── Bond yield ────────────────────────────────────────────────────────────────
# 10Y Swiss government bond yield (SNB, end-of-period, monthly)
# Level is already a rate (%) so differences are the stationary transform
message("Extracting bond yield series...")
bond <- list(
  lvl    = extract_series("swgbond__lvl", key_dates_idx = "swgbond_"),
  dif_1m = extract_series("swgbond__dif_1m"),
  dif_3m = extract_series("swgbond__dif_3m"),
  dif_1y = extract_series("swgbond__dif_1y")
)

# ── VAR-ready quarterly dataset ───────────────────────────────────────────────
# Monthly bond → quarterly (sum monthly diffs; average level)
bond_q_dif <- bond$dif_3m %>%
  mutate(date = as.Date(as.yearqtr(date))) %>%
  group_by(date) %>%
  summarise(bond_dif = sum(value, na.rm = TRUE), .groups = "drop")

bond_q_lvl <- bond$lvl %>%
  mutate(date = as.Date(as.yearqtr(date))) %>%
  group_by(date) %>%
  summarise(bond_lvl = mean(value, na.rm = TRUE), .groups = "drop")

var_data <- gdp$pct_3m     %>% rename(gdp_pct_3m = value) %>%
  inner_join(gdp$pct_1y    %>% rename(gdp_pct_1y = value),  by = "date") %>%
  inner_join(cpi$q_pct_3m  %>% rename(cpi_pct_3m = value),  by = "date") %>%
  inner_join(cpi$q_pct_1y  %>% rename(cpi_pct_1y = value),  by = "date") %>%
  inner_join(bond_q_dif,                                     by = "date") %>%
  inner_join(bond_q_lvl,                                     by = "date") %>%
  arrange(date)

# COVID dummies 2020 Q1 – 2021 Q4 (passed as exogen in VAR)
covid_dates  <- as.Date(c(
  "2020-01-01","2020-04-01","2020-07-01","2020-10-01",
  "2021-01-01","2021-04-01","2021-07-01","2021-10-01"))
dummy_labels <- c("d_2020q1","d_2020q2","d_2020q3","d_2020q4",
                  "d_2021q1","d_2021q2","d_2021q3","d_2021q4")
for (i in seq_along(covid_dates))
  var_data[[dummy_labels[i]]] <- as.integer(var_data$date == covid_dates[i])

message("VAR dataset: ", nrow(var_data), " quarters (",
        format(min(var_data$date)), " – ", format(max(var_data$date)), ")")

# ── Series inventory ──────────────────────────────────────────────────────────
message("\n── Series loaded ─────────────────────────────────────────────")
invisible(lapply(c("gdp","cpi","bond"), function(grp_name) {
  grp <- get(grp_name)
  for (nm in names(grp)) {
    s <- grp[[nm]]
    message(sprintf("  %-5s $%-15s %4d obs   %s – %s",
                    grp_name, nm, nrow(s),
                    format(min(s$date)), format(max(s$date))))
  }
}))
message("──────────────────────────────────────────────────────────────\n")