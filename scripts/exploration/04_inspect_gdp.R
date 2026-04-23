# ============================================================
# exploration/04_inspect_gdp.R
# Inspect the Swiss GDP target variable
# ============================================================
#
# PURPOSE:
#   Understand the GDP series available in the JSON and confirm
#   the correct target variable before modelling.
#
# TARGET VARIABLE CHOICE:
#   We use ch_kof_modelinput_gdpos (GDP excluding sporting
#   events) rather than ch_seco_gdp_real_gdp_ssa (headline
#   GDP) for three reasons:
#   1. FIFA/UEFA revenues create large one-off spikes in
#      headline GDP that have nothing to do with the business
#      cycle — they cannot be predicted by economic indicators.
#   2. Our predictors (CLI, surveys, exchange rates) track
#      the real economic cycle, not sporting revenues.
#   3. KOF themselves use the sports-excluded series in their
#      own nowcasting model — we follow their approach.
#
# AVAILABLE GDP KEYS IN JSON:
#   ch_seco_gdp_real_gdp_ssa — headline real GDP (no _lvl)
#   ch_kof_modelinput_gdpos  — GDP ex-sports (no _lvl)
#   Both exist only as growth rate variants: pct_3m, pct_1y,
#   dif_3m, dif_1y (+ detrended versions). No level series.
#
# ============================================================

library(here)
library(jsonlite)

json_text   <- paste(readLines(here("data", "swiss_nowcast_data.json"),
                               encoding = "UTF-8"), collapse = "")
nowcast_raw <- fromJSON(json_text, flatten = TRUE)


# ============================================================
# 1. QoQ GROWTH RATE — main target variable
# ============================================================

# pct_3m = percentage change over the past 3 months
# = quarter-on-quarter growth rate.
# This is the standard metric reported by KOF and SECO.

dates_qoq  <- as.Date(
  nowcast_raw$dates[["ch_kof_modelinput_gdpos_pct_3m"]],
  format = "%d.%m.%Y")
values_qoq <- as.numeric(
  nowcast_raw[["ch_kof_modelinput_gdpos_pct_3m"]])

cat("--- Swiss GDP ex-sports: QoQ growth (pct_3m) ---\n")
cat("Date range:   ", format(min(dates_qoq), "%Y-%m"),
    "to", format(max(dates_qoq), "%Y-%m"), "\n")
cat("Observations:", length(values_qoq), "quarters\n\n")

cat("Last 12 quarters:\n")
print(data.frame(
  date        = tail(dates_qoq, 12),
  gdp_qoq_pct = round(tail(values_qoq, 12), 3)
))


# ============================================================
# 2. YoY GROWTH RATE — for reference
# ============================================================

# pct_1y = year-on-year growth. Less volatile than QoQ and
# easier to communicate, but lags the cycle more. Reported
# here for reference only — we model in QoQ.

dates_yoy  <- as.Date(
  nowcast_raw$dates[["ch_kof_modelinput_gdpos_pct_1y"]],
  format = "%d.%m.%Y")
values_yoy <- as.numeric(
  nowcast_raw[["ch_kof_modelinput_gdpos_pct_1y"]])

cat("\n--- Swiss GDP ex-sports: YoY growth (pct_1y) ---\n")
cat("Last 8 quarters:\n")
print(data.frame(
  date        = tail(dates_yoy, 8),
  gdp_yoy_pct = round(tail(values_yoy, 8), 3)
))


# ============================================================
# 3. CURRENT NOWCASTING SITUATION
# ============================================================

cat("\n--- Current nowcasting situation ---\n")
cat("Last observed quarter: Q ending",
    format(max(dates_qoq), "%Y-%m"), "\n")
cat("Target quarter:        Q1 2026 (ending 2026-03)\n")
cat("KOF benchmark:         +0.30% QoQ (headline GDP)\n")
cat("\nNote: the KOF benchmark is for headline GDP.\n")
cat("Our model targets GDP ex-sports — the comparison\n")
cat("is approximate but valid since no major sporting\n")
cat("event falls in Q1 2026.\n")
