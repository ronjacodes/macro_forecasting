# ============================================================
# exploration/02_inspect_json.R
# Deep inspection of JSON naming convention and date structure
# ============================================================
#
# PURPOSE:
#   Understand how series keys, suffixes, and dates are
#   structured in the JSON before building the data pipeline.
#
# WHAT WE LEARNED:
#   - "dates" is a named list — each series has its own date
#     vector in "DD.MM.YYYY" format. Access via:
#     nowcast_raw$dates[[key]]
#   - Each base series has 14 transformation variants:
#     _lvl, _lvl_detrended,
#     _pct_1m/3m/1y, _pct_1m/3m/1y_detrended,
#     _dif_1m/3m/1y, _dif_1m/3m/1y_detrended
#   - Series lengths vary — not all start/end at same date.
#   - Some series (ch_kof_modelinput_*, ch_seco_gdp_*) are
#     quarterly, stepping by 3 months instead of 1.
#
# ============================================================

library(readxl)
library(jsonlite)
library(here)

metadata    <- read_excel(here("data", "swiss_nowcast_metadata.xlsx"))
json_text   <- paste(readLines(here("data", "swiss_nowcast_data.json"),
                               encoding = "UTF-8"), collapse = "")
nowcast_raw <- fromJSON(json_text, flatten = TRUE)
all_keys    <- setdiff(names(nowcast_raw), "dates")


# ============================================================
# 1. DATES FIELD
# ============================================================

# "dates" is not a single shared vector — it is a named list
# with one character vector per series. Each series has its
# own date index of potentially different length and range.

cat("--- dates field ---\n")
cat("Class:", class(nowcast_raw$dates), "\n")
cat("Length (number of series with dates):",
    length(nowcast_raw$dates), "\n")
cat("\nFirst 5 date vectors (name + first date):\n")
for (i in 1:5) {
  nm  <- names(nowcast_raw$dates)[i]
  d1  <- nowcast_raw$dates[[nm]][1]
  cat(sprintf("  %-40s  %s\n", nm, d1))
}


# ============================================================
# 2. NAMING CONVENTION — suffixes
# ============================================================

# For one known base key, find all its variants in the JSON
example_base <- "bdeuscciq"
variants <- all_keys[startsWith(all_keys, example_base)]
cat("\n\n--- Variants of", example_base, "in JSON ---\n")
print(variants)

# Extract all unique suffixes by stripping known base key prefixes
suffixes <- unique(unlist(lapply(metadata$key, function(b) {
  v <- all_keys[startsWith(all_keys, paste0(b, "_"))]
  sub(paste0("^", b, "_"), "", v)
})))
cat("\nAll unique suffixes found across dataset:\n")
print(sort(suffixes))


# ============================================================
# 3. ONE SERIES END-TO-END
# ============================================================

# Verify that values and dates align correctly for a known key
key <- "bdeuscciq_lvl"
dates  <- as.Date(nowcast_raw$dates[[key]], format = "%d.%m.%Y")
values <- as.numeric(nowcast_raw[[key]])

cat("\n\n--- End-to-end check:", key, "---\n")
cat("Length of values:", length(values), "\n")
cat("Length of dates: ", length(dates),  "\n")
cat("Date range:", format(min(dates), "%Y-%m"),
    "to", format(max(dates), "%Y-%m"), "\n")
cat("First 5 observations:\n")
print(data.frame(date = head(dates, 5), value = head(values, 5)))
