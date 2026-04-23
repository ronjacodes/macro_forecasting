# ============================================================
# exploration/01_inspect.R
# Initial inspection of raw data files
# ============================================================
#
# PURPOSE:
#   First look at the two data files before any processing.
#   Establishes the basic structure of the metadata and JSON,
#   which informed all subsequent data handling decisions.
#
# WHAT WE LEARNED:
#   - Metadata: 366 rows, 4 columns (key, description, source,
#     frequency). 303 monthly, 50 quarterly, 12 daily, 1 NA.
#   - JSON: 4341 top-level fields (4340 series + "dates").
#     Each series is a numeric vector. "dates" is a named list
#     with one character vector per series — not a shared index.
#   - JSON keys are base key + suffix (e.g. "bdeuscciq_lvl").
#     Metadata keys are base only (e.g. "bdeuscciq"). This
#     mismatch is resolved in subsequent scripts.
#
# ============================================================

library(readxl)
library(jsonlite)
library(here)


# ============================================================
# 1. METADATA
# ============================================================

metadata <- read_excel(here("data", "swiss_nowcast_metadata.xlsx"))

cat("--- Metadata ---\n")
cat("Dimensions:", dim(metadata), "\n")
cat("Columns:   ", paste(names(metadata), collapse = ", "), "\n\n")
print(head(metadata, 20))

cat("\nFrequency breakdown:\n")
print(table(metadata$frequency, useNA = "ifany"))


# ============================================================
# 2. JSON — top-level structure
# ============================================================

json_text   <- paste(readLines(here("data", "swiss_nowcast_data.json"),
                               encoding = "UTF-8"), collapse = "")
nowcast_raw <- fromJSON(json_text, flatten = TRUE)

cat("\n\n--- JSON top-level structure ---\n")
cat("Number of top-level fields:", length(nowcast_raw), "\n")
cat("Class:", class(nowcast_raw), "\n")

cat("\nFirst 10 field names:\n")
print(names(nowcast_raw)[1:10])

cat("\nLast 10 field names:\n")
print(tail(names(nowcast_raw), 10))

# Peek at one series to understand the data format
first_key <- names(nowcast_raw)[1]
cat("\n--- Peek at first field:", first_key, "---\n")
cat("Class:", class(nowcast_raw[[first_key]]), "\n")
cat("First 10 values:\n")
print(head(nowcast_raw[[first_key]], 10))
