# ============================================================
# DATA CLEANING & EXPLORATION
# Hierarchical VAR-X Nowcasting of Swiss Overnight Stays
# ============================================================
#
# RESEARCH QUESTION:
#   Forecast Swiss overnight stays hierarchically across cantons
#   and visitor origins using VAR-X (VAR with exogenous nowcast
#   predictors).
#
# WHY VAR-X AND NOT PLAIN VAR:
#   With 252 monthly observations, a pure VAR with many variables
#   is not feasible because parameter count explodes with k^2 * p.
#   VAR-X treats nowcast indicators as exogenous (not modelled
#   jointly), reducing parameters to k*m. This lets us use
#   ~20-30 meaningful predictors without overfitting.
#
# DATA SOURCES:
#   1. ons_data_hierarchy_FS26.xlsx  — Swiss overnight stays
#   2. swiss_nowcast_data.json       — Nowcast predictor series
#   3. swiss_nowcast_metadata.xlsx   — Variable descriptions
#   4. kof_forecasts.xlsx            — KOF macro forecasts
#
# HIERARCHY:
#   Level 0: total.total             (1 series)
#   Level 1: [canton].total         (13 series)
#   Level 2: total.[origin]         (10 series)
#   Level 3: [canton].[origin]     (130 series)
# ============================================================


# ============================================================
# 0. SETUP
# ============================================================

packages <- c("readxl", "jsonlite", "here", "dplyr", "tidyr",
              "ggplot2", "zoo", "lubridate", "scales", "corrplot",
              "tseries", "gridExtra")

# Load packages (Installation is handled manually via renv::install)
lapply(packages, library, character.only = TRUE)

cantons <- c("ag","be","bs","fr","ge","gr","ju","lu","os","ti","vd","vs","zh")
destinations <- c("ch","europa","de","fr","it","amerika","asien","china","golf","rest")


# ============================================================
# PART 1: ONS DATA — SWISS OVERNIGHT STAYS
# ============================================================

cat("\n============================================================\n")
cat("PART 1: ONS DATA\n")
cat("============================================================\n")

# --- 1.1 Load & parse ---------------------------------------
ons_raw      <- read_excel(here("data", "ons_data_hierarchy_FS26.xlsx"),
                           sheet = "hierarchy")
ons_raw$date <- as.Date(paste0(ons_raw$date, "-01"))
ons_raw      <- ons_raw[order(ons_raw$date), ]

cat("Dimensions:   ", dim(ons_raw), "\n")
cat("Date range:    ", format(min(ons_raw$date), "%Y-%m"),
    "to", format(max(ons_raw$date), "%Y-%m"), "\n")
cat("Observations: ", nrow(ons_raw), "months\n")

# --- 1.2 Hierarchy structure --------------------------------
canton_total_cols <- paste0(CANTONS, ".total")
origin_cols       <- paste0("total.", DESTINATIONS)

cat("\n--- Hierarchy ---\n")
cat("Level 0: total.total                  (1 series)\n")
cat("Level 1: [canton].total              (13 series)\n")
cat("Level 2: total.[origin]              (10 series)\n")
cat("Level 3: [canton].[origin]          (130 series)\n")
cat("Total columns:", ncol(ons_raw) - 1, "\n")

# --- 1.3 Missing values -------------------------------------
na_count <- sum(is.na(ons_raw[, -1]))
cat("\nMissing values:", na_count,
    ifelse(na_count == 0, "— data is complete ✓\n", "\n"))

# --- 1.4 Hierarchy consistency ------------------------------
cat("\n--- Hierarchy Consistency ---\n")

canton_sum  <- rowSums(ons_raw[, canton_total_cols])
diff_canton <- abs(ons_raw$total.total - canton_sum)
cat("Canton totals sum to national:",
    ifelse(max(diff_canton) == 0, "YES ✓\n", paste("NO, max diff:", max(diff_canton), "\n")))

origin_sum  <- rowSums(ons_raw[, origin_cols])
diff_origin <- abs(ons_raw$total.total - origin_sum)
cat("Origin totals sum to national:",
    ifelse(max(diff_origin) == 0, "YES ✓\n",
           paste("NO, max diff:", round(max(diff_origin)),
                 "— expected: 'de','fr','it' are subsets of 'europa'\n")))

# --- 1.5 Summary statistics ---------------------------------
cat("\n--- Summary: National Total ---\n")
print(summary(ons_raw$total.total))

cat("\nAverage monthly stays by canton (thousands):\n")
print(round(sort(colMeans(ons_raw[, canton_total_cols]),
                 decreasing = TRUE) / 1000, 1))

cat("\nAverage monthly stays by visitor origin (thousands):\n")
print(round(sort(colMeans(ons_raw[, origin_cols]),
                 decreasing = TRUE) / 1000, 1))

# --- 1.6 Annual totals + COVID impact -----------------------
cat("\n--- Annual Totals (millions) ---\n")
annual <- ons_raw %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(total_mio  = round(sum(total.total) / 1e6, 2),
            pct_change = round((sum(total.total) /
                                  lag(sum(total.total)) - 1) * 100, 1),
            .groups = "drop")
print(annual, n = 25)

# --- 1.7 Plot: National total over time ---------------------
p1 <- ggplot(ons_raw, aes(x = date, y = total.total / 1e6)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  annotate("rect",
           xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-01"),
           ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1) +
  annotate("text", x = as.Date("2021-01-01"), y = 1.2,
           label = "COVID-19", color = "darkred", size = 3.5, fontface = "bold") +
  scale_y_continuous(labels = comma) +
  labs(title    = "Swiss Total Overnight Stays",
       subtitle = "Monthly, January 2005 – January 2026",
       x = NULL, y = "Overnight Stays (millions)") +
  theme_minimal()
print(p1)

# --- 1.8 Plot: By visitor origin ----------------------------
origin_long <- ons_raw %>%
  select(date, all_of(origin_cols)) %>%
  pivot_longer(-date, names_to = "origin", values_to = "stays") %>%
  mutate(origin = gsub("total\\.", "", origin))

p2 <- ggplot(origin_long, aes(x = date, y = stays / 1000, color = origin)) +
  geom_line(linewidth = 0.6) +
  scale_y_continuous(labels = comma) +
  labs(title    = "Overnight Stays by Visitor Origin",
       subtitle = "Monthly, 2005–2026",
       x = NULL, y = "Stays (thousands)", color = "Origin") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p2)

# --- 1.9 Plot: Top 5 cantons --------------------------------
canton_means  <- sort(colMeans(ons_raw[, canton_total_cols]), decreasing = TRUE)
top5_cantons  <- names(canton_means)[1:5]

canton_long <- ons_raw %>%
  select(date, all_of(top5_cantons)) %>%
  pivot_longer(-date, names_to = "canton", values_to = "stays") %>%
  mutate(canton = gsub("\\.total", "", toupper(canton)))

p3 <- ggplot(canton_long, aes(x = date, y = stays / 1000, color = canton)) +
  geom_line(linewidth = 0.7) +
  scale_y_continuous(labels = comma) +
  labs(title    = "Top 5 Cantons by Overnight Stays",
       subtitle = "Monthly, 2005–2026",
       x = NULL, y = "Stays (thousands)", color = "Canton") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p3)

# --- 1.10 Seasonality (COVID years excluded) ----------------
ons_raw$month <- month(ons_raw$date, label = TRUE)

p4 <- ggplot(ons_raw %>% filter(!(year(date) %in% c(2020, 2021))),
             aes(x = month, y = total.total / 1e6)) +
  geom_boxplot(fill = "steelblue", alpha = 0.5, outlier.color = "red") +
  scale_y_continuous(labels = comma) +
  labs(title    = "Seasonality: Swiss Total Overnight Stays",
       subtitle = "COVID years 2020–2021 excluded",
       x = NULL, y = "Stays (millions)") +
  theme_minimal()
print(p4)

# --- 1.11 Stationarity check on key series ------------------
cat("\n--- ADF Stationarity Tests (national + top 5 cantons) ---\n")
test_cols <- c("total.total", top5_cantons)
for (col in test_cols) {
  pval <- adf.test(na.omit(ons_raw[[col]]))$p.value
  cat(sprintf("%-20s p = %.4f  %s\n", col, pval,
              ifelse(pval < 0.05, "stationary ✓", "NON-STATIONARY — needs differencing")))
}

# --- 1.12 Clean ONS -----------------------------------------
ons_clean <- ons_raw %>%
  select(-month) %>%
  mutate(covid_dummy   = as.integer(date >= as.Date("2020-03-01") &
                                      date <= as.Date("2021-12-01")),
         month_dummy   = as.integer(format(date, "%m")))  # 1–12 for seasonal dummies

cat("\n--- ONS Clean: ready ---\n")
cat("Dimensions:", dim(ons_clean), "\n")
cat("Added: covid_dummy, month_dummy\n")


# ============================================================
# PART 2: NOWCAST METADATA
# ============================================================

cat("\n============================================================\n")
cat("PART 2: NOWCAST METADATA\n")
cat("============================================================\n")

metadata <- read_excel(here("data", "swiss_nowcast_metadata.xlsx"))

cat("Dimensions:", dim(metadata), "\n")

cat("\n--- Variables by Frequency ---\n")
print(table(metadata$frequency, useNA = "ifany"))

cat("\n--- Top 10 Sources ---\n")
print(sort(table(metadata$source), decreasing = TRUE)[1:10])

# Country breakdown
metadata$country <- sapply(strsplit(metadata$description, ","), `[`, 1)
cat("\n--- Variables by Country/Region (monthly only) ---\n")
monthly_meta <- metadata %>% filter(frequency == "Month")
print(sort(table(monthly_meta$country), decreasing = TRUE))

cat("\nTotal monthly variables:", nrow(monthly_meta), "\n")
cat("Non-monthly (need interpolation or exclusion):",
    nrow(metadata) - nrow(monthly_meta), "\n")

View(metadata)


# ============================================================
# PART 3: NOWCAST PREDICTOR SERIES (JSON)
# ============================================================

cat("\n============================================================\n")
cat("PART 3: NOWCAST PREDICTOR SERIES\n")
cat("============================================================\n")

# --- 3.1 Load -----------------------------------------------
json_text   <- paste(readLines(here("data", "swiss_nowcast_data.json"),
                               encoding = "UTF-8"), collapse = "")
nowcast_raw <- fromJSON(json_text, flatten = TRUE)

series_keys <- setdiff(names(nowcast_raw), "dates")
cat("Total series in JSON:", length(series_keys), "\n")
cat("Matched to metadata: ", sum(series_keys %in% metadata$key),
    "/", length(series_keys), "\n")

# --- 3.2 Monthly keys available -----------------------------
monthly_keys <- metadata$key[metadata$frequency == "Month" &
                               metadata$key %in% series_keys]
cat("Monthly series for modelling:", length(monthly_keys), "\n")

# --- 3.3 Build tidy long data frame -------------------------
cat("Building tidy data frame of all monthly series...\n")
nowcast_list <- lapply(monthly_keys, function(k) {
  vals  <- nowcast_raw[[k]]
  dates <- as.Date(nowcast_raw$dates[[k]])
  data.frame(date = dates, value = as.numeric(vals),
             key = k, stringsAsFactors = FALSE)
})
nowcast_long <- bind_rows(nowcast_list)

# --- 3.4 Coverage analysis ----------------------------------
cat("\n--- Series Coverage ---\n")
ons_start <- as.Date("2005-01-01")
ons_end   <- as.Date("2026-01-01")

coverage <- nowcast_long %>%
  group_by(key) %>%
  summarise(
    start    = min(date),
    end      = max(date),
    n_obs    = n(),
    n_na     = sum(is.na(value)),
    pct_na   = round(mean(is.na(value)) * 100, 1),
    .groups  = "drop"
  ) %>%
  left_join(metadata[, c("key","description","source","country")],
            by = "key") %>%
  mutate(
    covers_ons    = start <= ons_start & end >= ons_end,
    description_s = substr(description, 1, 55)
  )

cat("Series covering full ONS period (2005–2026):",
    sum(coverage$covers_ons, na.rm = TRUE), "/", nrow(coverage), "\n")
cat("Series with any NAs:", sum(coverage$n_na > 0), "\n")
cat("Series with >10% NAs:", sum(coverage$pct_na > 10), "\n")

# --- 3.5 PREDICTOR SELECTION RATIONALE ---------------------
# VAR-X parameter budget with 252 observations:
#   k=1 endogenous, p=2 lags → 3 endogenous params per equation
#   Each exogenous variable adds 1 parameter per equation
#   Safe budget: obs/params >= 5 → max ~47 exogenous vars
#   Conservative choice: ~25 exogenous predictors
#   Selected to cover: Swiss leading indicators, key visitor
#   origin countries (DE, FR, IT, CH domestic, US, CN), and
#   financial/monetary conditions

cat("\n--- PREDICTOR SELECTION ---\n")
cat("VAR-X budget: ~25 exogenous predictors (obs/param >= 10)\n\n")

# Full curated list of ~25 predictors with justification
predictor_selection <- data.frame(
  key = c(
    # --- Swiss leading indicators (7) ---
    "swobs085q",    # CH OECD composite → overall Swiss cycle
    "swcnfbusq",    # KOF manufacturing survey → business cycle
    "swpmidltq",    # PMI delivery times → supply side
    "swun_p_totq",  # Unemployment rate → labour market
    "swxsfec_",     # CHF/EUR exchange rate → price competitiveness
    "swxrusd_",     # CHF/USD exchange rate → US visitors
    "swrettotg",    # Retail trade turnover → domestic demand
    
    # --- Germany (3) — top visitor origin ---
    "bdeuscciq",    # DE consumer confidence
    "bdeusicir",    # DE industrial confidence
    "bdiptot_g",    # DE industrial production
    
    # --- France (2) ---
    "frcnfconq",    # FR consumer confidence
    "friptot_g",    # FR industrial production
    
    # --- Italy (2) ---
    "iteuscciq",    # IT consumer confidence
    "itiptot_g",    # IT industrial production
    
    # --- Euro Area (3) ---
    "ekeusesig",    # EA economic sentiment
    "emm3____b",    # EA M3 money supply → financial conditions
    "empmia_hq",    # EA composite PMI
    
    # --- United States (3) ---
    "uscnfconq",    # US consumer confidence
    "usipman_g",    # US industrial production
    "usobs085q",    # US OECD composite
    
    # --- China (2) ---
    "chcnfconr",    # CN consumer confidence
    "chpmim_hq",    # CN manufacturing PMI
    
    # --- Global (1) ---
    "wdcampimf"     # World commodity prices → energy/cost proxy
  ),
  group = c(
    rep("Switzerland", 7),
    rep("Germany",     3),
    rep("France",      2),
    rep("Italy",       2),
    rep("Euro Area",   3),
    rep("USA",         3),
    rep("China",       2),
    rep("Global",      1)
  ),
  stringsAsFactors = FALSE
)

# Check availability
predictor_selection <- predictor_selection %>%
  left_join(metadata[, c("key","description","frequency")], by = "key") %>%
  mutate(
    in_json      = key %in% series_keys,
    description_s = substr(description, 1, 55)
  )

cat("Predictors found in JSON:", sum(predictor_selection$in_json),
    "/", nrow(predictor_selection), "\n\n")

# Print grouped
for (grp in unique(predictor_selection$group)) {
  cat(sprintf("  [%s]\n", grp))
  sub <- predictor_selection %>% filter(group == grp)
  for (i in 1:nrow(sub)) {
    status <- ifelse(sub$in_json[i], "✓", "✗ NOT FOUND")
    cat(sprintf("    %-28s %s  %s\n",
                sub$key[i], status, sub$description_s[i]))
  }
  cat("\n")
}

# --- 3.6 Build aligned predictor matrix ---------------------
available_preds <- predictor_selection$key[predictor_selection$in_json]

cat("Building predictor matrix aligned to ONS dates...\n")

# Pivot to wide format, filter to ONS date range
predictor_wide <- nowcast_long %>%
  filter(key %in% available_preds,
         date >= ons_start,
         date <= ons_end) %>%
  pivot_wider(names_from = key, values_from = value) %>%
  arrange(date)

cat("Predictor matrix dimensions:", dim(predictor_wide), "\n")

# Check NAs per predictor
pred_na <- colSums(is.na(predictor_wide[, -1]))
if (any(pred_na > 0)) {
  cat("\nPredictors with missing values:\n")
  print(pred_na[pred_na > 0])
  cat("→ Will forward-fill in cleaning step\n")
}

# Forward-fill NAs
predictor_clean <- predictor_wide %>%
  mutate(across(-date, ~ zoo::na.locf(.x, na.rm = FALSE))) %>%
  mutate(across(-date, ~ zoo::na.locf(.x, fromLast = TRUE,
                                      na.rm = FALSE)))  # backward fill any start NAs

cat("\nPredictor matrix clean: ready\n")

# --- 3.7 Plot: Key predictors -------------------------------
cat("\nPlotting key predictor series...\n")

plot_keys  <- available_preds[1:min(6, length(available_preds))]
plot_data  <- nowcast_long %>%
  filter(key %in% plot_keys,
         date >= as.Date("2005-01-01"),
         !is.na(value)) %>%
  left_join(metadata[, c("key","description")], by = "key") %>%
  mutate(label = substr(description, 1, 35))

ggplot(plot_data, aes(x = date, y = value)) +
  geom_line(color = "steelblue", linewidth = 0.6) +
  facet_wrap(~label, scales = "free_y", ncol = 2) +
  labs(title    = "Sample Nowcast Predictors",
       subtitle = "Monthly series, 2005–present",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(strip.text = element_text(size = 7))

# --- 3.8 Correlation with ONS total -------------------------
cat("\n--- Correlation of Predictors with ONS Total ---\n")

# Merge predictors with ONS total on common dates
merged <- ons_clean %>%
  select(date, total.total) %>%
  inner_join(predictor_clean, by = "date")

if (nrow(merged) > 10) {
  cors <- cor(merged$total.total,
              merged[, available_preds[available_preds %in% names(merged)]],
              use = "pairwise.complete.obs")
  cors_df <- data.frame(
    key  = colnames(cors),
    cor  = as.numeric(cors)
  ) %>%
    left_join(metadata[, c("key","description")], by = "key") %>%
    mutate(description = substr(description, 1, 50)) %>%
    arrange(desc(abs(cor)))
  
  cat("Top 10 correlated predictors with total.total:\n")
  print(cors_df[1:min(10, nrow(cors_df)), c("key","cor","description")])
}


# ============================================================
# PART 4: KOF FORECAST DATA
# ============================================================

cat("\n============================================================\n")
cat("PART 4: KOF FORECAST DATA\n")
cat("============================================================\n")

# --- 4.1 Download (only if not yet present) -----------------
kof_data_url <- "https://datenservice.kof.ethz.ch/api/v1/public/sets/vja_public_q?mime=xlsx"
kof_meta_url <- "https://datenservice.kof.ethz.ch/api/v1/public/metadata/collections/vja_public_q?mime=xlsx&locale=de"

if (!file.exists(here("data", "kof_forecasts.xlsx"))) {
  cat("Downloading KOF data...\n")
  download.file(kof_data_url,
                destfile = here("data", "kof_forecasts.xlsx"), mode = "wb")
  download.file(kof_meta_url,
                destfile = here("data", "kof_metadata.xlsx"),  mode = "wb")
  cat("Downloaded successfully.\n")
} else {
  cat("KOF files already present.\n")
}

# --- 4.2 Load -----------------------------------------------
kof_data <- read_excel(here("data", "kof_forecasts.xlsx"))
kof_meta <- read_excel(here("data", "kof_metadata.xlsx"))

cat("KOF data dimensions:    ", dim(kof_data), "\n")
cat("KOF metadata dimensions:", dim(kof_meta), "\n")
cat("\nKOF data columns (first 10):\n")
print(names(kof_data)[1:min(10, ncol(kof_data))])
cat("\nKOF metadata columns:\n")
print(names(kof_meta))
cat("\nKOF data preview:\n")
print(head(kof_data[, 1:min(6, ncol(kof_data))]))

View(kof_data)
View(kof_meta)

cat("\nNote: KOF data is quarterly — will need to be\n")
cat("interpolated to monthly or used as quarterly dummies.\n")


# ============================================================
# PART 5: SAVE CLEANED OBJECTS FOR MODELLING
# ============================================================

cat("\n============================================================\n")
cat("PART 5: SAVE CLEANED OBJECTS\n")
cat("============================================================\n")

# Save all cleaned objects as a single .RData file
save(ons_clean,
     predictor_clean,
     predictor_selection,
     available_preds,
     metadata,
     kof_data,
     kof_meta,
     CANTONS,
     DESTINATIONS,
     canton_total_cols,
     origin_cols,
     file = here("data", "cleaned_data.RData"))

cat("Saved: data/cleaned_data.RData\n")
cat("Contains: ons_clean, predictor_clean, predictor_selection,\n")
cat("          available_preds, metadata, kof_data, kof_meta\n")
cat("\nLoad in any other script with:\n")
cat("  load(here('data', 'cleaned_data.RData'))\n")


# ============================================================
# PART 6: SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n\n")

cat("ONS overnight stays:\n")
cat(sprintf("  %d months × %d series, 2005-01 to 2026-01\n",
            nrow(ons_clean), ncol(ons_clean) - 3))
cat("  No missing values ✓\n")
cat("  Strong summer seasonality → seasonal dummies needed\n")
cat("  COVID shock 2020-03 to 2021-12 → dummy added ✓\n\n")

cat("Predictor matrix:\n")
cat(sprintf("  %d series × %d months\n",
            length(available_preds), nrow(predictor_clean)))
cat(sprintf("  Covers: CH, DE, FR, IT, EA, US, CN, global\n\n"))

cat("KOF forecasts:\n")
cat("  Quarterly — needs interpolation to monthly\n\n")

cat("MODELLING PLAN:\n")
cat("  Step 1: VAR-X at Level 0 (national) with ~25 exogenous preds\n")
cat("  Step 2: VAR-X at Level 1 (13 cantons) — parsimony needed\n")
cat("  Step 3: VAR-X at Level 2 (10 origins) — parsimony needed\n")
cat("  Step 4: Hierarchical reconciliation (bottom-up or MinT)\n")
cat("  Step 5: Forecast evaluation vs benchmark\n")
cat("\n============================================================\n")

