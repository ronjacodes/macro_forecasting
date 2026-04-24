# ============================================================
# analysis/01_variable_selection.R
# Predictor selection for the bridge VAR
# ============================================================
#
# PURPOSE:
#   Select the monthly indicators to include as predictors
#   in the bridge VAR. Documents the selection criteria,
#   checks availability and coverage, confirms stationarity,
#   and plots the final candidates.
#
# SELECTION CRITERIA:
#   1. Economic relevance — does the variable plausibly
#      drive or signal Swiss GDP growth?
#   2. Stationarity — VARs require stationary inputs. Non-
#      stationary series use their precomputed growth rate
#      variant (_pct_1m) from the JSON instead of _lvl.
#   3. Coverage — must start by January 2000 to cover the
#      full estimation sample including the 2008-09 crisis.
#   4. Parsimony — one representative per category to avoid
#      multicollinearity and parameter proliferation.
#
# FINAL SELECTION (5 predictors):
#
#   Swiss domestic activity:
#     swobs085q_lvl    — CH OECD composite leading indicator.
#                        Designed to predict turning points
#                        6-9 months ahead. Stationary in lvl.
#     swcnfbusq_lvl    — KOF business situation survey.
#                        Starts Jan 2004 — constrains the
#                        estimation sample start date.
#                        Stationary in lvl.
#
#   Foreign demand (exports ~70% of Swiss GDP):
#     bdiptot_g_pct_1m — German industrial production MoM%.
#                        Germany is CH's largest trade partner.
#                        Non-stationary in lvl → use pct_1m.
#     ekeusesig_lvl    — EA Economic Sentiment Indicator.
#                        Preferred over empmia_hq (EA PMI)
#                        which only starts 2010. Stationary.
#
#   Financial conditions:
#     swxsfec__pct_1m  — CHF/EUR exchange rate MoM%.
#                        Stronger franc hurts export competit-
#                        iveness. Non-stationary → use pct_1m.
#
#   DROPPED: Swiss PMI — no composite series in the data,
#   only sub-components. OECD CLI captures same information.
#
# KEY FINDINGS:
#   All 5 candidates found in JSON, all stationary. Coverage
#   starts Feb 2000 at earliest (bdiptot_g, swxsfec_) — one
#   month off target, acceptable. swcnfbusq starts Jan 2004,
#   constraining the estimation sample to 2004 Q1.
#
#   Granger causality (bivariate VAR(2), p < 0.05):
#     swobs085q  p = 0.0002  Granger-causes GDP ✓
#     bdiptot_g  p = 0.0003  Granger-causes GDP ✓
#     swxsfec_   p = 0.0188  Granger-causes GDP ✓
#     swcnfbusq  p = 0.3710  no significant result
#     ekeusesig  p = 0.6781  no significant result
#
#   swcnfbusq and ekeusesig are retained despite failing the
#   bivariate Granger test. The bivariate test does not
#   capture their joint contribution in the full VAR, and
#   both have strong economic rationale. Their correlation
#   with swobs085q (CLI) likely absorbs their individual
#   predictive content in the bivariate setting.
#
# ============================================================

library(here)
library(jsonlite)
library(dplyr)
library(ggplot2)
library(tseries)
library(vars)

select <- dplyr::select

source(here("R", "data_utils.R"))
source(here("R", "model_utils.R"))

json_text   <- paste(readLines(here("data", "swiss_nowcast_data.json"),
                               encoding = "UTF-8"), collapse = "")
nowcast_raw <- fromJSON(json_text, flatten = TRUE)


# ============================================================
# 1. DEFINE CANDIDATES
# ============================================================

candidates <- data.frame(
  base_key = c("swobs085q", "swcnfbusq",
               "bdiptot_g", "ekeusesig", "swxsfec_"),
  json_key = c("swobs085q_lvl",    # stationary in levels
               "swcnfbusq_lvl",    # stationary balance statistic
               "bdiptot_g_pct_1m", # non-stationary → MoM growth
               "ekeusesig_lvl",    # stationary in levels
               "swxsfec__pct_1m"), # non-stationary → MoM growth
  label    = c("CH OECD CLI", "KOF Business Situation",
               "DE Ind. Production (MoM%)", "EA Economic Sentiment",
               "CHF/EUR (MoM%)"),
  stringsAsFactors = FALSE
)

series_list <- lapply(seq_len(nrow(candidates)), function(i) {
  s <- pull_series(candidates$json_key[i])
  if (is.null(s)) return(NULL)
  s$label    <- candidates$label[i]
  s$base_key <- candidates$base_key[i]
  s
})
names(series_list) <- candidates$base_key


# ============================================================
# 2. AVAILABILITY AND COVERAGE CHECK
# ============================================================

# Series must start by 2000-01 to cover the full estimation
# sample including the 2008-09 financial crisis, which is
# critical for the model to learn recession dynamics.

cat("--- Availability and coverage check ---\n\n")

for (i in seq_len(nrow(candidates))) {
  k  <- candidates$base_key[i]
  jk <- candidates$json_key[i]
  s  <- series_list[[k]]
  
  if (is.null(s)) {
    cat(sprintf("%-12s  NOT FOUND in JSON (key: %s)\n", k, jk))
    next
  }
  
  start <- min(s$date)
  end   <- max(s$date)
  n_na  <- sum(is.na(s$value))
  ok    <- start <= as.Date("2000-01-01")
  
  cat(sprintf("%-12s  %s to %s  NAs: %d  Covers 2000: %s\n",
              k,
              format(start, "%Y-%m"),
              format(end,   "%Y-%m"),
              n_na,
              ifelse(ok, "YES \u2713", "NO — note sample start")))
}


# ============================================================
# 3. STATIONARITY — ADF tests
# ============================================================

# German IP and CHF/EUR are already in growth rate form here,
# so they should pass even though their levels would not.

cat("\n--- ADF stationarity tests ---\n")
cat("(H0: unit root — p < 0.05 \u2192 stationary)\n\n")

for (i in seq_len(nrow(candidates))) {
  k <- candidates$base_key[i]
  s <- series_list[[k]]
  if (is.null(s)) next
  x <- s$value[!is.na(s$value)]
  p <- tryCatch(adf.test(x)$p.value, error = function(e) NA)
  cat(sprintf("%-12s  p = %.4f  %s\n", k, p,
              ifelse(!is.na(p) & p < 0.05,
                     "stationary \u2713",
                     "NON-STATIONARY — reconsider transformation")))
}


# ============================================================
# 4. GRANGER CAUSALITY
# ============================================================

# Granger causality tests whether past values of a predictor
# contain information about future GDP growth beyond what
# GDP's own lags already provide. A significant result
# (p < 0.05) is statistical confirmation that the variable
# has predictive content — supporting the economic argument
# for including it.
#
# We test using a bivariate VAR(2) for each predictor vs GDP,
# estimated on the full available sample for each pair.
# H0: predictor does NOT Granger-cause GDP growth.

cat("\n--- Granger causality tests ---\n")
cat("(H0: predictor does not Granger-cause GDP)\n")
cat("(p < 0.05 \u2192 reject H0 \u2192 predictor has predictive content)\n\n")

gdp <- pull_series("ch_kof_modelinput_gdpos_pct_3m") %>%
  rename(gdp = value)

for (i in seq_len(nrow(candidates))) {
  k <- candidates$base_key[i]
  s <- series_list[[k]]
  if (is.null(s)) next
  
  # Merge GDP and predictor on common quarterly dates
  # GDP is quarterly so we aggregate the predictor first
  pred_q <- to_quarterly(s, k)
  merged  <- gdp %>%
    rename(qdate = date) %>%
    inner_join(pred_q, by = "qdate") %>%
    filter(!is.na(gdp), !is.na(.data[[k]]))
  
  if (nrow(merged) < 20) {
    cat(sprintf("%-12s  insufficient data\n", k))
    next
  }
  
  # Bivariate VAR(2) — small enough to be fast, large enough
  # to capture typical quarterly dynamics
  bvar <- tryCatch(
    VAR(merged[, c("gdp", k)], p = 2, type = "const"),
    error = function(e) NULL
  )
  if (is.null(bvar)) next
  
  # causality() tests whether the predictor block of
  # coefficients in the GDP equation is jointly zero
  gc <- tryCatch(
    causality(bvar, cause = k)$Granger,
    error = function(e) NULL
  )
  if (is.null(gc)) next
  
  p_val <- gc$p.value
  cat(sprintf("%-12s  F = %6.3f  p = %.4f  %s\n",
              k, gc$statistic, p_val,
              ifelse(p_val < 0.05,
                     "Granger-causes GDP \u2713",
                     "no significant predictive content")))
}


# ============================================================
# 5. PLOT CANDIDATES
# ============================================================

# Shaded regions mark the 2008-09 recession and 2020 COVID
# crash — a good nowcasting predictor should show clear
# movement in both episodes.

plot_df <- bind_rows(series_list) %>%
  filter(!is.na(value), date >= as.Date("2000-01-01"))

ggplot(plot_df, aes(x = date, y = value)) +
  geom_line(color = "steelblue", linewidth = 0.5) +
  geom_rect(
    data = data.frame(
      xmin = as.Date(c("2008-09-01", "2020-03-01")),
      xmax = as.Date(c("2009-06-01", "2020-09-01")),
      ymin = -Inf, ymax = Inf),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = "red", alpha = 0.08) +
  facet_wrap(~ label, scales = "free_y", ncol = 2) +
  labs(title    = "Final predictor candidates",
       subtitle = "Shaded: 2008-09 recession and 2020 COVID shock",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 9))