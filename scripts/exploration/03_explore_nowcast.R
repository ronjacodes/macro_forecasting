# ============================================================
# exploration/03_explore_nowcast.R
# Explore the full nowcast dataset: plots, coverage, stationarity
# ============================================================
#
# PURPOSE:
#   Get a broad overview of all 133 monthly level series —
#   what they look like, how long they are, and whether they
#   are stationary. This informs predictor selection.
#
# REQUIRES:
#   data/cleaned_nowcast.RData — produced by 01_load_nowcast.R
#   (nowcast_long, coverage, metadata)
#
# KEY FINDINGS:
#   - 92 / 133 monthly series are stationary at the 5% level
#   - 41 non-stationary: mainly exchange rates, money supply,
#     stock index, some survey series
#   - 25 series start after 2005 — constrain sample if used
#   - Switzerland dominates with 72 series
#   - COVID crash (2020 Q2) visible in all business surveys
#
# ============================================================

library(here)
library(dplyr)
library(ggplot2)
library(tseries)

load(here("data", "cleaned_nowcast.RData"))

# Focus on monthly series — quarterly and daily series are
# handled separately and are not used as bridge predictors
monthly_long <- nowcast_long %>% filter(frequency == "Month")
monthly_keys <- unique(monthly_long$base_key)
cat("Monthly series:", length(monthly_keys), "\n")


# ============================================================
# 1. PLOTS BY COUNTRY / REGION
# ============================================================

# Extract country label from description (text before first comma)
monthly_long <- monthly_long %>%
  mutate(country = sub(",.*", "", description))

cat("\nSeries count by country/region:\n")
print(sort(table(
  monthly_long %>% distinct(base_key, country) %>% pull(country)
), decreasing = TRUE))

countries <- monthly_long %>%
  distinct(base_key, country) %>%
  pull(country) %>% unique() %>% sort()

# For countries with many series (>6) use facets with free
# y-axis — overlaying them produces unreadable spaghetti plots
# because series are on incompatible scales (e.g. SMI at
# 10,000+ vs survey balances at ±50).
for (ctry in countries) {
  plot_data <- monthly_long %>% filter(country == ctry)
  n_series  <- length(unique(plot_data$base_key))

  if (n_series > 6) {
    p <- ggplot(plot_data, aes(x = date, y = value)) +
      geom_line(linewidth = 0.4, color = "steelblue") +
      facet_wrap(~ base_key, scales = "free_y") +
      labs(title    = ctry,
           subtitle = paste(n_series,
                            "series (level) — free y-axis per panel"),
           x = NULL, y = NULL) +
      theme_minimal(base_size = 7) +
      theme(strip.text = element_text(size = 6),
            plot.title = element_text(face = "bold"))
  } else {
    p <- ggplot(plot_data,
                aes(x = date, y = value,
                    color = base_key, group = base_key)) +
      geom_line(linewidth = 0.5, alpha = 0.8) +
      labs(title    = ctry,
           subtitle = paste(n_series, "series (level)"),
           x = NULL, y = NULL, color = NULL) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom",
            legend.text     = element_text(size = 7),
            plot.title      = element_text(face = "bold"))
  }
  print(p)
}


# ============================================================
# 2. COVERAGE PLOT
# ============================================================

# Shows when each series starts and ends — useful for deciding
# the estimation sample start date. Series starting later
# force a shorter sample if included.

monthly_coverage <- coverage %>% filter(frequency == "Month")

ggplot(monthly_coverage,
       aes(x = start, xend = end,
           y  = reorder(base_key, start),
           yend = reorder(base_key, start))) +
  geom_segment(linewidth = 1.5, color = "steelblue", alpha = 0.7) +
  labs(title    = "Series coverage",
       subtitle = "Monthly level series — sorted by start date",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.y = element_text(size = 5))


# ============================================================
# 3. STATIONARITY — ADF tests
# ============================================================

# VAR models require stationary variables. We test all monthly
# level series with the Augmented Dickey-Fuller test.
# H0: series has a unit root (non-stationary).
# p < 0.05 → reject H0 → series is stationary.
#
# Non-stationary series can still be used if we switch to a
# stationary transformation (e.g. _pct_1m or _dif_1y variant).

cat("\n--- ADF stationarity tests (level series) ---\n")

adf_results <- lapply(monthly_keys, function(k) {
  x <- monthly_long %>%
    filter(base_key == k) %>%
    arrange(date) %>%
    pull(value)
  x <- x[!is.na(x)]
  if (length(x) < 20) {
    return(data.frame(base_key = k, adf_p = NA,
                      stationary = NA, n = length(x)))
  }
  p <- tryCatch(adf.test(x)$p.value, error = function(e) NA)
  data.frame(base_key   = k,
             adf_p      = round(p, 4),
             stationary = ifelse(!is.na(p), p < 0.05, NA),
             n          = length(x))
})

adf_results <- bind_rows(adf_results) %>%
  left_join(metadata %>% select(key, description, source),
            by = c("base_key" = "key")) %>%
  arrange(adf_p)

cat("Stationary (p < 0.05):    ",
    sum(adf_results$stationary, na.rm = TRUE), "\n")
cat("Non-stationary (p >= 0.05):",
    sum(!adf_results$stationary, na.rm = TRUE), "\n")
cat("Could not test:            ",
    sum(is.na(adf_results$stationary)), "\n")

cat("\nNon-stationary series:\n")
ns <- adf_results %>%
  filter(!stationary) %>%
  select(base_key, adf_p, n, description)
print(as.data.frame(ns), row.names = FALSE)


# ============================================================
# 4. SHORT SERIES — flag series starting after 2005
# ============================================================

# Our target estimation sample starts 2004 Q1. Series starting
# after 2005 would shorten the sample if included and may miss
# the 2008-09 financial crisis — a key episode for the model
# to learn recession dynamics.

cat("\n--- Series starting after 2005-01 ---\n")
short <- monthly_coverage %>%
  filter(start > as.Date("2005-01-01")) %>%
  left_join(metadata %>% select(key, description),
            by = c("base_key" = "key")) %>%
  arrange(start) %>%
  select(base_key, start, end, n_obs, description)

print(short, n = 50)
