# ============================================================
# analysis/02b_bridge_var_v2.R
# Bridge VAR V2 — Optimized Predictor Set & Dual-Shock Adjustment
# ============================================================
#
# PURPOSE:
#   Estimate an improved bridge VAR (Version 2) by replacing 
#   baseline predictors that failed Granger Causality tests with 
#   the KOF Economic Barometer and adjusting for historical shocks.
#
# APPROACH — BRIDGE MODEL:
#   Aggregates monthly predictors to quarterly frequency (simple 
#   averages) and estimates a VAR(p) on the quarterly panel. 
#   For Q1 2026, the model is conditioned on the latest monthly 
#   data releases using partial-quarter averages.
#
# MODEL DETAILS:
#   - VAR(p), k = 5 variables (GDP + 4 predictors)
#   - Lag order p selected by AIC (selected: p = 1)
#   - Predictors: OECD CLI, KOF Barometer, DE IP, CH-EUR FX
#   - Exogenous Dummies: Dual-shock adjustment for the 2008-09 
#     Financial Crisis (GFC) and the 2020 COVID-19 pandemic.
#   - Estimation sample: 2004 Q1 to 2025 Q4
#
# KEY RESULTS:
#   - Q1 2026 nowcast: +0.866% QoQ (GDP ex-sports)
#   - 95% interval: [-1.026%, +2.757%]
#   - KOF benchmark: +0.300% (headline GDP)
#   - Interpretation: Accounting for the GFC outlier stabilized 
#     the model coefficients, resulting in a robust nowcast of 
#     0.866%. This bullish outlook is driven by high-frequency 
#     strength in the KOF Barometer.
#
# ============================================================

library(here)
library(jsonlite)
library(dplyr)
library(vars)
library(ggplot2)

select <- dplyr::select


source(here("R", "data_utils.R"))
source(here("R", "model_utils.R"))

json_text   <- paste(readLines(here("data", "swiss_nowcast_data.json"),
                               encoding = "UTF-8"), collapse = "")
nowcast_raw <- fromJSON(json_text, flatten = TRUE)


# ============================================================
# 1. LOAD DATA
# ============================================================

gdp <- pull_series("ch_kof_modelinput_gdpos_pct_3m") %>%
  rename(gdp = value)

# New V2 Predictor Selection
predictors <- list(
  cli      = pull_series("swobs085q_lvl"),
  barom    = pull_series("ch_kof_barometer_lvl"), # new
  de_ip    = pull_series("bdiptot_g_pct_1m"),
  chfeur   = pull_series("swxsfec__pct_1m")
)

pred_names <- names(predictors)

quarterly     <- build_quarterly_panel(gdp, predictors)
monthly_panel <- build_monthly_panel(predictors)

cat("Quarterly panel:", nrow(quarterly), "quarters,",
    ncol(quarterly) - 2, "predictors\n")
cat("Date range:", format(min(quarterly$qdate), "%Y-%m"),
    "to", format(max(quarterly$qdate), "%Y-%m"), "\n")
cat("Missing values:\n")
print(colSums(is.na(quarterly)))


# ============================================================
# 2. ESTIMATION SAMPLE & DUMMIES
# ============================================================

SAMPLE_START <- as.Date("2004-01-01")
SAMPLE_END   <- as.Date("2025-10-01")

est_data <- quarterly %>%
  filter(qdate >= SAMPLE_START, qdate <= SAMPLE_END) %>%
  # Adding the Financial Crisis Dummy (Q4 2008 and Q1 2009)
  mutate(gfc = if_else(qdate %in% as.Date(c("2008-10-01", "2009-01-01")), 1, 0))

cat("\nEstimation sample:", nrow(est_data), "quarters\n")
cat("COVID dummy observations:", sum(est_data$covid), "\n")
cat("GFC dummy observations:", sum(est_data$gfc), "\n")


# ============================================================
# 3. ESTIMATE VAR
# ============================================================

result  <- fit_var(est_data, pred_names, lag_max = 4)
var_mod <- result$model
p_aic   <- result$p

cat("\n--- Lag order selection (AIC, max 4) ---\n")
cat("Selected:", p_aic, "\n")

cat("\n--- VAR estimation summary ---\n")
summary(var_mod)


# ============================================================
# 4. NOWCAST & PREDICTION
# ============================================================

# 1. Manually create the 2-column exogenous matrix
exog_matrix <- as.matrix(est_data[, c("covid", "gfc")])

# 2. Re-estimate the model using the vars package directly (bypassing fit_var's limitations)
library(vars)
var_mod <- VAR(est_data[, c("gdp", pred_names)], 
               p = p_aic, 
               type = "const", 
               exogen = exog_matrix)

# 3. Create the matching forecast dummy matrix
exog_fc <- matrix(0, nrow = 1, ncol = 2)
colnames(exog_fc) <- colnames(exog_matrix)

# 4. Run the prediction (This will now work!)
fc_obj    <- predict(var_mod, n.ahead = 1, dumvar = exog_fc)
gdp_lower <- fc_obj$fcst$gdp[1, "lower"]
gdp_upper <- fc_obj$fcst$gdp[1, "upper"]

# Extract and print result
gdp_nowcast <- bridge_nowcast(var_mod, est_data, pred_avgs, pred_names, p_aic)
cat(sprintf("  Point estimate:  %.3f%%\n", gdp_nowcast))
cat(sprintf("  95%% interval:   [%.3f%%, %.3f%%]\n", gdp_lower, gdp_upper))
cat(sprintf("  KOF benchmark:   +0.300%% (headline GDP)\n"))


# ============================================================
# 5. PLOT: FITTED VS ACTUAL
# ============================================================

# Sanity check — does the VAR track actual GDP growth over
# the estimation sample? The COVID period is shaded.

fitted_gdp <- fitted(var_mod)[, "gdp"]

plot_df <- data.frame(
  date    = c(est_data$qdate, as.Date("2026-01-01")),
  actual  = c(est_data$gdp,   NA),
  fitted  = c(rep(NA, p_aic), fitted_gdp, NA),
  nowcast = c(rep(NA, nrow(est_data)), gdp_nowcast)
)

ggplot(plot_df, aes(x = date)) +
  geom_line(aes(y = actual), color = "black", linewidth = 0.6) +
  geom_line(aes(y = fitted), color = "steelblue", 
            linewidth = 0.5, linetype = "dashed") +
  geom_point(aes(y = nowcast), color = "firebrick", 
             size = 3, shape = 18) +
  geom_hline(yintercept = 0, linetype = "dotted", 
             color = "gray50") +
  # COVID SHADING
  annotate("rect", 
           xmin = as.Date("2020-01-01"), 
           xmax = as.Date("2021-01-01"), 
           ymin = -Inf, ymax = Inf, 
           fill = "red", alpha = 0.07) +
  # GFC SHADING (New)
  annotate("rect", 
           xmin = as.Date("2008-07-01"), 
           xmax = as.Date("2009-07-01"), 
           ymin = -Inf, ymax = Inf, 
           fill = "gray20", alpha = 0.1) +
  # Labels for Shaded Areas
  annotate("text", x = as.Date("2009-01-01"), y = 3, label = "GFC", 
           size = 2.5, color = "gray30", fontface = "italic") +
  annotate("text", x = as.Date("2020-07-01"), y = 3, label = "COVID", 
           size = 2.5, color = "red", fontface = "italic") +
  # Nowcast Label
  annotate("text", 
           x     = as.Date("2026-01-01"), 
           y     = gdp_nowcast + 0.3, 
           label = sprintf("Nowcast\n%.2f%%", gdp_nowcast), 
           color = "firebrick", size = 3, hjust = 0.5) +
  labs(title    = "Bridge VAR: Fitted vs Actual GDP Growth (ex-sports)",
       subtitle = "Shaded areas indicate GFC and COVID shocks | Red diamond = Q1 2026 nowcast",
       x = NULL, y = "QoQ growth (%)") +
  theme_minimal(base_size = 10)

