# ============================================================
# analysis/02_bridge_var.R
# Bridge VAR — estimation and Q1 2026 nowcast
# ============================================================
#
# PURPOSE:
#   Estimate the baseline bridge VAR on the full sample and
#   produce a nowcast for Q1 2026.
#
# APPROACH — BRIDGE MODEL:
#   Solves the mixed-frequency problem by aggregating monthly
#   predictors to quarterly frequency (simple average within
#   each quarter) before estimation. A standard VAR is then
#   estimated on the quarterly panel. The "bridge" is the
#   link between the available monthly data and the quarterly
#   GDP target — we use partial-quarter averages for Q1 2026
#   to condition the nowcast on the latest monthly releases.
#
# MODEL DETAILS:
#   - VAR(p), k = 6 variables (GDP + 5 predictors)
#   - Lag order p selected by AIC (selected: p = 2)
#   - COVID dummy for 2020 Q2 and Q3 as exogenous variable
#     (absorbs the pandemic shock without being modelled)
#   - Estimation sample: 2004 Q1 to 2025 Q4
#     (constrained by swcnfbusq starting January 2004)
#
# KEY RESULTS:
#   - Q1 2026 nowcast: +1.07% QoQ (GDP ex-sports)
#   - 95% interval: [-0.74%, +2.88%]
#   - KOF benchmark: +0.30% (headline GDP)
#   - Gap reflects KOF's access to more timely/granular data
#     and their richer dynamic factor model. KOF's estimate
#     lies within our uncertainty interval.
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

predictors <- list(
  cli    = pull_series("swobs085q_lvl"),
  kof    = pull_series("swcnfbusq_lvl"),
  de_ip  = pull_series("bdiptot_g_pct_1m"),
  ea_esi = pull_series("ekeusesig_lvl"),
  chfeur = pull_series("swxsfec__pct_1m")
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
# 2. ESTIMATION SAMPLE
# ============================================================

SAMPLE_START <- as.Date("2004-01-01")
SAMPLE_END   <- as.Date("2025-10-01")

est_data <- quarterly %>%
  filter(qdate >= SAMPLE_START, qdate <= SAMPLE_END)

cat("\nEstimation sample:", nrow(est_data), "quarters\n")
cat("COVID dummy observations:", sum(est_data$covid), "\n")
cat("Missing values:\n")
print(colSums(is.na(est_data)))


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
# 4. Q1 2026 NOWCAST
# ============================================================

# Average available monthly data for Q1 2026 (Jan + Feb for
# most series, Jan only for some). Use these as current-quarter
# predictor values in the bridge step.

pred_avgs <- get_partial_avg(monthly_panel, pred_names,
                             as.Date("2026-01-01"), n_months = 2)

cat("\n--- Q1 2026 predictor averages (Jan-Feb 2026) ---\n")
print(as.data.frame(pred_avgs))

gdp_nowcast <- bridge_nowcast(var_mod, est_data,
                              pred_avgs, pred_names, p_aic)

# Confidence interval from VAR predict() for reference
exog_fc  <- matrix(0, nrow = 1, ncol = 1,
                   dimnames = list(NULL, "covid"))
fc_obj   <- predict(var_mod, n.ahead = 1, dumvar = exog_fc)
gdp_lower <- fc_obj$fcst$gdp[1, "lower"]
gdp_upper <- fc_obj$fcst$gdp[1, "upper"]

cat("\n--- Q1 2026 GDP nowcast (ex-sports, QoQ%) ---\n")
cat(sprintf("  Point estimate:  %.3f%%\n", gdp_nowcast))
cat(sprintf("  95%% interval:   [%.3f%%, %.3f%%]\n",
            gdp_lower, gdp_upper))
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
  annotate("rect",
           xmin = as.Date("2020-01-01"),
           xmax = as.Date("2021-01-01"),
           ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.07) +
  annotate("text",
           x     = as.Date("2026-01-01"),
           y     = gdp_nowcast + 0.3,
           label = sprintf("Nowcast\n%.2f%%", gdp_nowcast),
           color = "firebrick", size = 3, hjust = 0.5) +
  labs(title    = "Bridge VAR: fitted vs actual GDP growth (ex-sports)",
       subtitle = paste("Black = actual  |  Dashed blue = fitted",
                        "  |  Red diamond = Q1 2026 nowcast"),
       x = NULL, y = "QoQ growth (%)") +
  theme_minimal(base_size = 10)
