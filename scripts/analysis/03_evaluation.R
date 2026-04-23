# ============================================================
# analysis/03_evaluation.R
# Pseudo out-of-sample forecast evaluation — bridge VAR
# ============================================================
#
# PURPOSE:
#   Evaluate the bridge VAR's forecasting accuracy in a
#   realistic pseudo out-of-sample setting and compare
#   against a naive AR(1) benchmark.
#
# APPROACH:
#   For each quarter in the evaluation window (2015 Q1 to
#   2025 Q4, 44 quarters) we:
#     1. Estimate the VAR on all prior data (fully recursive
#        — model is re-estimated at each step)
#     2. Produce three nowcast vintages via the bridge step:
#          M1: average of first month of the quarter
#          M2: average of first two months
#          M3: average of all three months
#     3. Compare to actual GDP and AR(1) benchmark
#
# WHY THE GDP EQUATION DIRECTLY (bridge_nowcast):
#   vars::predict() forecasts all variables jointly from
#   lagged values only — it cannot condition on new within-
#   quarter monthly data. bridge_nowcast() extracts the GDP
#   equation coefficients and plugs in the partial-quarter
#   predictor averages, making M1/M2/M3 genuinely different.
#
# KEY RESULTS:
#   Full sample (incl. COVID):  rel. RMSE = 0.889  VAR wins
#   Excluding COVID:            rel. RMSE = 2.315  AR(1) wins
#   COVID quarters only:        rel. RMSE = 0.526  VAR wins
#
#   The VAR's advantage is concentrated in the COVID episode.
#   In tranquil periods the AR(1) is more accurate — a known
#   result in the macro forecasting literature (GDP growth is
#   hard to predict in stable times). Extensions (BVAR, MFVAR)
#   aim to improve performance outside crisis periods.
#
# ============================================================

library(here)
library(jsonlite)
library(dplyr)
library(tidyr)
library(vars)
library(ggplot2)

select <- dplyr::select

source(here("R", "data_utils.R"))
source(here("R", "model_utils.R"))

json_text   <- paste(readLines(here("data", "swiss_nowcast_data.json"),
                               encoding = "UTF-8"), collapse = "")
nowcast_raw <- fromJSON(json_text, flatten = TRUE)


# ============================================================
# 1. DATA
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

pred_names    <- names(predictors)
quarterly     <- build_quarterly_panel(gdp, predictors)
monthly_panel <- build_monthly_panel(predictors)


# ============================================================
# 2. PSEUDO OUT-OF-SAMPLE LOOP
# ============================================================

EVAL_START <- as.Date("2015-01-01")
EVAL_END   <- as.Date("2025-10-01")
MIN_TRAIN  <- as.Date("2004-01-01")

eval_quarters <- quarterly %>%
  filter(qdate >= EVAL_START, qdate <= EVAL_END) %>%
  pull(qdate)

cat("Evaluation quarters:", length(eval_quarters), "\n\n")

results <- list()

for (tgt in eval_quarters) {
  tgt_date <- as.Date(tgt, origin = "1970-01-01")

  train <- quarterly %>%
    filter(qdate >= MIN_TRAIN, qdate < tgt_date) %>%
    filter(!is.na(gdp) & !is.na(cli) & !is.na(kof) &
             !is.na(de_ip) & !is.na(ea_esi) & !is.na(chfeur))

  if (nrow(train) < 20) next

  actual_gdp <- quarterly$gdp[quarterly$qdate == tgt_date]
  if (is.na(actual_gdp)) next

  # Estimate VAR and AR(1) on training data
  fit    <- fit_var(train, pred_names, lag_max = 2)
  var_mod <- fit$model
  p       <- fit$p
  if (is.null(var_mod)) next

  ar_fc <- ar1_forecast(train$gdp)

  # Bridge nowcasts — M1, M2, M3 differ because each uses a
  # different number of months from the target quarter
  var_nowcasts <- sapply(1:3, function(m) {
    avgs <- tryCatch(
      get_partial_avg(monthly_panel, pred_names, tgt_date, m),
      error = function(e) NULL
    )
    if (is.null(avgs) || nrow(avgs) == 0) return(NA)
    tryCatch(
      bridge_nowcast(var_mod, train, avgs, pred_names, p),
      error = function(e) NA
    )
  })

  results[[as.character(tgt_date)]] <- data.frame(
    qdate  = tgt_date,
    actual = actual_gdp,
    var_m1 = var_nowcasts[1],
    var_m2 = var_nowcasts[2],
    var_m3 = var_nowcasts[3],
    ar1    = ar_fc,
    covid  = as.integer(tgt_date %in%
               as.Date(c("2020-04-01", "2020-07-01")))
  )
}

eval_df <- bind_rows(results)
cat("Quarters with results:    ", nrow(eval_df), "\n")
cat("Quarters with M3 nowcast:", sum(!is.na(eval_df$var_m3)), "\n\n")
print(eval_df)


# ============================================================
# 3. ACCURACY METRICS
# ============================================================

compute_metrics(eval_df, "Full sample 2015 Q1 - 2025 Q4")
compute_metrics(filter(eval_df, covid == 0), "Excluding COVID")
compute_metrics(filter(eval_df, covid == 1), "COVID quarters only")


# ============================================================
# 4. KOF-STYLE PLOT: NOWCAST EVOLUTION BY VINTAGE
# ============================================================

# Shows how the nowcast for each quarter converges as more
# monthly data arrives. Inspired by the KOF Nowcasting Lab
# forecast history chart. COVID quarters excluded for clarity.

plot_df <- eval_df %>%
  filter(covid == 0) %>%
  tidyr::pivot_longer(c(var_m1, var_m2, var_m3),
                      names_to  = "vintage",
                      values_to = "nowcast") %>%
  mutate(vintage = recode(vintage,
                          var_m1 = "Month 1",
                          var_m2 = "Month 2",
                          var_m3 = "Month 3"))

recent <- sort(unique(eval_df$qdate[eval_df$covid == 0]),
               decreasing = TRUE)[1:12]

ggplot(plot_df %>% filter(qdate %in% recent),
       aes(x = qdate)) +
  geom_line(aes(y = actual), color = "black", linewidth = 0.9) +
  geom_line(aes(y = nowcast, color = vintage, group = vintage),
            linewidth = 0.6) +
  geom_point(aes(y = nowcast, color = vintage), size = 2) +
  scale_color_manual(
    values = c("Month 1" = "#91bfdb",
               "Month 2" = "#4575b4",
               "Month 3" = "#d73027"),
    name = "Information vintage"
  ) +
  geom_hline(yintercept = 0, linetype = "dotted",
             color = "gray50") +
  labs(title    = "Nowcast evolution by information vintage",
       subtitle = paste("Black = actual GDP ex-sports (QoQ%)",
                        "\nColours = bridge VAR after 1, 2, 3 months"),
       x = NULL, y = "QoQ growth (%)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")


# ============================================================
# 5. FORECAST ERROR PLOT
# ============================================================

# Shows whether errors are systematic or random over time.
# Persistent sign in errors suggests a structural miss.

ggplot(eval_df %>%
         mutate(err_var = actual - var_m3,
                err_ar  = actual - ar1),
       aes(x = qdate)) +
  geom_col(aes(y = err_var), fill = "steelblue",
           alpha = 0.7, width = 60) +
  geom_line(aes(y = err_ar), color = "firebrick",
            linewidth = 0.6, linetype = "dashed") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(title    = "Forecast errors over time",
       subtitle = paste("Blue bars = VAR M3 error",
                        " | Red dashed = AR(1) error",
                        "\nPositive = underestimate,",
                        "negative = overestimate"),
       x = NULL, y = "Forecast error (pp)") +
  theme_minimal(base_size = 10)
