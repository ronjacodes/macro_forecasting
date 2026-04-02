# ============================================================
# VAR ANALYSIS — Swiss Export Dynamics
# ons_data_hierarchy_FS26.xlsx
# Monthly data: 2005-01 to 2026-01 (253 observations)
# ============================================================

# --- Libraries ----------------------------------------------
packages <- c("readxl", "vars", "tseries", "urca", "ggplot2",
              "dplyr", "tidyr", "lubridate", "here", "zoo")
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE)


# ============================================================
# 1. LOAD DATA
# ============================================================
df <- read_excel(here("data", "ons_data_hierarchy_FS26.xlsx"),
                 sheet = "hierarchy")

# Parse date column (format: "2005-01")
df$date <- as.Date(paste0(df$date, "-01"))
df <- df[order(df$date), ]

cat("Loaded:", nrow(df), "observations from",
    format(min(df$date), "%Y-%m"), "to",
    format(max(df$date), "%Y-%m"), "\n")
cat("Total columns:", ncol(df), "\n")


# ============================================================
# 2. SELECT VARIABLES FOR VAR
# ============================================================
# Focused selection: Switzerland total + key trading partners
# Keeping it small (6 vars) for a stable VAR with 253 obs

var_cols <- c(
  "total.ch",       # Switzerland total exports
  "total.de",       # Germany
  "total.europa",   # Europe
  "total.china",    # China
  "total.amerika",  # Americas
  "total.asien"     # Asia
)

# Subset and check for NAs
ts_df <- df[, c("date", var_cols)]
cat("\nMissing values per variable:\n")
print(colSums(is.na(ts_df[, var_cols])))

# Forward fill any gaps
ts_df[, var_cols] <- lapply(ts_df[, var_cols], zoo::na.locf, na.rm = FALSE)


# ============================================================
# 3. CONVERT TO TIME SERIES
# ============================================================
ts_data <- ts(
  ts_df[, var_cols],
  start     = c(2005, 1),
  frequency = 12
)

# Plot raw series
plot(ts_data, main = "Swiss Exports by Region (raw)", col = "steelblue")


# ============================================================
# 4. STATIONARITY TESTS (ADF)
# ============================================================
cat("\n===== ADF Tests — Levels =====\n")
adf_levels <- sapply(var_cols, function(v) {
  adf.test(na.omit(ts_data[, v]))$p.value
})
print(round(adf_levels, 4))

# First difference
ts_diff <- diff(ts_data)

cat("\n===== ADF Tests — First Difference =====\n")
adf_diff <- sapply(var_cols, function(v) {
  adf.test(na.omit(ts_diff[, v]))$p.value
})
print(round(adf_diff, 4))

# Use differenced data if series are non-stationary in levels
# (expected for export series — they trend upward over time)
ts_model <- ts_diff

plot(ts_model, main = "Swiss Exports by Region (differenced)", col = "steelblue")


# ============================================================
# 5. LAG SELECTION
# ============================================================
lag_select <- VARselect(ts_model, lag.max = 12, type = "const")
cat("\n===== Lag Selection =====\n")
print(lag_select$selection)

# AIC tends to overfit with small samples — consider BIC (SC) instead
optimal_lag <- lag_select$selection["SC(n)"]   # BIC
cat("Optimal lag (BIC):", optimal_lag, "\n")


# ============================================================
# 6. ESTIMATE VAR MODEL
# ============================================================
var_model <- VAR(ts_model, p = optimal_lag, type = "const")
summary(var_model)


# ============================================================
# 7. DIAGNOSTIC TESTS
# ============================================================
cat("\n===== Serial Correlation (want p > 0.05) =====\n")
print(serial.test(var_model, lags.pt = 12, type = "PT.asymptotic"))

cat("\n===== Normality Test =====\n")
print(normality.test(var_model, multivariate.only = TRUE))

cat("\n===== ARCH Test =====\n")
print(arch.test(var_model, lags.multi = 5))

cat("\n===== Stability (all roots should be < 1) =====\n")
roots_vals <- roots(var_model)
print(round(roots_vals, 4))
cat("Model stable:", all(roots_vals < 1), "\n")


# ============================================================
# 8. GRANGER CAUSALITY
# ============================================================
cat("\n===== Granger Causality =====\n")
for (v in var_cols) {
  cat("\n-- Does", v, "Granger-cause others? --\n")
  result <- tryCatch(causality(var_model, cause = v), error = function(e) NULL)
  if (!is.null(result)) print(result$Granger)
}


# ============================================================
# 9. IMPULSE RESPONSE FUNCTIONS (IRF)
# ============================================================
# Key question: how does a shock to German exports
# affect Swiss total exports and other regions?

irf_result <- irf(
  var_model,
  impulse  = "total.de",       # shock origin: Germany
  response = var_cols,          # response: all variables
  n.ahead  = 24,                # 24 months ahead
  boot     = TRUE,
  ci       = 0.95
)
plot(irf_result, main = "IRF: Shock to German Exports")

# Shock to China
irf_china <- irf(
  var_model,
  impulse  = "total.china",
  response = var_cols,
  n.ahead  = 24,
  boot     = TRUE,
  ci       = 0.95
)
plot(irf_china, main = "IRF: Shock to China Exports")


# ============================================================
# 10. FORECAST ERROR VARIANCE DECOMPOSITION (FEVD)
# ============================================================
fevd_result <- fevd(var_model, n.ahead = 24)

cat("\n===== FEVD at 12-month horizon =====\n")
print(lapply(fevd_result, function(x) round(x[12, ], 3)))

plot(fevd_result, main = "Forecast Error Variance Decomposition")


# ============================================================
# 11. FORECASTING (12 months ahead)
# ============================================================
n_ahead <- 12
fc <- predict(var_model, n.ahead = n_ahead, ci = 0.95)
fanchart(fc, main = "VAR Forecast — Swiss Exports (12 months)")

# Tidy forecast into data frame
last_date <- max(ts_df$date)

fc_df <- bind_rows(lapply(var_cols, function(v) {
  f <- fc$fcst[[v]]
  data.frame(
    variable = v,
    date     = seq.Date(last_date %m+% months(1),
                        by = "month", length.out = n_ahead),
    forecast = f[, "fcst"],
    lower    = f[, "lower"],
    upper    = f[, "upper"]
  )
}))

# ggplot forecast
ggplot(fc_df, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "steelblue") +
  geom_line(aes(y = forecast), color = "steelblue", linewidth = 1) +
  facet_wrap(~variable, scales = "free_y") +
  labs(title = "Swiss Export Forecasts by Region (12 months ahead)",
       x = NULL, y = "Change in Exports (differenced)") +
  theme_minimal()


# ============================================================
# 12. SAVE OUTPUTS
# ============================================================
write.csv(fc_df, here("data", "var_forecasts.csv"), row.names = FALSE)
saveRDS(var_model, here("data", "var_model.rds"))
cat("\nDone. Outputs saved to data/ folder.\n")
