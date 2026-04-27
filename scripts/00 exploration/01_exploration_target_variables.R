# ============================================================================
# 01_exploration.R
# Data exploration for the Macroeconomic Forecasting seminar (FS26, ETH/KOF)
#
# PROJECT CONTEXT:
#   We estimate VAR models to forecast three Swiss macroeconomic variables:
#   (1) Real GDP growth (QoQ %, KOF sports-event corrected)
#   (2) CPI inflation (QoQ %, headline)
#   (3) 10Y government bond yield (QoQ change in pp, first-differenced for stationarity)
#   Models vary by lag order (p=1..5), sample period, variable set, and prior
#   (frequentist VAR vs Bayesian VAR with Minnesota prior).
#   Evaluation via expanding-window pseudo out-of-sample RMSE.
#
# THIS SCRIPT:
#   - Multi-panel overview plots per variable (QoQ, YoY, scatter)
#   - Summary statistics (N, mean, SD, min, max, ADF p-value)
#   - Formal outlier detection (3×IQR rule)
#   - NA checks
#   - Granger causality tests between the three VAR target variables
#
# DATA TIMELINES (confirmed):
#   GDP QoQ  — KOF sports-corrected : Q2 2000 – Q4 2025  (103 quarterly obs)
#   GDP QoQ  — SECO SA              : Q2 2000 – Q4 2025  (103 quarterly obs)
#   GDP YoY  — both variants        : Q1 2001 – Q4 2025  (100 quarterly obs)
#   CPI MoM  — headline             : Feb 2000 – Jan 2026 (312 monthly obs)
#   CPI QoQ  — headline             : Q2 2000 – Q1 2026  (310 monthly obs)
#   CPI YoY  — headline             : Jan 2001 – Jan 2026 (301 monthly obs)
#   CPI MoM  — core                 : Aug 2000 – Jan 2026 (306 monthly obs)
#   CPI YoY  — core                 : Jul 2001 – Jan 2026 (295 monthly obs)  [NOT stationary]
#   Bond lvl — 10Y yield            : Jan 2000 – Feb 2026 (314 monthly obs)  [NOT stationary]
#   Bond dif — MoM change           : Feb 2000 – Feb 2026 (313 monthly obs)
#   Bond dif — QoQ change (summed)  : Q2 2000 – Q1 2026  (104 quarterly obs)
#   VAR dataset (QoQ intersection)  : Q2 2000 – Q4 2025  (103 quarterly obs)
#
# KEY FINDINGS FROM EXPLORATION:
#   Stationarity:
#     - GDP QoQ: stationary (ADF p≤0.01) ✓
#     - CPI QoQ headline: stationary (ADF p≤0.01) ✓
#     - CPI YoY core: NOT stationary (ADF p=0.28) — use QoQ instead
#     - Bond yield level: NOT stationary (ADF p=0.59) — use first difference ✓
#     - Bond yield QoQ change: stationary (ADF p≤0.01) ✓
#   Outliers (3×IQR, GDP only — CPI and bond yield clean):
#     - Q4 2008: GDP −2.1% QoQ (Lehman shock)
#     - Q2 2020: GDP −6.35% QoQ / −6.85% YoY (COVID lockdown, ~6 SD from mean)
#     - Q3 2020: GDP +6.11% QoQ (immediate rebound)
#     - Q2 2021: GDP +10.5% YoY (base effect vs COVID trough)
#     → COVID dummies for 2020 Q1 – 2021 Q4 are essential
#   No NAs in any series.
#   CPI choice: headline QoQ (not core) — SNB target, stationary, matches frequency
#
# Granger causality results (p=4 lags, full sample Q2 2000 – Q4 2025):
#   GDP  → CPI  : p = 0.845  — not significant
#   GDP  → Bond : p = 0.427  — not significant
#   CPI  → GDP  : p = 0.326  — not significant
#   CPI  → Bond : p = 0.206  — not significant
#   Bond → GDP  : p = 0.042  — ** significant at 5%
#   Bond → CPI  : p = 0.004  — *** significant at 1%
# Interpretation: Bond yield changes lead both GDP and CPI — consistent
#   with monetary transmission. GDP and CPI do not Granger-cause each other
#   or bond yields at conventional significance levels. This supports including
#   the bond yield as an informative predictor in the VAR.
#
# Prereq: source("00_setup.R")
# Outputs (commented out — uncomment to save):
#   output/figures/01a_gdp.png, 01b_cpi.png, 01c_bond.png
#   output/tables/01_summary_stats.csv, 01_outliers.csv
# ============================================================================

# ── Helper: two-line panel ────────────────────────────────────────────────────
# Draws one panel with two overlapping lines (s2 / col2 / label2 are optional)
panel <- function(s1, s2 = NULL,
                  col1, col2 = NULL,
                  label1, label2 = NULL,
                  title = "", y_lab = "%",
                  date_min = as.Date("2000-01-01"),
                  zero_line = TRUE) {
  
  s1f <- s1 %>% filter(date >= date_min)
  
  p <- ggplot() +
    { if (zero_line)
      geom_hline(yintercept = 0, linetype = "dashed",
                 color = COL_GREY, linewidth = 0.35) }
  
  # SA / core overlay drawn first (behind)
  if (!is.null(s2) && !is.null(col2) && !is.null(label2)) {
    s2f <- s2 %>% filter(date >= date_min)
    p <- p + geom_line(data = s2f,
                       aes(x = date, y = value, colour = label2),
                       linewidth = 0.65, alpha = 0.85)
  }
  
  p <- p +
    geom_line(data = s1f,
              aes(x = date, y = value, colour = label1),
              linewidth = 0.75) +
    scale_colour_manual(
      values = setNames(c(col1, col2),   c(label1, label2)),
      breaks = c(label1, label2),        # legend order: primary first
      name   = NULL
    ) +
    labs(title = title, x = NULL, y = y_lab) +
    scale_x_date(date_breaks = "4 years", date_labels = "%Y")
  
  p
}

# ── Helper: ADF p-value (handles tseries truncation at 0.01) ─────────────────
adf_p <- function(x) {
  tryCatch(adf.test(na.omit(x))$p.value, error = function(e) NA_real_)
}

# ── Helper: summary stats for one tibble with a value column ─────────────────
sumstat <- function(df, label) {
  tibble(
    Series  = label,
    N       = sum(!is.na(df$value)),
    From    = format(min(df$date, na.rm = TRUE)),
    To      = format(max(df$date, na.rm = TRUE)),
    Mean    = round(mean(df$value,   na.rm = TRUE), 3),
    SD      = round(sd(df$value,     na.rm = TRUE), 3),
    Min     = round(min(df$value,    na.rm = TRUE), 3),
    Max     = round(max(df$value,    na.rm = TRUE), 3),
    `ADF p` = round(adf_p(df$value), 3)
  )
}

# ============================================================================
# 1.  GDP — three-panel: QoQ, YoY, YoY diff
#     Primary (dark blue) = KOF sports-corrected
#     Overlay (light blue) = SECO seasonally adjusted
# ============================================================================
message("Plotting GDP...")

p_gdp_qoq <- panel(
  s1 = gdp$pct_3m,    col1 = COL_GDP,    label1 = "KOF sports-corrected",
  s2 = gdp$sa_pct_3m, col2 = COL_GDP_SA, label2 = "SECO SA",
  title = "Quarter-on-quarter %", y_lab = "%"
)

p_gdp_yoy <- panel(
  s1 = gdp$pct_1y,    col1 = COL_GDP,    label1 = "KOF sports-corrected",
  s2 = gdp$sa_pct_1y, col2 = COL_GDP_SA, label2 = "SECO SA",
  title = "Year-on-year %", y_lab = "%"
)

# Third panel: KOF vs SECO SA scatter to show how close the two series are
gdp_compare <- inner_join(
  gdp$pct_3m    %>% rename(kof = value),
  gdp$sa_pct_3m %>% rename(seco = value),
  by = "date"
)
p_gdp_scatter <- ggplot(gdp_compare, aes(x = seco, y = kof)) +
  geom_point(colour = COL_GDP, alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = COL_GREY) +
  labs(title = "KOF vs SECO SA (QoQ %)",
       x = "SECO SA QoQ %", y = "KOF sports-corrected QoQ %")

fig_gdp <- (p_gdp_qoq | p_gdp_yoy | p_gdp_scatter) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig_gdp <- fig_gdp +
  plot_annotation(
    title    = "Swiss GDP — target variable overview",
    subtitle = paste0("KOF sports-event corrected (dark blue) vs SECO seasonally adjusted (light blue)",
                      " | Source: SECO / KOF | Quarterly"),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "gray40")
    )
  )

# ggsave(file.path(PATH_FIG, "01a_gdp.png"),
#        fig_gdp, width = 14, height = 5, dpi = 150)
print(fig_gdp)

# ============================================================================
# 2.  CPI — four-panel: level, MoM, QoQ, YoY
#     Primary (red)    = headline CPI
#     Overlay (salmon) = core CPI (excl. food, energy, seasonal ≈ SA proxy)
# ============================================================================
message("Plotting CPI...")

# No level series available — growth rates only
p_cpi_mom <- panel(
  s1 = cpi$pct_1m,      col1 = COL_CPI,    label1 = "Headline MoM %",
  s2 = cpi$core_pct_1m, col2 = COL_CPI_SA, label2 = "Core MoM %",
  title = "Month-on-month %", y_lab = "%"
)

p_cpi_qoq <- panel(
  s1 = cpi$pct_3m,      col1 = COL_CPI,    label1 = "Headline QoQ %",
  s2 = cpi$core_pct_3m, col2 = COL_CPI_SA, label2 = "Core QoQ %",
  title = "Quarter-on-quarter %", y_lab = "%"
)

p_cpi_yoy <- panel(
  s1 = cpi$pct_1y,      col1 = COL_CPI,    label1 = "Headline YoY %",
  s2 = cpi$core_pct_1y, col2 = COL_CPI_SA, label2 = "Core YoY %",
  title = "Year-on-year %", y_lab = "%"
)

fig_cpi <- (p_cpi_mom | p_cpi_qoq | p_cpi_yoy) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig_cpi <- fig_cpi +
  plot_annotation(
    title    = "Swiss CPI — target variable overview",
    subtitle = "Headline (red) vs Core CPI (salmon) | Source: Swiss Federal Statistical Office | Monthly",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "gray40")
    )
  )

# ggsave(file.path(PATH_FIG, "01b_cpi.png"),
#        fig_cpi, width = 12, height = 8, dpi = 150)
print(fig_cpi)

# ============================================================================
# 3.  Bond yield — four-panel: level, MoM change, QoQ change, YoY change
#     No SA alternative (yield is already a rate; no seasonal adjustment needed)
# ============================================================================
message("Plotting bond yield...")

p_bond_lvl <- ggplot(bond$lvl %>% filter(date >= "2000-01-01"),
                     aes(x = date, y = value)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = COL_GREY, linewidth = 0.35) +
  geom_line(colour = COL_BOND, linewidth = 0.75) +
  labs(title = "Level (yield, %)", x = NULL, y = "%") +
  scale_x_date(date_breaks = "4 years", date_labels = "%Y")

p_bond_mom <- panel(
  s1 = bond$dif_1m, col1 = COL_BOND, label1 = "MoM change (pp)",
  title = "Month-on-month change", y_lab = "pp"
)

# QoQ = sum of 3 monthly changes within each quarter
bond_qoq <- bond$dif_3m %>%
  mutate(date = as.Date(as.yearqtr(date))) %>%
  group_by(date) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

p_bond_qoq <- panel(
  s1 = bond_qoq, col1 = COL_BOND, label1 = "QoQ change (pp, summed)",
  title = "Quarter-on-quarter change", y_lab = "pp"
)

p_bond_yoy <- panel(
  s1 = bond$dif_1y, col1 = COL_BOND, label1 = "YoY change (pp)",
  title = "Year-on-year change", y_lab = "pp"
)

fig_bond <- (p_bond_lvl | p_bond_mom) / (p_bond_qoq | p_bond_yoy) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig_bond <- fig_bond +
  plot_annotation(
    title    = "Swiss 10Y government bond yield — target variable overview",
    subtitle = "End-of-period yield | Source: Swiss National Bank | Monthly",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "gray40")
    )
  )

# ggsave(file.path(PATH_FIG, "01c_bond.png"),
#        fig_bond, width = 12, height = 8, dpi = 150)
print(fig_bond)

# ============================================================================
# 4.  Summary statistics
# ============================================================================
message("Computing summary statistics...")

stats_table <- bind_rows(
  # GDP (dif_* series excluded — stored in absolute CHF, not growth rates)
  sumstat(gdp$pct_3m,    "GDP QoQ % — KOF sports-corrected"),
  sumstat(gdp$pct_1y,    "GDP YoY % — KOF sports-corrected"),
  sumstat(gdp$sa_pct_3m, "GDP QoQ % — SECO SA"),
  sumstat(gdp$sa_pct_1y, "GDP YoY % — SECO SA"),
  # CPI
  sumstat(cpi$pct_1m,      "CPI MoM % — headline"),
  sumstat(cpi$pct_3m,      "CPI QoQ % — headline"),
  sumstat(cpi$pct_1y,      "CPI YoY % — headline"),
  sumstat(cpi$core_pct_1m, "CPI MoM % — core"),
  sumstat(cpi$core_pct_3m, "CPI QoQ % — core"),
  sumstat(cpi$core_pct_1y, "CPI YoY % — core"),
  # Bond
  sumstat(bond$lvl,    "Bond yield level (%)"),
  sumstat(bond$dif_1m, "Bond yield MoM change (pp)"),
  sumstat(bond$dif_1y, "Bond yield YoY change (pp)"),
  sumstat(bond_qoq,    "Bond yield QoQ change (pp, summed)")
)

print(stats_table, n = Inf, width = Inf)

# write.csv(stats_table,
#           file.path(PATH_TAB, "01_summary_stats.csv"),
#           row.names = FALSE)

message("\nNote: ADF p truncated at 0.01 by tseries package.")

# ============================================================================
# 5.  Formal outlier detection — 3×IQR rule
# ============================================================================
message("Detecting outliers...")

detect_outliers <- function(df, label, k = 3) {
  x       <- df$value
  Q1      <- quantile(x, 0.25, na.rm = TRUE)
  Q3      <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower   <- Q1 - k * IQR_val
  upper   <- Q3 + k * IQR_val
  df %>%
    filter(value < lower | value > upper) %>%
    mutate(series      = label,
           lower_bound = round(lower, 3),
           upper_bound = round(upper, 3),
           value       = round(value, 3))
}

outliers <- bind_rows(
  detect_outliers(gdp$pct_3m,    "GDP QoQ % — KOF"),
  detect_outliers(gdp$pct_1y,    "GDP YoY % — KOF"),
  detect_outliers(gdp$sa_pct_3m, "GDP QoQ % — SECO SA"),
  detect_outliers(gdp$sa_pct_1y, "GDP YoY % — SECO SA"),
  detect_outliers(cpi$pct_1m,    "CPI MoM % — headline"),
  detect_outliers(cpi$pct_3m,    "CPI QoQ % — headline"),
  detect_outliers(cpi$pct_1y,    "CPI YoY % — headline"),
  detect_outliers(cpi$core_pct_1m, "CPI MoM % — core"),
  detect_outliers(cpi$core_pct_3m, "CPI QoQ % — core"),
  detect_outliers(cpi$core_pct_1y, "CPI YoY % — core"),
  detect_outliers(bond$lvl,      "Bond yield level (%)"),
  detect_outliers(bond$dif_1m,   "Bond yield MoM change (pp)"),
  detect_outliers(bond_qoq,      "Bond yield QoQ change (pp)")
)

cat("\n=== Outliers identified (3×IQR rule) ===\n")
print(outliers %>% select(series, date, value, lower_bound, upper_bound),
      n = Inf, width = Inf)

# write.csv(outliers %>% select(series, date, value, lower_bound, upper_bound),
#           file.path(PATH_TAB, "01_outliers.csv"),
#           row.names = FALSE)

# ── NA check ─────────────────────────────────────────────────────────────────
cat("\n=== NA counts across all series ===\n")
na_found <- FALSE
for (grp_name in c("gdp", "cpi", "bond")) {
  grp <- get(grp_name)
  for (nm in names(grp)) {
    n_na <- sum(is.na(grp[[nm]]$value))
    if (n_na > 0) {
      cat(sprintf("  %s$%s: %d NAs\n", grp_name, nm, n_na))
      na_found <- TRUE
    }
  }
}
if (!na_found) cat("  No NAs found in any series.\n")

message("Exploration complete. Figures → output/figures/ | Tables → output/tables/")

# ============================================================================
# 6.  Granger causality tests
#     H0: X does NOT Granger-cause Y (i.e. lags of X don't help predict Y)
#     Reject H0 (p < 0.05) → X Granger-causes Y
#     Using p=4 lags (1 year) as standard choice; VAR must be estimated first
# ============================================================================
message("Running Granger causality tests...")

# Build the quarterly VAR dataset with QoQ series (our actual target variables)
granger_data <- gdp$pct_3m %>% rename(gdp = value) %>%
  inner_join(cpi$q_pct_3m  %>% rename(cpi  = value), by = "date") %>%
  inner_join(bond_q_dif     %>% rename(bond = bond_dif), by = "date") %>%
  arrange(date)

# Helper: run grangertest and extract p-value cleanly
granger_p <- function(y, x, order = 4) {
  df <- data.frame(y = y, x = x)
  tryCatch(
    lmtest::grangertest(y ~ x, order = order, data = df)$`Pr(>F)`[2],
    error = function(e) NA_real_
  )
}

# All six directional pairs
granger_results <- tribble(
  ~cause, ~effect, ~p_value,
  "GDP",  "CPI",  granger_p(granger_data$cpi,  granger_data$gdp),
  "GDP",  "Bond", granger_p(granger_data$bond, granger_data$gdp),
  "CPI",  "GDP",  granger_p(granger_data$gdp,  granger_data$cpi),
  "CPI",  "Bond", granger_p(granger_data$bond, granger_data$cpi),
  "Bond", "GDP",  granger_p(granger_data$gdp,  granger_data$bond),
  "Bond", "CPI",  granger_p(granger_data$cpi,  granger_data$bond)
) %>%
  mutate(
    p_value    = round(p_value, 4),
    significant = case_when(
      p_value <= 0.01 ~ "*** (1%)",
      p_value <= 0.05 ~ "**  (5%)",
      p_value <= 0.10 ~ "*   (10%)",
      TRUE            ~ "not significant"
    ),
    interpretation = paste0(cause, " Granger-causes ", effect, ": ", significant)
  )

cat("\n=== Granger Causality Tests (p=4 lags) ===\n")
cat("H0: row variable does NOT Granger-cause column variable\n\n")
print(granger_results %>% select(cause, effect, p_value, significant), n = Inf)
cat("\nInterpretations:\n")
cat(paste0("  ", granger_results$interpretation), sep = "\n")
cat("\nNote: Granger causality ≠ structural causality.\n")
cat("      Results depend on lag length — robustness to p=1,2,8 worth checking.\n")