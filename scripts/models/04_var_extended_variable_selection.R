# ============================================================================
# scripts/models/04_var_extended_selection.R
#
# PURPOSE:
#   Identify and evaluate candidate predictor variables to add to the
#   baseline 3-variable VAR (GDP, CPI, bond yield change).
#   Selection evaluates predictive content for ALL THREE target variables,
#   not just GDP — a good predictor for CPI or bond yields is equally valuable.
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/01_var_baseline.R")   — var_input, samples
#
# CANDIDATES (all confirmed in JSON):
#   NOTE: Oil price (oilbren) is NOT evaluated here as an endogenous candidate.
#   Per Prof. Rathke: oil price should be included as EXOGENOUS only (VARX).
#   It will be passed via exogen= alongside COVID dummies in 05_var_extended.R.
#
#   FX — export competitiveness, imported inflation
#     swxsfec_          : CHF/EUR (QoQ %) — dominant trade relationship
#
#   Euro Area — confirmed by prof as important (circled in meeting notes)
#     ekeusesig         : EA Economic Sentiment Indicator (level)
#     bdiptot_g         : German IP (QoQ %) — collinear with EA ESI, for comparison
#     empmim_hq         : EA PMI manufacturing (level) — collinear, for comparison
#
#   World demand (Welt Nachfrage) — global cycle proxy
#     mswrld_d_         : MSCI World equity index (QoQ %) — best available proxy
#                         No trade-weighted world GDP series in dataset
#
#   Labour market (Arbeitsmarkt) — flagged by prof as important
#     swun_p_totq       : Swiss unemployment rate, SA (level) — stationary ✓
#                         Also a candidate as 4th target variable
#
#   Price measures / deflators (Preisvasse, Deflatoren)
#     swimpprce         : Import prices (QoQ %) — FX/oil transmission into prices
#     swproprce         : Producer Price Index (QoQ %) — upstream price pressure
#     swppingdf         : Producer & Import PPI, domestic (QoQ %) — combined measure
#
#   Swiss leading indicators — domestic forward-looking
#     ch_kof_barometer  : KOF economic barometer (level)
#     swpurchsq         : Swiss PMI total (level)
#     swobs085q         : OECD composite leading indicator CH (level)
#
# ECONOMIC RATIONALE PER TARGET:
#   GDP  : EA/DE indicators, Swiss leading indicators, FX (competitiveness)
#   CPI  : Oil/commodities (input costs), FX (imported inflation),
#          EA ESI (demand-pull from trading partners)
#   Bond : EA indicators (SNB follows ECB), oil (inflation expectations),
#          US ISM (global risk appetite), FX (safe-haven flows)
#
# TRANSFORMATION CHOICE:
#   Survey / index levels (PMI, ESI, barometer): _lvl — already stationary
#   Growth rates / prices (IP, FX, oil): _pct_3m — QoQ % change
#   swcnfbusq: _dif_3m (no _pct available)
#   All aggregated to quarterly (mean of 3 months within quarter)
#
# SELECTION CRITERIA (applied to all three targets):
#   1. Stationarity (ADF test on quarterly aggregated series)
#   2. Correlation with each target variable (contemporaneous + lag 1)
#   3. Granger causality: X → GDP, X → CPI, X → bond (p < 0.10)
#   4. Data coverage (overlap with VAR sample Q2 2000 – Q4 2025)
#   5. Economic rationale (no purely statistical selection)
#
# PLANNED VARIABLE SETS (for 05_var_extended.R):
#   Oil price (oilbren) included as EXOGENOUS in ALL models via exogen=.
#   Variable sets below refer to ENDOGENOUS variables only.
#   All sets estimated at p=1,2,4 × full/post08/post15 samples.
#
#   4-variable models (baseline + 1):
#   set_fx    : baseline + swxsfec_            (FX channel)
#   set_ea    : baseline + ekeusesig           (EA demand channel)
#   set_world : baseline + mswrld_d_           (global demand channel)
#   set_unemp : baseline + swun_p_totq         (labour market channel)
#   set_imp   : baseline + swimpprce           (import price / deflator channel)
#   set_sw    : baseline + ch_kof_barometer    (Swiss leading indicator)
#
#   5-variable models (baseline + 2):
#   set_small : baseline + swxsfec_ + ekeusesig
#                → FX + EA demand (core external channels)
#   set_sw2   : baseline + swxsfec_ + ch_kof_barometer
#                → FX + Swiss leading
#
#   6-variable models (baseline + 3):
#   set_medium: baseline + swxsfec_ + ekeusesig + ch_kof_barometer
#                → external + domestic leading
#   set_price : baseline + swxsfec_ + ekeusesig + swimpprce
#                → external + price channel
#
#   7-variable models (baseline + 4):
#   set_large : baseline + swxsfec_ + ekeusesig + ch_kof_barometer + swun_p_totq
#                → external + domestic + labour market
#
#   8-variable models (baseline + 5):
#   set_full  : baseline + swxsfec_ + ekeusesig + mswrld_d_ +
#               ch_kof_barometer + swun_p_totq
#                → all main channels
#
#   (revise after reviewing selection results below)
#
# FINDINGS (update after re-running with new candidates):
#   Previous run (13 candidates): all stationary, Brent caused all 3 targets
#   New candidates added: MSCI World, Swiss unemployment, import prices,
#                         producer prices, Producer & Import PPI
#   Oil price moved to exogenous — not evaluated here
#   CHF/USD, IMF commodity, US ISM, KOF business survey dropped
#   EA collinearity: keep EA ESI as primary, DE IP and EA PMI for comparison
#   [fill in correlation and Granger results after running]
#
# Outputs (commented — uncomment to save):
#   output/figures/04_candidate_overview.png
#   output/tables/04_selection_summary.csv
# ============================================================================

library(dplyr)

# ── 0. Target variables ───────────────────────────────────────────────────────
targets <- list(
  gdp_g    = list(label = "GDP growth (QoQ %)",       data = var_input %>% select(date, value = gdp_g)),
  cpi_g    = list(label = "CPI inflation (QoQ %)",    data = var_input %>% select(date, value = cpi_g)),
  bond_dif = list(label = "Bond yield change (QoQ pp)", data = var_input %>% select(date, value = bond_dif))
)

# ── 1. Define candidates ──────────────────────────────────────────────────────
# NOTE: Oil price (oilbren) is NOT included here as a candidate for
# endogenous selection — per Prof. Rathke's instruction, oil price is
# always included as EXOGENOUS only (passed via exogen= in VARX).
# It will be added to all models automatically in 05_var_extended.R.

candidates <- tribble(
  ~key,               ~transform, ~label,                             ~group,
  # FX — export competitiveness, imported inflation
  "swxsfec_",         "pct_3m",   "CHF/EUR (QoQ %)",                  "FX",
  # Euro Area — Switzerland's main export market (confirmed by prof)
  "ekeusesig",        "lvl",      "EA Economic Sentiment (level)",    "Euro Area",
  # World demand proxy — global cycle
  "mswrld_d_",        "pct_3m",   "MSCI World (QoQ %)",               "World demand",
  # Labour market — domestic activity, potential 4th target variable
  "swun_p_totq",      "lvl",      "Swiss unemployment rate (level)",  "Labour market",
  # Price measures / deflators
  "swimpprce",        "pct_3m",   "Import prices (QoQ %)",            "Prices/Deflators",
  "swproprce",        "pct_3m",   "Producer Price Index (QoQ %)",     "Prices/Deflators",
  "swppingdf",        "pct_3m",   "Producer & Import PPI (QoQ %)",    "Prices/Deflators",
  # Swiss leading indicators — domestic forward-looking
  "ch_kof_barometer", "lvl",      "KOF barometer (level)",            "Swiss leading",
  "swpurchsq",        "lvl",      "Swiss PMI total (level)",          "Swiss leading",
  "swobs085q",        "lvl",      "OECD CLI Switzerland (level)",     "Swiss leading",
  # Euro Area alternatives (collinear — keep for comparison)
  "bdiptot_g",        "pct_3m",   "German IP (QoQ %)",                "Euro Area",
  "empmim_hq",        "lvl",      "EA PMI manufacturing (level)",     "Euro Area"
)

cat("Candidates:", nrow(candidates), "| Targets: GDP, CPI, Bond yield\n")
cat("Note: Brent oil excluded from endogenous candidates — treated as exogenous in all VARX models\n\n")

# ── 2. Extract and aggregate to quarterly ─────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 1: EXTRACT AND AGGREGATE TO QUARTERLY\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

extract_candidate <- function(key, transform) {
  json_key  <- paste0(key, "_", transform)
  vals      <- nowcast_raw[[json_key]]
  dates_chr <- nowcast_raw$dates[[json_key]]
  if (is.null(vals) || is.null(dates_chr)) {
    warning("Key not found: ", json_key); return(NULL)
  }
  tibble(date = dmy(dates_chr), value = as.numeric(vals)) %>%
    filter(!is.na(value)) %>%
    mutate(date = as.Date(as.yearqtr(date))) %>%
    group_by(date) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
}

cand_data <- list()
for (i in seq_len(nrow(candidates))) {
  row <- candidates[i, ]
  s   <- extract_candidate(row$key, row$transform)
  if (!is.null(s)) {
    cand_data[[row$key]] <- s
    cat(sprintf("  %-22s %-8s %3d quarterly obs  %s – %s\n",
                row$key, row$transform, nrow(s),
                format(min(s$date)), format(max(s$date))))
  }
}

# ── 3. Stationarity ───────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 2: STATIONARITY (ADF TESTS)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

adf_res <- candidates %>%
  rowwise() %>%
  mutate(
    adf_p = tryCatch(
      adf.test(na.omit(cand_data[[key]]$value))$p.value,
      error = function(e) NA_real_
    ),
    stationary = case_when(
      is.na(adf_p)  ~ "n/a",
      adf_p <= 0.05 ~ "✓",
      adf_p <= 0.10 ~ "~",
      TRUE          ~ "✗"
    )
  ) %>% ungroup()

adf_res %>%
  select(group, label, transform, adf_p, stationary) %>%
  mutate(adf_p = round(adf_p, 3)) %>%
  arrange(group) %>%
  print(n = Inf, width = Inf)

# ── 4. Correlations with all three targets ────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 3: CORRELATIONS WITH ALL THREE TARGET VARIABLES\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Contemporaneous (lag=0) and 1-quarter lead (lag=1: candidate leads target)\n\n")

cor_res <- candidates %>%
  rowwise() %>%
  mutate(
    # GDP
    cor_gdp_0  = {
      m <- inner_join(targets$gdp_g$data,    cand_data[[key]], by = "date")
      round(cor(m$value.x, m$value.y, use = "complete.obs"), 3)
    },
    cor_gdp_1  = {
      s <- cand_data[[key]] %>% mutate(date = date + months(3))
      m <- inner_join(targets$gdp_g$data, s, by = "date")
      round(cor(m$value.x, m$value.y, use = "complete.obs"), 3)
    },
    # CPI
    cor_cpi_0  = {
      m <- inner_join(targets$cpi_g$data,    cand_data[[key]], by = "date")
      round(cor(m$value.x, m$value.y, use = "complete.obs"), 3)
    },
    cor_cpi_1  = {
      s <- cand_data[[key]] %>% mutate(date = date + months(3))
      m <- inner_join(targets$cpi_g$data, s, by = "date")
      round(cor(m$value.x, m$value.y, use = "complete.obs"), 3)
    },
    # Bond
    cor_bond_0 = {
      m <- inner_join(targets$bond_dif$data, cand_data[[key]], by = "date")
      round(cor(m$value.x, m$value.y, use = "complete.obs"), 3)
    },
    cor_bond_1 = {
      s <- cand_data[[key]] %>% mutate(date = date + months(3))
      m <- inner_join(targets$bond_dif$data, s, by = "date")
      round(cor(m$value.x, m$value.y, use = "complete.obs"), 3)
    }
  ) %>% ungroup()

cat("── Contemporaneous correlations ────────────────────────────────────\n")
cor_res %>%
  select(group, label, cor_gdp_0, cor_cpi_0, cor_bond_0) %>%
  arrange(desc(pmax(abs(cor_gdp_0), abs(cor_cpi_0), abs(cor_bond_0)))) %>%
  print(n = Inf, width = Inf)

cat("\n── 1-quarter lead correlations (candidate leads target) ─────────────\n")
cor_res %>%
  select(group, label, cor_gdp_1, cor_cpi_1, cor_bond_1) %>%
  arrange(desc(pmax(abs(cor_gdp_1), abs(cor_cpi_1), abs(cor_bond_1)))) %>%
  print(n = Inf, width = Inf)

# ── 5. Granger causality: X → each target ────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 4: GRANGER CAUSALITY X → GDP / CPI / BOND\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Min p-value across lags p=1,2,4 | *** <1%  ** <5%  * <10%\n\n")

granger_p <- function(target_data, cand_series, order) {
  m <- inner_join(target_data, cand_series, by = "date") %>%
    filter(!is.na(value.x), !is.na(value.y))
  if (nrow(m) < order * 2 + 10) return(NA_real_)
  tryCatch(
    lmtest::grangertest(value.x ~ value.y, order = order, data = m)$`Pr(>F)`[2],
    error = function(e) NA_real_
  )
}

granger_res <- candidates %>%
  rowwise() %>%
  mutate(
    # GDP
    p_gdp  = pmin(
      granger_p(targets$gdp_g$data,    cand_data[[key]], 1),
      granger_p(targets$gdp_g$data,    cand_data[[key]], 2),
      granger_p(targets$gdp_g$data,    cand_data[[key]], 4),
      na.rm = TRUE),
    # CPI
    p_cpi  = pmin(
      granger_p(targets$cpi_g$data,    cand_data[[key]], 1),
      granger_p(targets$cpi_g$data,    cand_data[[key]], 2),
      granger_p(targets$cpi_g$data,    cand_data[[key]], 4),
      na.rm = TRUE),
    # Bond
    p_bond = pmin(
      granger_p(targets$bond_dif$data, cand_data[[key]], 1),
      granger_p(targets$bond_dif$data, cand_data[[key]], 2),
      granger_p(targets$bond_dif$data, cand_data[[key]], 4),
      na.rm = TRUE)
  ) %>%
  mutate(across(c(p_gdp, p_cpi, p_bond), ~ round(.x, 3))) %>%
  mutate(
    sig_gdp  = case_when(p_gdp  <= 0.01 ~ "***", p_gdp  <= 0.05 ~ "**",
                         p_gdp  <= 0.10 ~ "*",   TRUE ~ "—"),
    sig_cpi  = case_when(p_cpi  <= 0.01 ~ "***", p_cpi  <= 0.05 ~ "**",
                         p_cpi  <= 0.10 ~ "*",   TRUE ~ "—"),
    sig_bond = case_when(p_bond <= 0.01 ~ "***", p_bond <= 0.05 ~ "**",
                         p_bond <= 0.10 ~ "*",   TRUE ~ "—"),
    n_sig    = (sig_gdp != "—") + (sig_cpi != "—") + (sig_bond != "—")
  ) %>% ungroup()

granger_res %>%
  select(group, label, p_gdp, sig_gdp, p_cpi, sig_cpi, p_bond, sig_bond, n_sig) %>%
  arrange(desc(n_sig), p_gdp) %>%
  print(n = Inf, width = Inf)

# ── 6. Overview plots ─────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 5: OVERVIEW PLOTS — CANDIDATE vs ALL THREE TARGETS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

group_colours <- c(
  "FX"               = "#7b3294",
  "Euro Area"        = "#e66101",
  "World demand"     = "#d01c8b",
  "Labour market"    = "#4dac26",
  "Prices/Deflators" = "#8c510a",
  "Swiss leading"    = COL_GDP
)

# One figure per candidate group: each row = one candidate,
# columns = GDP / CPI / Bond overlay
for (grp in unique(candidates$group)) {
  keys_grp   <- candidates %>% filter(group == grp) %>% pull(key)
  labels_grp <- candidates %>% filter(group == grp) %>% pull(label)
  col_grp    <- group_colours[grp]
  
  panels <- list()
  for (i in seq_along(keys_grp)) {
    k  <- keys_grp[i]
    lb <- labels_grp[i]
    s  <- cand_data[[k]]
    if (is.null(s)) next
    
    for (tname in names(targets)) {
      tgt <- targets[[tname]]$data %>% filter(date >= min(s$date))
      cnd <- s %>% filter(date >= min(tgt$date))
      
      # Z-score both for overlay
      tgt_z <- tgt %>%
        mutate(z = (value - mean(value, na.rm=TRUE)) / sd(value, na.rm=TRUE))
      cnd_z <- cnd %>%
        mutate(z = (value - mean(value, na.rm=TRUE)) / sd(value, na.rm=TRUE))
      
      target_col <- c(gdp_g = COL_GDP, cpi_g = COL_CPI, bond_dif = COL_BOND)[tname]
      
      panels[[paste(k, tname, sep="_")]] <- ggplot() +
        geom_hline(yintercept = 0, linetype = "dashed",
                   color = COL_GREY, linewidth = 0.3) +
        geom_line(data = tgt_z, aes(x = date, y = z),
                  colour = target_col, linewidth = 0.6, alpha = 0.6) +
        geom_line(data = cnd_z, aes(x = date, y = z),
                  colour = col_grp, linewidth = 0.8) +
        scale_x_date(date_breaks = "4 years", date_labels = "%Y") +
        labs(
          title    = sprintf("%s", lb),
          subtitle = sprintf("vs %s", targets[[tname]]$label),
          x = NULL, y = "Z-score"
        )
    }
  }
  
  # layout: rows = candidates, cols = GDP / CPI / Bond
  n_cands <- length(keys_grp)
  fig <- wrap_plots(panels, nrow = n_cands, ncol = 3) +
    plot_annotation(
      title    = sprintf("Candidate predictors — %s", grp),
      subtitle = "Z-scored | Coloured = candidate | Faded = target (blue=GDP, red=CPI, green=Bond)",
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "gray40")
      )
    )
  print(fig)
  
  # # Save (uncomment)
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("04_candidates_%s.png",
  #                          gsub("/| ", "_", grp))),
  #        fig, width = 14, height = 4 * n_cands, dpi = 150)
}

# ── 7. Final selection summary ────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 6: SELECTION SUMMARY\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

selection_summary <- adf_res %>%
  select(key, group, label, transform, adf_p, stationary) %>%
  left_join(cor_res    %>% select(key, cor_gdp_0, cor_gdp_1,
                                  cor_cpi_0, cor_cpi_1,
                                  cor_bond_0, cor_bond_1), by = "key") %>%
  left_join(granger_res %>% select(key, sig_gdp, sig_cpi, sig_bond, n_sig), by = "key") %>%
  mutate(
    max_cor = pmax(abs(cor_gdp_0), abs(cor_cpi_0), abs(cor_bond_0),
                   abs(cor_gdp_1), abs(cor_cpi_1), abs(cor_bond_1), na.rm = TRUE),
    recommendation = case_when(
      adf_p > 0.15                       ~ "exclude — non-stationary",
      max_cor < 0.10 & n_sig == 0        ~ "exclude — no predictive content",
      n_sig >= 2 | max_cor >= 0.30       ~ "include — strong candidate",
      TRUE                               ~ "consider — moderate evidence"
    )
  ) %>%
  arrange(desc(n_sig), desc(max_cor))

cat("Summary (sorted by number of significant Granger targets):\n\n")
selection_summary %>%
  select(group, label, stationary, cor_gdp_0, cor_cpi_0, cor_bond_0,
         sig_gdp, sig_cpi, sig_bond, recommendation) %>%
  print(n = Inf, width = Inf)

cat("\n── Strong candidates (include) ──────────────────────────────────────\n")
selection_summary %>%
  filter(grepl("include", recommendation)) %>%
  select(group, label, key) %>%
  print(n = Inf)

cat("\n── Moderate candidates (consider) ───────────────────────────────────\n")
selection_summary %>%
  filter(grepl("consider", recommendation)) %>%
  select(group, label, key) %>%
  print(n = Inf)

# # Save (uncomment)
# write.csv(selection_summary,
#           file.path(here("output","tables"), "04_selection_summary.csv"),
#           row.names = FALSE)

cat("\nSelection complete. Update planned variable sets in header and\n")
cat("proceed to 05_var_extended.R.\n")