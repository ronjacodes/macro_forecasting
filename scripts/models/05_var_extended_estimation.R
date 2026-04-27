# ============================================================================
# scripts/models/05_var_extended_estimation.R
#
# PURPOSE:
#   Estimate extended VARX models adding candidate predictors one at a time
#   and in combinations to the 3-variable baseline (GDP, CPI, bond yield).
#   Oil price (Brent) is included as EXOGENOUS in ALL models, alongside
#   the COVID dummies — per Prof. Rathke's instruction.
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/01_var_baseline.R")   — var_input, samples,
#                                                   exog_full, covid_cols,
#                                                   subset_sample(), make_ts(),
#                                                   run_diagnostics()
#   source("scripts/models/04_var_extended_selection.R") — cand_data
#
# ENDOGENOUS VARIABLES (3 baseline + up to 6 additional):
#   Always: gdp_g, cpi_g, bond_dif
#   Added:  swxsfec_ (CHF/EUR), ekeusesig (EA ESI), mswrld_d_ (MSCI World),
#           swimpprce (import prices), ch_kof_barometer (KOF barometer),
#           swun_p_totq (unemployment)
#
# EXOGENOUS VARIABLES:
#   Always: 8 COVID dummies (carried over from baseline)
#   When oil is included in a set: oilbren (Brent crude, QoQ %) added to exogen=
#   Per prof: when oil IS included, it must be exogenous — not endogenous
#   Sets without oil: use only COVID dummies as exogen (same as baseline)
#
# MODEL GRID:
#   Variable sets × lag orders × sample periods
#
#   4-variable sets (baseline + 1 endogenous):
#     set_fx    : + swxsfec_            FX / imported inflation
#     set_ea    : + ekeusesig           EA demand
#     set_world : + mswrld_d_           global demand (Welt Nachfrage)
#     set_imp   : + swimpprce           import price / deflator channel
#     set_sw    : + ch_kof_barometer    Swiss domestic leading
#     set_unemp : + swun_p_totq         labour market (Arbeitsmarkt)
#     set_oil   : [oil exog only]       Brent as exogenous — no new endogenous
#
#   5-variable sets (baseline + 2):
#     set_ext   : + swxsfec_ + ekeusesig              FX + EA
#     set_dom   : + ch_kof_barometer + swimpprce       domestic + prices
#     set_oil_fx: + swxsfec_  [oil exog]               FX + oil exogenous
#
#   6-variable sets (baseline + 3):
#     set_medium   : + swxsfec_ + ekeusesig + ch_kof_barometer
#     set_price    : + swxsfec_ + ekeusesig + swimpprce
#     set_oil_ext  : + swxsfec_ + ekeusesig  [oil exog]
#
#   7-variable sets (baseline + 4):
#     set_large    : + swxsfec_ + ekeusesig + mswrld_d_ + ch_kof_barometer
#     set_oil_med  : + swxsfec_ + ekeusesig + ch_kof_barometer  [oil exog]
#
#   8-variable sets (baseline + 5):
#     set_full     : + swxsfec_ + ekeusesig + mswrld_d_ + swimpprce + ch_kof_barometer
#     set_oil_large: + swxsfec_ + ekeusesig + mswrld_d_ + ch_kof_barometer  [oil exog]
#
#   9-variable sets (baseline + 6):
#     set_max      : all 6 additional variables
#     set_oil_full : + swxsfec_ + ekeusesig + mswrld_d_ + swimpprce + ch_kof_barometer  [oil exog]
#
#   Lag orders : p = 1, 2, 4
#   Sample periods: full, post08, post15
#   → up to ~18 sets × 3 lags × 3 samples = ~162 models
#      (many skipped for df constraints)
#
# STRUCTURE:
#   Section 0 — Prepare exogenous matrix (oil + COVID)
#   Section 1 — Extract and align extended variables
#   Section 2 — Define variable sets
#   Section 3 — Estimate all models
#   Section 4 — Summary table and Q1 2026 forecasts
#   Section 5 — Forecast plots for best models per set
#
# KEY FINDINGS:
#
#   ── Model size vs. fit ───────────────────────────────────────────────────
#   Adding variables improves in-sample fit but gains diminish quickly:
#     3-var baseline:  mean σ(GDP) = 0.389
#     4-var (+1 var):  mean σ(GDP) = 0.405  (some 4-var beat baseline)
#     5-var (+2 var):  mean σ(GDP) = 0.400
#     6-var (+3 var):  mean σ(GDP) = 0.440  (larger models hurt — df cost)
#     7-8-var:         mean σ(GDP) = 0.425–0.426
#   → Post-GFC / post-2015 samples + p=1 dominate the top 10 rankings
#   → Larger models mostly only estimable on full sample at p=1 (df constraint)
#
#   ── Best models per size ─────────────────────────────────────────────────
#   Best 3-var: set_oil_post15_p1      σ(GDP)=0.355  Q1 GDP=+0.43%
#   Best 4-var: set_sw_post08_p1       σ(GDP)=0.336  Q1 GDP=+0.57%
#               set_world_post15_p1    σ(GDP)=0.341  Q1 GDP=+0.50%
#               set_ea_post15_p1       σ(GDP)=0.344  Q1 GDP=+0.37%
#   Best 5-var: set_dom_post08_p1      σ(GDP)=0.335  Q1 GDP=+0.61%
#   Best 7-var: set_oil_lrg_full_p1    σ(GDP)=0.416  Q1 GDP=+0.54%
#   Best 8-var: set_oil_full_full_p1   σ(GDP)=0.417  Q1 GDP=+0.56%
#
#   ── Effect of adding oil as exogenous ────────────────────────────────────
#   Comparing matched pairs (same endogenous, with vs without oil exog):
#   Oil exog clearly improves fit at larger model sizes (7-8-var)
#   set_oil_lrg vs set_large: σ(GDP) 0.416 vs 0.435 — oil helps
#   set_oil_full vs set_full:  σ(GDP) 0.417 vs 0.434 — oil helps
#   At small sizes (3-4-var) the difference is minimal
#
#   ── Which variable adds most value? ──────────────────────────────────────
#   KOF barometer (set_sw):  best single-variable addition for GDP
#   MSCI World (set_world):  second best, helps CPI too
#   EA ESI (set_ea):         strong for GDP, useful for CPI at longer horizon
#   KOF + import prices (set_dom): best 5-var combination
#   CHF/EUR: useful for CPI channel, modest GDP improvement
#   Unemployment: improves GDP in-sample but persistent serial correlation issues
#
#   ── CPI and bond yield forecasts Q1 2026 ─────────────────────────────────
#   CPI Q1 2026: all models forecast near 0% QoQ (−0.10% to +0.10%)
#     → Consistent with continued Swiss disinflation, near SNB target
#   Bond Q1 2026: range −0.40pp to +0.07pp
#     → Larger models (with MSCI World) forecast small positive yield move
#     → Smaller models forecast slight yield decline
#
#   ── Q1 2026 GDP range ────────────────────────────────────────────────────
#   All stable models: +0.35% to +0.71% QoQ
#   Models with KOF barometer forecast higher (+0.50–0.61%)
#   Models with EA ESI alone forecast lower (+0.37–0.44%)
#   All above Nowcasting Lab benchmark (+0.30%) — consistent with baseline
#
# Outputs (commented — uncomment to save):
#   output/tables/05_extended_summary.csv
#   output/figures/05_forecast_<set>.png
# ============================================================================

library(dplyr)

# ── 0. Prepare exogenous matrices ────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 0: PREPARE EXOGENOUS VARIABLES\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Two exog matrices:\n")
cat("  exog_base    : 8 COVID dummies only (for sets WITHOUT oil)\n")
cat("  exog_with_oil: 8 COVID dummies + Brent QoQ % (for sets WITH oil)\n\n")

# exog_base is already available as exog_full from 01_var_baseline.R
exog_base <- exog_full

# Extract Brent crude QoQ % and aggregate to quarterly
oil_monthly <- tibble(
  date  = dmy(nowcast_raw$dates[["oilbren_pct_3m"]]),
  oil_g = as.numeric(nowcast_raw[["oilbren_pct_3m"]])
) %>%
  filter(!is.na(oil_g)) %>%
  mutate(date = as.Date(as.yearqtr(date))) %>%
  group_by(date) %>%
  summarise(oil_g = mean(oil_g, na.rm = TRUE), .groups = "drop")

cat(sprintf("Brent oil QoQ %%: %d quarterly obs  %s – %s\n",
            nrow(oil_monthly),
            format(min(oil_monthly$date)),
            format(max(oil_monthly$date))))

oil_aligned <- var_input %>%
  select(date) %>%
  left_join(oil_monthly, by = "date") %>%
  pull(oil_g)

n_missing_oil <- sum(is.na(oil_aligned))
if (n_missing_oil > 0) {
  cat(sprintf("  WARNING: %d missing oil values — filling with 0\n", n_missing_oil))
  oil_aligned[is.na(oil_aligned)] <- 0
}

# Extended exog: COVID dummies + oil (used only for oil-inclusive sets)
exog_with_oil     <- cbind(exog_full, oil_g = oil_aligned)
exog_cols_base    <- covid_cols
exog_cols_with_oil <- c(covid_cols, "oil_g")

cat(sprintf("exog_base:     %d rows × %d cols\n",
            nrow(exog_base), ncol(exog_base)))
cat(sprintf("exog_with_oil: %d rows × %d cols\n\n",
            nrow(exog_with_oil), ncol(exog_with_oil)))

# ── 1. Extract and align extended endogenous variables ────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 1: EXTRACT EXTENDED ENDOGENOUS VARIABLES\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

ext_vars <- tribble(
  ~key,               ~transform, ~col_name,  ~label,
  "swxsfec_",         "pct_3m",   "fx_eur",   "CHF/EUR (QoQ %)",
  "ekeusesig",        "lvl",      "ea_esi",   "EA Economic Sentiment",
  "mswrld_d_",        "pct_3m",   "msci_wld", "MSCI World (QoQ %)",
  "swimpprce",        "pct_3m",   "imp_prc",  "Import prices (QoQ %)",
  "ch_kof_barometer", "lvl",      "kof_bar",  "KOF barometer",
  "swun_p_totq",      "lvl",      "unemp",    "Unemployment rate"
)

ext_data <- var_input %>% select(date)

for (i in seq_len(nrow(ext_vars))) {
  row  <- ext_vars[i, ]
  jkey <- paste0(row$key, "_", row$transform)
  vals <- nowcast_raw[[jkey]]
  dts  <- nowcast_raw$dates[[jkey]]
  
  if (is.null(vals) || is.null(dts)) {
    cat(sprintf("  WARNING: %s not found\n", jkey)); next
  }
  
  s <- tibble(date = dmy(dts), v = as.numeric(vals)) %>%
    filter(!is.na(v)) %>%
    mutate(date = as.Date(as.yearqtr(date))) %>%
    group_by(date) %>%
    summarise(v = mean(v, na.rm = TRUE), .groups = "drop") %>%
    rename(!!row$col_name := v)
  
  ext_data <- left_join(ext_data, s, by = "date")
  cat(sprintf("  %-10s %-8s  %3d obs aligned\n",
              row$col_name, row$transform,
              sum(!is.na(ext_data[[row$col_name]]))))
}

cat(sprintf("\next_data: %d rows × %d cols\n\n",
            nrow(ext_data), ncol(ext_data)))

# ── 2. Define variable sets ───────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 2: VARIABLE SETS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

BASE_COLS <- c("gdp_g", "cpi_g", "bond_dif")

# Each set has: id, label, extra (endogenous), oil (logical — use oil as exog?)
# n_endo = 3 (baseline) + length(extra)
var_sets <- list(
  # ── No oil ────────────────────────────────────────────────────────────────
  # 4-variable: baseline + 1
  list(id = "set_fx",      oil = FALSE, n_endo = 4,
       label = "Baseline + CHF/EUR",
       extra = "fx_eur"),
  list(id = "set_ea",      oil = FALSE, n_endo = 4,
       label = "Baseline + EA ESI",
       extra = "ea_esi"),
  list(id = "set_world",   oil = FALSE, n_endo = 4,
       label = "Baseline + MSCI World",
       extra = "msci_wld"),
  list(id = "set_imp",     oil = FALSE, n_endo = 4,
       label = "Baseline + Import prices",
       extra = "imp_prc"),
  list(id = "set_sw",      oil = FALSE, n_endo = 4,
       label = "Baseline + KOF barometer",
       extra = "kof_bar"),
  list(id = "set_unemp",   oil = FALSE, n_endo = 4,
       label = "Baseline + Unemployment",
       extra = "unemp"),
  # 5-variable: baseline + 2
  list(id = "set_ext",     oil = FALSE, n_endo = 5,
       label = "Baseline + CHF/EUR + EA ESI",
       extra = c("fx_eur", "ea_esi")),
  list(id = "set_dom",     oil = FALSE, n_endo = 5,
       label = "Baseline + KOF barometer + Import prices",
       extra = c("kof_bar", "imp_prc")),
  # 6-variable: baseline + 3
  list(id = "set_medium",  oil = FALSE, n_endo = 6,
       label = "Baseline + CHF/EUR + EA ESI + KOF barometer",
       extra = c("fx_eur", "ea_esi", "kof_bar")),
  list(id = "set_price",   oil = FALSE, n_endo = 6,
       label = "Baseline + CHF/EUR + EA ESI + Import prices",
       extra = c("fx_eur", "ea_esi", "imp_prc")),
  # 7-variable: baseline + 4
  list(id = "set_large",   oil = FALSE, n_endo = 7,
       label = "Baseline + CHF/EUR + EA ESI + MSCI World + KOF barometer",
       extra = c("fx_eur", "ea_esi", "msci_wld", "kof_bar")),
  # 8-variable: baseline + 5
  list(id = "set_full",    oil = FALSE, n_endo = 8,
       label = "Baseline + CHF/EUR + EA ESI + MSCI + Import prices + KOF",
       extra = c("fx_eur", "ea_esi", "msci_wld", "imp_prc", "kof_bar")),
  # 9-variable: baseline + 6
  list(id = "set_max",     oil = FALSE, n_endo = 9,
       label = "All 6 additional variables",
       extra = c("fx_eur", "ea_esi", "msci_wld", "imp_prc", "kof_bar", "unemp")),
  
  # ── With oil as exogenous ─────────────────────────────────────────────────
  # Mirrors of selected sets above, with Brent added to exogen=
  list(id = "set_oil",     oil = TRUE,  n_endo = 3,
       label = "Baseline only + oil exog",
       extra = character(0)),
  list(id = "set_oil_fx",  oil = TRUE,  n_endo = 4,
       label = "Baseline + CHF/EUR  [oil exog]",
       extra = "fx_eur"),
  list(id = "set_oil_ext", oil = TRUE,  n_endo = 5,
       label = "Baseline + CHF/EUR + EA ESI  [oil exog]",
       extra = c("fx_eur", "ea_esi")),
  list(id = "set_oil_med", oil = TRUE,  n_endo = 6,
       label = "Baseline + CHF/EUR + EA ESI + KOF  [oil exog]",
       extra = c("fx_eur", "ea_esi", "kof_bar")),
  list(id = "set_oil_lrg", oil = TRUE,  n_endo = 7,
       label = "Baseline + CHF/EUR + EA ESI + MSCI + KOF  [oil exog]",
       extra = c("fx_eur", "ea_esi", "msci_wld", "kof_bar")),
  list(id = "set_oil_full",oil = TRUE,  n_endo = 8,
       label = "Baseline + CHF/EUR + EA ESI + MSCI + imp + KOF  [oil exog]",
       extra = c("fx_eur", "ea_esi", "msci_wld", "imp_prc", "kof_bar"))
)

LAG_ORDERS <- c(1, 2, 4)

cat("Variable sets:\n")
for (vs in var_sets) {
  cat(sprintf("  %-14s [%d-var]  oil=%-5s  extra: %s\n",
              vs$id, vs$n_endo,
              ifelse(vs$oil, "TRUE", "FALSE"),
              if (length(vs$extra) > 0) paste(vs$extra, collapse=", ") else "(none)"))
}
cat(sprintf("\nLag orders: %s | Samples: full, post08, post15\n\n",
            paste(LAG_ORDERS, collapse = ", ")))

# ── 3. Helper functions ───────────────────────────────────────────────────────

# Build ts + aligned data for a given set + sample
make_ext_ts <- function(endo_cols, start_date, ts_start) {
  extra_cols <- endo_cols[!endo_cols %in% BASE_COLS]
  
  dat <- var_input %>%
    select(date, all_of(BASE_COLS)) %>%
    { if (length(extra_cols) > 0)
      left_join(., ext_data %>% select(date, all_of(extra_cols)), by = "date")
      else . } %>%
    filter(date >= start_date) %>%
    arrange(date)
  
  # Drop rows with NAs in any endogenous column
  complete_rows <- complete.cases(dat[, endo_cols])
  dat <- dat[complete_rows, ]
  if (nrow(dat) < 10) return(NULL)
  
  list(
    dat  = dat,
    ts   = ts(dat[, endo_cols], start = ts_start, frequency = 4),
    n    = nrow(dat)
  )
}

# ── 4. Estimate all models ────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 3: ESTIMATE EXTENDED VAR MODELS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

N_AHEAD           <- 8
# Two forecast exog matrices (all zeros — no COVID, flat oil in forecast)
exog_fc_base      <- matrix(0, nrow = N_AHEAD, ncol = ncol(exog_base),
                            dimnames = list(NULL, exog_cols_base))
exog_fc_with_oil  <- matrix(0, nrow = N_AHEAD, ncol = ncol(exog_with_oil),
                            dimnames = list(NULL, exog_cols_with_oil))

ext_results  <- list()
summary_rows <- list()
n_estimated  <- 0
n_skipped    <- 0

for (vs in var_sets) {
  ext_results[[vs$id]] <- list()
  endo_cols <- c(BASE_COLS, vs$extra)
  k         <- length(endo_cols)
  
  cat(sprintf("\n── %s [%d-var]: %s ──\n", vs$id, k, vs$label))
  
  for (sname in names(samples)) {
    s      <- samples[[sname]]
    ext_results[[vs$id]][[sname]] <- list()
    
    ts_obj <- tryCatch(
      make_ext_ts(endo_cols, s$start_date, s$ts_start),
      error = function(e) NULL
    )
    if (is.null(ts_obj)) {
      cat(sprintf("  [%s] SKIPPED — data error\n", sname))
      n_skipped <- n_skipped + length(LAG_ORDERS); next
    }
    
    # Select correct exog matrix for this set
    exog_mat  <- if (vs$oil) exog_with_oil else exog_base
    exog_s    <- exog_mat[var_input$date %in% ts_obj$dat$date, , drop = FALSE]
    exog_fc   <- if (vs$oil) exog_fc_with_oil else exog_fc_base
    n_exog    <- ncol(exog_mat)
    
    for (p in LAG_ORDERS) {
      pkey     <- paste0("p", p)
      id       <- sprintf("%s_%s_p%d", vs$id, sname, p)
      n_params <- k * (k * p + 1) + n_exog
      
      if (ts_obj$n < n_params + 10) {
        cat(sprintf("  [%s p=%d] SKIPPED — n=%d < params+10=%d\n",
                    sname, p, ts_obj$n, n_params + 10))
        n_skipped <- n_skipped + 1
        ext_results[[vs$id]][[sname]][[pkey]] <- NULL
        next
      }
      
      mod <- tryCatch(
        VAR(ts_obj$ts, p = p, type = "const", exogen = exog_s),
        error = function(e) NULL
      )
      if (is.null(mod)) {
        cat(sprintf("  [%s p=%d] SKIPPED — estimation error\n", sname, p))
        n_skipped <- n_skipped + 1
        ext_results[[vs$id]][[sname]][[pkey]] <- NULL
        next
      }
      
      diag <- tryCatch(run_diagnostics(mod),
                       error = function(e) list(
                         serial_p=NA, normal_p=NA, max_root=NA,
                         stable=NA, adj_r2_gdp=NA, sigma_gdp=NA))
      
      fc <- tryCatch(
        predict(mod, n.ahead = N_AHEAD, ci = 0.95, dumvar = exog_fc),
        error = function(e) NULL
      )
      fc_gdp <- if (!is.null(fc)) fc$fcst$gdp_g else
        matrix(NA, N_AHEAD, 3, dimnames = list(NULL, c("fcst","lower","upper")))
      
      ext_results[[vs$id]][[sname]][[pkey]] <- list(
        id = id, set = vs$id, sample = sname, p = p,
        n_endo = k, endo = endo_cols,
        model = mod, diag = diag, fc_gdp = fc_gdp, n_obs = ts_obj$n
      )
      
      # Also extract CPI and bond sigma
      sigma_cpi  <- tryCatch(
        round(summary(mod)$varresult$cpi_g$sigma,    4), error = function(e) NA_real_)
      sigma_bond <- tryCatch(
        round(summary(mod)$varresult$bond_dif$sigma, 4), error = function(e) NA_real_)
      
      # Q1 2026 forecasts for CPI and bond (from same predict() call)
      fc_cpi  <- if (!is.null(fc)) fc$fcst$cpi_g    else matrix(NA, N_AHEAD, 3)
      fc_bond <- if (!is.null(fc)) fc$fcst$bond_dif else matrix(NA, N_AHEAD, 3)
      
      summary_rows[[length(summary_rows) + 1]] <- tibble(
        id         = id,        set        = vs$id,
        label      = vs$label,  n_endo     = k,
        oil_exog   = vs$oil,
        sample     = sname,     p          = p,
        n_obs      = ts_obj$n,
        adj_r2_gdp = round(diag$adj_r2_gdp, 3),
        sigma_gdp  = round(diag$sigma_gdp, 4),
        sigma_cpi  = sigma_cpi,
        sigma_bond = sigma_bond,
        serial_p   = round(diag$serial_p, 3),
        normal_p   = round(diag$normal_p, 3),
        max_root   = round(diag$max_root, 3),
        stable     = diag$stable,
        fc_gdp_q1  = round(fc_gdp[1, "fcst"],  3),
        fc_gdp_lo  = round(fc_gdp[1, "lower"], 3),
        fc_gdp_hi  = round(fc_gdp[1, "upper"], 3),
        fc_cpi_q1  = round(fc_cpi[1, "fcst"],  3),
        fc_bond_q1 = round(fc_bond[1, "fcst"], 3)
      )
      
      n_estimated <- n_estimated + 1
      
      cat(sprintf("  [%s p=%d] n=%d | σ(gdp)=%.4f σ(cpi)=%.4f σ(bond)=%.4f | serial p=%.3f %s | root=%.3f %s | Q1 gdp=%+.3f%% cpi=%+.3f%% bond=%+.3f\n",
                  sname, p, ts_obj$n,
                  diag$sigma_gdp, sigma_cpi, sigma_bond,
                  diag$serial_p, ifelse(!is.na(diag$serial_p) && diag$serial_p > 0.05, "✓","✗"),
                  diag$max_root, ifelse(!is.na(diag$stable) && diag$stable, "✓","✗"),
                  fc_gdp[1, "fcst"], fc_cpi[1, "fcst"], fc_bond[1, "fcst"]))
    }
  }
}

cat(sprintf("\n\nEstimation complete: %d models estimated, %d skipped\n\n",
            n_estimated, n_skipped))

# ── 5. Summary tables ─────────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 4: SUMMARY TABLES\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

summary_tbl <- bind_rows(summary_rows)

# Best model per set
cat("── Best model per variable set (lowest GDP σ, stable + serial OK) ───\n\n")
best_per_set <- summary_tbl %>%
  filter(stable == TRUE, !is.na(serial_p), serial_p > 0.05) %>%
  group_by(set) %>%
  slice_min(sigma_gdp, n = 1) %>%
  ungroup() %>%
  select(set, n_endo, oil_exog, sample, p, n_obs,
         adj_r2_gdp, sigma_gdp, sigma_cpi, sigma_bond,
         serial_p, max_root,
         fc_gdp_q1, fc_cpi_q1, fc_bond_q1,
         fc_gdp_lo, fc_gdp_hi) %>%
  arrange(n_endo, sigma_gdp)

print(best_per_set, n = Inf, width = Inf)

# Average by model size
cat("\n── Average GDP sigma by model size ──────────────────────────────────\n\n")
summary_tbl %>%
  filter(stable == TRUE, !is.na(serial_p), serial_p > 0.05) %>%
  group_by(n_endo) %>%
  summarise(
    n_models        = n(),
    mean_sigma_gdp  = round(mean(sigma_gdp,  na.rm = TRUE), 4),
    mean_sigma_cpi  = round(mean(sigma_cpi,  na.rm = TRUE), 4),
    mean_sigma_bond = round(mean(sigma_bond, na.rm = TRUE), 4),
    min_sigma_gdp   = round(min(sigma_gdp,   na.rm = TRUE), 4),
    mean_r2_gdp     = round(mean(adj_r2_gdp, na.rm = TRUE), 3),
    mean_q1_gdp     = round(mean(fc_gdp_q1,  na.rm = TRUE), 3),
    mean_q1_cpi     = round(mean(fc_cpi_q1,  na.rm = TRUE), 3),
    mean_q1_bond    = round(mean(fc_bond_q1, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(n_endo) %>%
  print(width = Inf)

# Top 10 overall
cat("\n── Top 10 models overall (lowest GDP σ) ─────────────────────────────\n\n")
summary_tbl %>%
  filter(stable == TRUE, !is.na(serial_p), serial_p > 0.05) %>%
  arrange(sigma_gdp) %>%
  select(id, n_endo, oil_exog, sigma_gdp, sigma_cpi, sigma_bond,
         serial_p, max_root, fc_gdp_q1, fc_cpi_q1, fc_bond_q1) %>%
  head(10) %>%
  print(width = Inf)

# ── 6. Forecast plots for best model per set ──────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 5: FORECAST PLOTS — BEST MODEL PER SET\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Helper: one panel for one variable
plot_one_var <- function(res_obj, var_name, var_label, col,
                         nowcast_line = NULL, nowcast_label = NULL) {
  fc_mat <- res_obj$model$model  # not used directly
  mod    <- res_obj$model
  sname  <- res_obj$sample
  s      <- samples[[sname]]
  
  # Extract forecast for this variable from the stored predict() output
  fc_full <- tryCatch(
    predict(mod, n.ahead = N_AHEAD, ci = 0.95,
            dumvar = if (res_obj$set %in% names(ext_results) &&
                         any(sapply(var_sets, function(v) v$id == res_obj$set && v$oil)))
              exog_fc_with_oil else exog_fc_base),
    error = function(e) NULL
  )
  if (is.null(fc_full)) return(NULL)
  fc <- fc_full$fcst[[var_name]]
  
  hist_dat <- var_input %>%
    filter(date >= s$start_date) %>%
    tail(20) %>%
    select(date, value = !!sym(var_name))
  
  last_date <- max(hist_dat$date)
  fc_dates  <- seq(last_date + months(3), by = "quarter", length.out = nrow(fc))
  fc_df     <- tibble(date  = fc_dates,
                      fcst  = fc[, "fcst"],
                      lower = fc[, "lower"],
                      upper = fc[, "upper"])
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = COL_GREY, linewidth = 0.3) +
    geom_ribbon(data = fc_df, aes(x = date, ymin = lower, ymax = upper),
                fill = col, alpha = 0.15) +
    geom_line(data = hist_dat, aes(x = date, y = value),
              colour = col, linewidth = 0.8) +
    geom_point(data = hist_dat, aes(x = date, y = value),
               colour = col, size = 1.5) +
    geom_line(data = fc_df, aes(x = date, y = fcst),
              colour = col, linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = fc_df$date[1], y = fc_df$fcst[1] + 0.08,
             label = sprintf("%+.2f", fc_df$fcst[1]),
             size = 3.2, colour = col, fontface = "bold") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(title = var_label, x = NULL, y = NULL)
  
  # Add reference line if provided (e.g. Nowcasting Lab for GDP)
  if (!is.null(nowcast_line)) {
    p <- p +
      geom_hline(yintercept = nowcast_line, linetype = "dotted",
                 color = "gray50", linewidth = 0.4) +
      annotate("text", x = min(fc_df$date), y = nowcast_line + 0.08,
               label = nowcast_label, size = 2.5, color = "gray40", hjust = 0)
  }
  p
}

# Helper: 3-panel forecast plot for all three targets
plot_ext_forecast <- function(res_obj) {
  sname <- res_obj$sample
  s     <- samples[[sname]]
  
  p_gdp  <- plot_one_var(res_obj, "gdp_g",    "GDP growth (QoQ %)",
                         COL_GDP,  nowcast_line = 0.30,
                         nowcast_label = "Nowcasting Lab: +0.30%")
  p_cpi  <- plot_one_var(res_obj, "cpi_g",    "CPI inflation (QoQ %)",   COL_CPI)
  p_bond <- plot_one_var(res_obj, "bond_dif", "Bond yield change (QoQ pp)", COL_BOND)
  
  panels <- Filter(Negate(is.null), list(p_gdp, p_cpi, p_bond))
  if (length(panels) == 0) return(invisible(NULL))
  
  wrap_plots(panels, ncol = 3) +
    plot_annotation(
      title    = sprintf("VAR forecast — %s", res_obj$id),
      subtitle = sprintf("%s | p=%d | %d-var VARX%s | Q1 2026 GDP: %+.2f%%",
                         s$label, res_obj$p, res_obj$n_endo,
                         if (any(sapply(var_sets,
                                        function(v) v$id == res_obj$set && v$oil)))
                           " (Brent exog.)" else "",
                         res_obj$fc_gdp[1, "fcst"]),
      theme = theme(
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "gray40")
      )
    )
}

for (i in seq_len(nrow(best_per_set))) {
  row <- best_per_set[i, ]
  res <- ext_results[[row$set]][[row$sample]][[paste0("p", row$p)]]
  if (!is.null(res)) {
    fig <- tryCatch(plot_ext_forecast(res), error = function(e) NULL)
    if (!is.null(fig)) print(fig)
  }
}

# # Save outputs (uncomment)
# write.csv(summary_tbl,
#           file.path(here("output","tables"), "05_extended_summary.csv"),
#           row.names = FALSE)

cat("\nScript complete.\n")
cat(sprintf("Objects: ext_results (%d sets) | summary_tbl (%d rows)\n",
            length(ext_results), nrow(summary_tbl)))
cat("Next: 06_var_extended_evaluation.R\n")