# ============================================================================
# 02_var_baseline.R
# Baseline VAR models — three target variables, varying lag order and sample
#
# VARIABLES (all QoQ, stationary):
#   gdp_g   : Real GDP growth QoQ %      (KOF sports-event corrected)
#   cpi_g   : CPI inflation QoQ %        (headline)
#   bond_dif: 10Y bond yield QoQ change  (pp, first difference of level)
#
# EXOGENOUS: 8 COVID dummies (2020 Q1 – 2021 Q4), passed via exogen=
#
# MODELS:
#   Dimension 1 — Lag order  : p = 1, 2, 3, 4, 5
#   Dimension 2 — Sample     :
#     "full"   Q2 2000 – Q4 2025 (103 obs) — full available history
#     "post08" Q1 2009 – Q4 2025  (68 obs) — post-GFC, more homogeneous regime
#     "post15" Q1 2015 – Q4 2025  (44 obs) — low/negative rate era onwards
#
#   Total: 5 × 3 = 15 models
#
# SAMPLE RATIONALE:
#   Structural breaks visible in the data:
#   - 2008 GFC: sharp regime shift in bond yields and GDP volatility
#   - 2015: SNB drops EUR floor + enters negative rate territory
#   - 2020: COVID — handled via dummies, not sample truncation
#   We estimate all three samples to assess stability of results.
#   Note: post15 has only ~44 obs — use with caution, few df with p≥3.
#
# FOR EACH MODEL WE REPORT:
#   - Lag selection criteria (AIC/HQ/SC/FPE) for reference
#   - Coefficient summary (GDP equation focus)
#   - Diagnostics: serial correlation, normality, stability (max root)
#   - 4-quarter ahead forecast with 95% CI
#   - GDP Q1 2026 point forecast vs Nowcasting Lab benchmark (+0.30%)
#
# LAG SELECTION RESULTS:
#   Full sample   (n=103): AIC=1 HQ=1 SC=1 FPE=1  → unanimous p=1
#   Post-GFC      (n=68) : AIC=3 HQ=1 SC=1 FPE=3  → split; p=1 or p=3
#   Post-2015     (n=44) : AIC=8 HQ=8 SC=8 FPE=6  → all criteria want many lags
#                           but sample too small — only p=1,2 feasible
#
# FINDINGS (all 9 estimated models pass serial correlation and stability):
#   All models stable     : max root range 0.44 – 0.74, all well below 1 ✓
#   Serial correlation    : all pass (p range 0.10 – 0.79) ✓
#   Normality             : full sample fails (COVID-driven fat tails) ✗
#                           post-GFC p=1 passes (p=0.097), post-2015 p=1,2 pass ✓
#   GDP Q1 2026 range     : +0.368% to +0.525% across all models
#   Best model (full)     : full_p1  — sigma=0.478, serial p=0.600, Q1=+0.391%
#   Best model (post-GFC) : post08_p3 — sigma=0.357, serial p=0.726, Q1=+0.388%
#   Best model (post-2015): post15_p1 — sigma=0.349, serial p=0.399, Q1=+0.435%
#   vs Nowcasting Lab (+0.30%): all models forecast +0.07pp to +0.23pp above
#   Shorter samples → tighter CIs: full=1.87pp, post08=1.40pp, post15=1.37pp
#   Interpretation: removing the volatile 2000s regime improves fit substantially
#   (adj-R² GDP: 0.80 full → 0.91 post08 → 0.94 post15; sigma drops from 0.48→0.35)
#   The post-2015 sample has very few df for p≥3 — treat with caution.
#
# Prereq: source("00_setup.R")
# Outputs (commented — uncomment to save):
#   output/models/var_baseline_<id>.rds
#   output/figures/02_var_<id>_forecast.png
#   output/tables/02_var_baseline_summary.csv
# ============================================================================

# ── 0. Prepare VAR input data ─────────────────────────────────────────────────
# Join the three QoQ target series on quarterly dates
var_input <- gdp$pct_3m %>% rename(gdp_g = value) %>%
  inner_join(cpi$q_pct_3m  %>% rename(cpi_g    = value),   by = "date") %>%
  inner_join(bond_q_dif     %>% rename(bond_dif = bond_dif), by = "date") %>%
  arrange(date)

cat("VAR input: ", nrow(var_input), "quarters |",
    format(min(var_input$date)), "–", format(max(var_input$date)), "\n")
cat("Columns:", paste(names(var_input), collapse = ", "), "\n\n")

# COVID dummy matrix (always full-sample rows — subset later per sample)
covid_cols  <- paste0("d_", c("2020q1","2020q2","2020q3","2020q4",
                              "2021q1","2021q2","2021q3","2021q4"))
exog_full   <- as.matrix(var_data[match(var_input$date, var_data$date), covid_cols])

# ── 1. Define sample periods ──────────────────────────────────────────────────
samples <- list(
  full   = list(label = "Full sample (Q2 2000 – Q4 2025)",
                start_date = as.Date("2000-04-01"),
                ts_start   = c(2000, 2)),
  post08 = list(label = "Post-GFC (Q1 2009 – Q4 2025)",
                start_date = as.Date("2009-01-01"),
                ts_start   = c(2009, 1)),
  post15 = list(label = "Post-2015 / negative rate era (Q1 2015 – Q4 2025)",
                start_date = as.Date("2015-01-01"),
                ts_start   = c(2015, 1))
)

# ── 2. Helper functions ───────────────────────────────────────────────────────

# Subset var_input and exog to a given start date
subset_sample <- function(start_date) {
  idx  <- var_input$date >= start_date
  dat  <- var_input[idx, ]
  exog <- exog_full[idx, , drop = FALSE]
  list(dat = dat, exog = exog)
}

# Build ts object from data frame
make_ts <- function(dat, ts_start) {
  ts(dat[, c("gdp_g","cpi_g","bond_dif")],
     start = ts_start, frequency = 4)
}

# Run diagnostics and return as named list
run_diagnostics <- function(mod) {
  list(
    serial_p  = tryCatch(
      serial.test(mod, lags.pt = 10, type = "PT.asymptotic")$serial$p.value,
      error = function(e) NA_real_),
    normal_p  = tryCatch(
      normality.test(mod)$jb.mul$JB$p.value,
      error = function(e) NA_real_),
    max_root  = max(roots(mod)),
    stable    = max(roots(mod)) < 1,
    adj_r2_gdp = summary(mod)$varresult$gdp_g$adj.r.squared,
    sigma_gdp  = summary(mod)$varresult$gdp_g$sigma
  )
}

# Print a compact model summary
print_model <- function(id, mod, diag, fc_gdp) {
  cat("\n", strrep("─", 70), "\n")
  cat("MODEL:", id, "\n")
  cat(strrep("─", 70), "\n")
  cat(sprintf("  GDP eq. adj-R²: %.3f  |  Residual SD: %.4f\n",
              diag$adj_r2_gdp, diag$sigma_gdp))
  cat(sprintf("  Serial corr p:  %.4f  %s\n", diag$serial_p,
              ifelse(diag$serial_p > 0.05, "✓", "✗ autocorrelation present")))
  cat(sprintf("  Normality p:    %.4f  %s\n", diag$normal_p,
              ifelse(diag$normal_p > 0.05, "✓", "✗ non-normal residuals")))
  cat(sprintf("  Max root:       %.4f  %s\n", diag$max_root,
              ifelse(diag$stable, "✓ stable", "✗ UNSTABLE")))
  cat(sprintf("  GDP Q1 2026:    %+.3f%%  [%+.3f%%, %+.3f%%]\n",
              fc_gdp[1,"fcst"], fc_gdp[1,"lower"], fc_gdp[1,"upper"]))
  cat(sprintf("  CI width Q1:    %.3f pp\n",
              fc_gdp[1,"upper"] - fc_gdp[1,"lower"]))
  cat(sprintf("  vs NowcastLab:  %+.3f pp difference\n",
              fc_gdp[1,"fcst"] - 0.30))
}

# ── 3. Lag selection (run once per sample, for reference) ─────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("LAG SELECTION (VARselect, max p=8)\n")
cat("══════════════════════════════════════════════════════════════════════\n")

for (sname in names(samples)) {
  s   <- samples[[sname]]
  s_data <- subset_sample(s$start_date)
  en  <- make_ts(s_data$dat, s$ts_start)
  sel <- VARselect(en, lag.max = 8, type = "const", exogen = s_data$exog)
  cat(sprintf("\n%s [n=%d]:\n", s$label, nrow(s_data$dat)))
  cat("  AIC:", sel$selection["AIC(n)"],
      " HQ:", sel$selection["HQ(n)"],
      " SC:", sel$selection["SC(n)"],
      " FPE:", sel$selection["FPE(n)"], "\n")
}

# ── 4. Estimate all models ────────────────────────────────────────────────────
# Storage: results[[sample_id]][[paste0("p", lag)]]
results <- list()

# Forecast horizon
n_ahead  <- 8   # 2 years ahead
exog_fc  <- matrix(0, nrow = n_ahead, ncol = 8,
                   dimnames = list(NULL, covid_cols))

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("BASELINE VAR ESTIMATES — ALL MODELS\n")
cat("══════════════════════════════════════════════════════════════════════\n")

for (sname in names(samples)) {
  s      <- samples[[sname]]
  s_data <- subset_sample(s$start_date)
  en     <- make_ts(s_data$dat, s$ts_start)
  n_obs  <- nrow(s_data$dat)
  results[[sname]] <- list()
  
  cat("\n\n╔══════════════════════════════════════════════════════════════════╗\n")
  cat(sprintf("║  SAMPLE: %-56s║\n", s$label))
  cat(sprintf("║  N = %-3d observations                                          ║\n", n_obs))
  cat("╚══════════════════════════════════════════════════════════════════╝\n")
  
  for (p in 1:5) {
    id <- sprintf("%s_p%d", sname, p)
    
    # Skip if too few df (rough rule: need at least 3*p^2 + 8 dummies < n_obs)
    min_obs_needed <- 3 * p * (3 * p + 1) / 2 + 8 + 10  # generous buffer
    if (n_obs < min_obs_needed) {
      cat(sprintf("\n  [p=%d] SKIPPED — insufficient df (n=%d, need ~%d)\n",
                  p, n_obs, min_obs_needed))
      results[[sname]][[paste0("p", p)]] <- NULL
      next
    }
    
    # ── Estimate ──────────────────────────────────────────────────────────────
    mod <- tryCatch(
      VAR(en, p = p, type = "const", exogen = s_data$exog),
      error = function(e) { message("ERROR fitting ", id, ": ", e$message); NULL }
    )
    if (is.null(mod)) { results[[sname]][[paste0("p", p)]] <- NULL; next }
    
    # ── Diagnostics ───────────────────────────────────────────────────────────
    diag <- run_diagnostics(mod)
    
    # ── Forecast ──────────────────────────────────────────────────────────────
    fc     <- predict(mod, n.ahead = n_ahead, ci = 0.95, dumvar = exog_fc)
    fc_gdp <- fc$fcst$gdp_g
    
    # ── Print ─────────────────────────────────────────────────────────────────
    print_model(id, mod, diag, fc_gdp)
    
    # ── Store ─────────────────────────────────────────────────────────────────
    results[[sname]][[paste0("p", p)]] <- list(
      id       = id,
      model    = mod,
      diag     = diag,
      forecast = fc,
      fc_gdp   = fc_gdp,
      n_obs    = n_obs,
      p        = p,
      sample   = sname
    )
    
    # ── Save model object (uncomment to persist) ───────────────────────────────
    # saveRDS(results[[sname]][[paste0("p", p)]],
    #         file.path(here("output","models"), paste0("var_baseline_", id, ".rds")))
  }
}

# ── 5. Summary comparison table ───────────────────────────────────────────────
cat("\n\n══════════════════════════════════════════════════════════════════════\n")
cat("SUMMARY COMPARISON TABLE — ALL MODELS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

summary_rows <- list()
for (sname in names(samples)) {
  for (p in 1:5) {
    res <- results[[sname]][[paste0("p", p)]]
    if (is.null(res)) next
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      model           = res$id,
      sample          = sname,
      p               = res$p,
      n               = res$n_obs,
      adj_r2_gdp      = round(res$diag$adj_r2_gdp, 3),
      sigma_gdp       = round(res$diag$sigma_gdp,  4),
      serial_p        = round(res$diag$serial_p,   4),
      normality_p     = round(res$diag$normal_p,   4),
      max_root        = round(res$diag$max_root,   4),
      stable          = res$diag$stable,
      gdp_q1_fcst     = round(res$fc_gdp[1,"fcst"],  3),
      gdp_q1_lower    = round(res$fc_gdp[1,"lower"], 3),
      gdp_q1_upper    = round(res$fc_gdp[1,"upper"], 3),
      gdp_q1_ci_width = round(res$fc_gdp[1,"upper"] - res$fc_gdp[1,"lower"], 3),
      vs_nowcast      = round(res$fc_gdp[1,"fcst"] - 0.30, 3)
    )
  }
}

summary_table <- bind_rows(summary_rows)
print(summary_table, n = Inf, width = Inf)

# Flag models passing all diagnostics
cat("\nModels passing all diagnostics (serial p>0.05, stable, max_root<0.9):\n")
good <- summary_table %>%
  filter(serial_p > 0.05, stable, max_root < 0.90)
if (nrow(good) > 0) {
  print(good %>% select(model, serial_p, max_root, gdp_q1_fcst))
} else {
  cat("  None pass all criteria.\n")
}

# # Save summary table (uncomment to save)
# write.csv(summary_table,
#           file.path(here("output","tables"), "02_var_baseline_summary.csv"),
#           row.names = FALSE)

# ── 6. Forecast plot — best model per sample ──────────────────────────────────
# Best = lowest sigma_gdp among models passing serial test
best_models <- summary_table %>%
  filter(serial_p > 0.05, stable) %>%
  group_by(sample) %>%
  slice_min(sigma_gdp, n = 1) %>%
  ungroup()

cat("\nBest model per sample (serial p>0.05, min sigma GDP):\n")
print(best_models %>% select(model, p, adj_r2_gdp, sigma_gdp,
                             serial_p, gdp_q1_fcst, gdp_q1_ci_width))

# Helper: build one forecast panel for any variable
make_fc_panel <- function(fc_mat, hist_vec, hist_dates, var_label,
                          col_hist, col_fc, y_lab = "QoQ %",
                          ref_line = NULL, ref_label = NULL) {
  hist_df <- tibble(date  = hist_dates, value = hist_vec) %>%
    filter(date >= as.Date("2022-01-01"))
  fore_df <- tibble(
    date  = seq(as.Date("2026-01-01"), by = "quarter", length.out = n_ahead),
    value = fc_mat[, "fcst"],
    lower = fc_mat[, "lower"],
    upper = fc_mat[, "upper"]
  )
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = COL_GREY, linewidth = 0.35)
  if (!is.null(ref_line)) {
    p <- p + geom_hline(yintercept = ref_line, linetype = "dotted",
                        color = COL_BOND, linewidth = 0.7)
  }
  p <- p +
    geom_line(data = hist_df, aes(x = date, y = value),
              colour = col_hist, linewidth = 0.8) +
    geom_point(data = hist_df, aes(x = date, y = value),
               colour = col_hist, size = 1.8) +
    geom_ribbon(data = fore_df, aes(x = date, ymin = lower, ymax = upper),
                fill = col_fc, alpha = 0.15) +
    geom_line(data = fore_df, aes(x = date, y = value),
              colour = col_fc, linewidth = 1) +
    geom_point(data = fore_df, aes(x = date, y = value),
               colour = col_fc, size = 2) +
    annotate("text", x = as.Date("2026-01-01"),
             y = fc_mat[1,"fcst"] + diff(range(hist_df$value)) * 0.08,
             label = sprintf("%+.2f%%", fc_mat[1,"fcst"]),
             colour = col_fc, size = 3.2) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(title = var_label, x = NULL, y = y_lab)
  if (!is.null(ref_label)) {
    p <- p + annotate("text", x = as.Date("2024-01-01"),
                      y = ref_line + diff(range(hist_df$value)) * 0.06,
                      label = ref_label, color = COL_BOND, size = 2.8, hjust = 0)
  }
  p
}

# Plot all 3 variables for each best model
for (i in seq_len(nrow(best_models))) {
  bm     <- best_models[i, ]
  res    <- results[[bm$sample]][[paste0("p", bm$p)]]
  fc     <- res$forecast
  s_data <- subset_sample(samples[[bm$sample]]$start_date)
  hdates <- s_data$dat$date
  
  # GDP panel
  p_gdp <- make_fc_panel(
    fc_mat    = fc$fcst$gdp_g,
    hist_vec  = s_data$dat$gdp_g,
    hist_dates = hdates,
    var_label = "GDP growth (QoQ %)",
    col_hist  = COL_GDP, col_fc = COL_GDP,
    ref_line  = 0.30, ref_label = "Nowcasting Lab: +0.30%"
  )
  
  # CPI panel
  p_cpi <- make_fc_panel(
    fc_mat    = fc$fcst$cpi_g,
    hist_vec  = s_data$dat$cpi_g,
    hist_dates = hdates,
    var_label = "CPI inflation (QoQ %)",
    col_hist  = COL_CPI, col_fc = COL_CPI
  )
  
  # Bond yield change panel
  p_bond <- make_fc_panel(
    fc_mat    = fc$fcst$bond_dif,
    hist_vec  = s_data$dat$bond_dif,
    hist_dates = hdates,
    var_label = "Bond yield change (QoQ pp)",
    col_hist  = COL_BOND, col_fc = COL_BOND,
    y_lab     = "pp"
  )
  
  # Combine 3 panels
  fig <- (p_gdp | p_cpi | p_bond) +
    plot_annotation(
      title    = sprintf("VAR forecast — %s", bm$model),
      subtitle = sprintf("%s | p=%d | GDP adj-R²=%.2f | Serial p=%.3f",
                         samples[[bm$sample]]$label, bm$p,
                         bm$adj_r2_gdp, bm$serial_p),
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "gray40")
      )
    )
  
  print(fig)
  
  # # Save (uncomment)
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("02_var_forecast_%s.png", bm$model)),
  #        fig, width = 14, height = 5, dpi = 150)
}