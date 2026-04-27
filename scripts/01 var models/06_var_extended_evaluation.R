# ============================================================================
# scripts/models/06_var_extended_evaluation.R
#
# PURPOSE:
#   Pseudo out-of-sample evaluation of all extended VARX models from
#   05_var_extended_estimation.R. Compares against:
#     (a) AR(1) univariate benchmark
#     (b) Best baseline VAR (post08_p3 — identified in 02_var_baseline_evaluation.R)
#   for all three target variables: GDP, CPI, bond yield.
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/01_var_baseline.R")
#   source("scripts/models/02_var_baseline_evaluation.R") — forecast_errors (AR1),
#                                                           EVAL_START, HORIZONS,
#                                                           COVID_EVAL_EXCL
#   source("scripts/models/05_var_extended_estimation.R") — var_sets, ext_data,
#                                                           exog_base, exog_with_oil,
#                                                           exog_cols_base,
#                                                           exog_cols_with_oil
#
# EVALUATION DESIGN:
#   Scheme    : Expanding window (same as baseline evaluation)
#   Start     : 2010 Q1
#   Horizons  : h = 1, 2, 4, 8
#   Targets   : gdp_g, cpi_g, bond_dif (all three — equally weighted)
#   COVID     : reported with and without 2020Q1–2021Q4 in error calc
#   Benchmarks: AR(1) per variable | best baseline VAR (post08_p3)
#
# STRATEGY:
#   We only evaluate the BEST model per variable set (identified in script 05)
#   to keep runtime manageable. This means we evaluate ~18 models rather than
#   all 72 estimated models. The best model per set is defined as:
#   lowest σ(GDP), stable, serial correlation OK — from best_per_set table.
#
# STRUCTURE:
#   Section 0 — Parameters and benchmark forecast errors (AR1 already available)
#   Section 1 — Extended VARX expanding-window forecast loop
#   Section 2 — RMSE tables: all three variables, with/without COVID
#   Section 3 — Relative RMSE vs AR(1) and vs baseline VAR
#   Section 4 — Diebold-Mariano tests: best extended vs AR(1)
#   Section 5 — RMSE-by-horizon plots for each target variable
#   Section 6 — Summary: best model per target variable per horizon
#
# KEY FINDINGS (excluding COVID):
#
#   ── GDP ──────────────────────────────────────────────────────────────────
#   Extended models do NOT beat AR(1) for GDP at any horizon.
#   Best extended h=1: set_unemp_post08_p2 at 98.5% of AR(1) — barely under
#   All other models: 103–133% of AR(1) at h=1
#   vs baseline VAR: set_unemp_post08_p2 beats baseline (95.7% at h=1)
#   → For GDP: AR(1) remains the toughest benchmark.
#     Unemployment adds marginal value; other variables hurt.
#
#   ── CPI ──────────────────────────────────────────────────────────────────
#   Extended models clearly beat AR(1) at ALL horizons h=2,4,8.
#   h=1: most models worse than AR(1) (77–121%), except set_oil_fx (77%)
#   h=2: set_world_post15_p1 at 49.5% of AR(1) — halves the error
#   h=4: set_ea_post15_p1 at 48.8% of AR(1)
#   h=8: set_price_post08_p1 at 32.2% of AR(1) — remarkable improvement
#   vs baseline VAR: almost all models beat baseline for CPI at h=1,2
#     (baseline had h=2 spike to 0.487, extended models avoid this)
#   → Best CPI models: set_ea_post15_p1, set_world_post15_p1,
#     set_price_post08_p1, set_oil_fx_post15_p1
#   → CPI forecastability improves dramatically with the right variables
#
#   ── Bond yield ───────────────────────────────────────────────────────────
#   Extended models clearly beat AR(1) at ALL horizons.
#   h=1: set_medium_full_p2 at 44.8% of AR(1) — more than halves error
#        set_oil_fx_post15_p1 at 47.3%
#   h=2: set_ea_post15_p1 at 47.1%
#   h=4: set_ea_post15_p1 at 39.9% — best result overall
#   h=8: set_oil_fx_post15_p1 at 25.1% — extraordinary improvement
#   vs baseline VAR: most models beat baseline for bond at h=1,2
#     (set_medium_full_p2: 37.4% of baseline at h=1)
#   → Bond yield highly predictable with extended variables — especially
#     EA ESI, CHF/EUR, and the post-2015 sample
#   → DM tests: set_medium_full_p2 significantly better than AR1 at h=2 (**)
#
#   ── Bottom line ──────────────────────────────────────────────────────────
#   GDP  : extended models do not beat AR(1) — AR(1) dominates
#   CPI  : extended models beat AR(1) at h=2,4,8 — large improvements
#   Bond : extended models beat AR(1) at ALL horizons — very large improvements
#   The post-2015 sample models (set_ea, set_world, set_oil_fx) perform best
#   for CPI and bond, suggesting monetary regime matters for these variables.
#   The key new finding vs baseline: adding EA ESI and global variables
#   dramatically improves CPI and bond forecasting beyond what the 3-variable
#   baseline VAR could achieve.
#
# Outputs (commented — uncomment to save):
#   output/tables/06_rmse_extended_incl_covid.csv
#   output/tables/06_rmse_extended_excl_covid.csv
#   output/tables/06_rmse_relative_extended.csv
#   output/tables/06_dm_extended.csv
#   output/figures/06_rmse_extended_<variable>.png
# ============================================================================

library(dplyr)
library(sandwich)

TARGET_VARS  <- c("gdp_g", "cpi_g", "bond_dif")
TARGET_LABELS <- c(gdp_g = "GDP growth (QoQ %)",
                   cpi_g = "CPI inflation (QoQ %)",
                   bond_dif = "Bond yield change (QoQ pp)")

# ── 0. Parameters ─────────────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 0: PARAMETERS AND BENCHMARKS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# These should already exist from 02_var_baseline_evaluation.R:
# EVAL_START, HORIZONS, COVID_EVAL_EXCL, forecast_errors (contains AR1_*)
# If not, redefine:
if (!exists("EVAL_START")) EVAL_START      <- as.Date("2010-01-01")
if (!exists("HORIZONS"))   HORIZONS        <- c(1, 2, 4, 8)
if (!exists("COVID_EVAL_EXCL"))
  COVID_EVAL_EXCL <- seq(as.Date("2020-01-01"),
                         as.Date("2021-10-01"), by = "quarter")

cat(sprintf("Evaluation: %s – %s | horizons h=%s\n",
            format(EVAL_START),
            format(max(var_input$date)),
            paste(HORIZONS, collapse = ",")))
cat(sprintf("COVID quarters excluded from error calc: %d\n\n",
            length(COVID_EVAL_EXCL)))

# Verify AR1 benchmarks available (from 02_var_baseline_evaluation.R)
ar1_available <- all(paste0("AR1_", TARGET_VARS) %in% names(forecast_errors))
if (!ar1_available) {
  stop("AR(1) benchmark errors not found. Run 02_var_baseline_evaluation.R first.")
}
cat("AR(1) benchmarks: available ✓\n")

# Also need baseline VAR errors — post08_p3 is our best baseline
baseline_id <- "VAR_post08_p3"
baseline_available <- baseline_id %in% names(forecast_errors)
cat(sprintf("Baseline VAR (%s): %s\n\n",
            baseline_id,
            if (baseline_available) "available ✓" else
              "NOT FOUND — will compute from scratch"))

# ── 1. Extended VARX expanding-window forecast loop ───────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 1: EXTENDED VARX EXPANDING-WINDOW FORECASTS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Evaluating best model per variable set (from best_per_set table)\n\n")

# Helper: build training data for a given set + sample + window
make_ext_train <- function(endo_cols, start_date, up_to_date) {
  extra_cols <- endo_cols[!endo_cols %in% BASE_COLS]
  
  dat <- var_input %>%
    select(date, all_of(BASE_COLS)) %>%
    { if (length(extra_cols) > 0)
      left_join(., ext_data %>% select(date, all_of(extra_cols)), by = "date")
      else . } %>%
    filter(date >= start_date, date < up_to_date) %>%
    arrange(date)
  
  # Drop rows with NAs in endogenous cols
  complete_rows <- complete.cases(dat[, endo_cols])
  dat[complete_rows, ]
}

# Extended VARX expanding-window forecast function
ext_var_expanding <- function(vs, p_var, eval_start, horizons) {
  s         <- samples[[vs$sample_used]]
  endo_cols <- c(BASE_COLS, vs$extra)
  k         <- length(endo_cols)
  rows      <- list()
  
  # Select exog matrix
  use_oil  <- vs$oil
  exog_mat <- if (use_oil) exog_with_oil else exog_base
  exog_cols_use <- if (use_oil) exog_cols_with_oil else exog_cols_base
  
  # Forecast exog (all zeros)
  exog_fc_h <- matrix(0, nrow = max(horizons), ncol = ncol(exog_mat),
                      dimnames = list(NULL, exog_cols_use))
  
  all_idx  <- which(var_input$date >= s$start_date)
  eval_idx <- which(var_input$date >= eval_start)
  eval_idx <- intersect(eval_idx, all_idx)
  
  for (i in eval_idx) {
    for (h in horizons) {
      target_i <- i + h - 1
      if (target_i > nrow(var_input)) next
      
      # Build training data
      train_dat <- make_ext_train(endo_cols, s$start_date,
                                  var_input$date[i])
      n_train <- nrow(train_dat)
      
      # Min df check
      n_params <- k * (k * p_var + 1) + ncol(exog_mat)
      if (n_train < n_params + 5) next
      
      # Exog for training
      exog_train <- exog_mat[var_input$date %in% train_dat$date, ,
                             drop = FALSE]
      
      # Build ts
      ts_start_y <- as.numeric(format(train_dat$date[1], "%Y"))
      ts_start_q <- as.numeric(format(train_dat$date[1], "%m")) %/% 3 + 1
      en <- ts(train_dat[, endo_cols],
               start = c(ts_start_y, ts_start_q), frequency = 4)
      
      # Fit
      mod <- tryCatch(
        VAR(en, p = p_var, type = "const", exogen = exog_train),
        error = function(e) NULL
      )
      if (is.null(mod)) next
      
      # Forecast
      fc <- tryCatch(
        predict(mod, n.ahead = max(horizons), ci = 0.95,
                dumvar = exog_fc_h),
        error = function(e) NULL
      )
      if (is.null(fc)) next
      
      # Store errors for all three targets
      for (vname in TARGET_VARS) {
        actual <- var_input[[vname]][target_i]
        fcst   <- fc$fcst[[vname]][h, "fcst"]
        rows[[length(rows) + 1]] <- tibble(
          date     = var_input$date[target_i],
          origin   = var_input$date[i],
          h        = h,
          variable = vname,
          fcst     = fcst,
          actual   = actual,
          error    = actual - fcst
        )
      }
    }
  }
  bind_rows(rows)
}

# Run evaluation for best model per set
ext_forecast_errors <- list()

# Also run baseline VAR (post08_p3) if not already available
if (!baseline_available) {
  cat("Computing baseline VAR (post08_p3) forecast errors...\n")
  # Reuse var_expanding from 02_var_baseline_evaluation.R
  fc_base <- var_expanding(sname = "post08", p_var = 3,
                           eval_start = EVAL_START, horizons = HORIZONS)
  if (nrow(fc_base) > 0) {
    forecast_errors[[baseline_id]] <- fc_base
    cat("  Baseline VAR computed ✓\n\n")
  }
}

n_eval <- nrow(best_per_set)
cat(sprintf("Evaluating %d models (best per set)...\n\n", n_eval))

for (i in seq_len(n_eval)) {
  row    <- best_per_set[i, ]
  vs_def <- Filter(function(v) v$id == row$set, var_sets)[[1]]
  
  # Attach sample info needed by the helper
  vs_def$sample_used <- row$sample
  
  mid <- sprintf("%s_%s_p%d", row$set, row$sample, row$p)
  cat(sprintf("  [%d/%d] %s (n_endo=%d, oil=%s)...\n",
              i, n_eval, mid, row$n_endo,
              ifelse(row$oil_exog, "TRUE", "FALSE")))
  
  fc_ext <- tryCatch(
    ext_var_expanding(vs_def, p_var = row$p,
                      eval_start = EVAL_START, horizons = HORIZONS),
    error = function(e) {
      message("  ERROR: ", e$message); NULL
    }
  )
  
  if (!is.null(fc_ext) && nrow(fc_ext) > 0) {
    ext_forecast_errors[[mid]] <- fc_ext %>% mutate(model = mid)
    cat(sprintf("    → %d forecast-error rows\n", nrow(fc_ext)))
  } else {
    cat("    → SKIPPED\n")
  }
}

cat(sprintf("\nEvaluation complete: %d models\n\n",
            length(ext_forecast_errors)))

# ── 2. RMSE tables ────────────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 2: RMSE TABLES — ALL THREE TARGETS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Combine all forecast errors: AR1 + baseline VAR + extended
all_fc <- bind_rows(
  # AR1 benchmarks
  bind_rows(lapply(TARGET_VARS, function(v) {
    forecast_errors[[paste0("AR1_", v)]] %>% mutate(model = "AR1")
  })),
  # Baseline VAR
  if (baseline_id %in% names(forecast_errors)) {
    forecast_errors[[baseline_id]] %>% mutate(model = "baseline_VAR")
  },
  # Extended models
  bind_rows(ext_forecast_errors)
)

compute_rmse_tbl <- function(fc_df, excl_covid = FALSE) {
  if (excl_covid) fc_df <- fc_df %>% filter(!date %in% COVID_EVAL_EXCL)
  fc_df %>%
    group_by(model, variable, h) %>%
    summarise(n    = n(),
              rmse = sqrt(mean(error^2, na.rm = TRUE)),
              mae  = mean(abs(error),   na.rm = TRUE),
              .groups = "drop")
}

rmse_incl <- compute_rmse_tbl(all_fc, excl_covid = FALSE)
rmse_excl <- compute_rmse_tbl(all_fc, excl_covid = TRUE)

for (covid_label in c("incl", "excl")) {
  tbl <- if (covid_label == "incl") rmse_incl else rmse_excl
  cat(sprintf("\n── RMSE (%s COVID) ──────────────────────────────────────────\n",
              ifelse(covid_label == "incl", "including", "excluding")))
  for (vname in TARGET_VARS) {
    cat(sprintf("\n  %s:\n", TARGET_LABELS[vname]))
    tbl %>%
      filter(variable == vname) %>%
      select(model, h, rmse) %>%
      pivot_wider(names_from = h, values_from = rmse,
                  names_prefix = "h=") %>%
      mutate(across(starts_with("h="), ~ round(.x, 4))) %>%
      arrange(model) %>%
      print(n = Inf, width = Inf)
  }
}

# ── 3. Relative RMSE vs AR(1) and vs baseline VAR ────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 3: RELATIVE RMSE (AR1 = 100 | baseline VAR = 100)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Values < 100 = beats the benchmark. Primary comparison: excluding COVID.\n\n")

make_rel_rmse <- function(rmse_tbl, benchmark_model) {
  bind_rows(lapply(TARGET_VARS, function(vname) {
    bm <- rmse_tbl %>%
      filter(model == benchmark_model, variable == vname) %>%
      select(h, bm_rmse = rmse)
    
    rmse_tbl %>%
      filter(variable == vname,
             !model %in% c("AR1", "baseline_VAR")) %>%
      left_join(bm, by = "h") %>%
      mutate(rel_rmse  = round(rmse / bm_rmse * 100, 1),
             variable  = vname,
             benchmark = benchmark_model) %>%
      select(variable, benchmark, model, h, rel_rmse)
  }))
}

# Relative to AR1
rel_ar1 <- make_rel_rmse(rmse_excl, "AR1")

cat("── Relative to AR(1) — excluding COVID ─────────────────────────────\n")
for (vname in TARGET_VARS) {
  cat(sprintf("\n  %s:\n", TARGET_LABELS[vname]))
  rel_ar1 %>%
    filter(variable == vname) %>%
    select(model, h, rel_rmse) %>%
    pivot_wider(names_from = h, values_from = rel_rmse,
                names_prefix = "h=") %>%
    arrange(`h=1`) %>%
    print(n = Inf, width = Inf)
}

# Relative to baseline VAR
if (baseline_id %in% names(forecast_errors)) {
  rel_base <- make_rel_rmse(rmse_excl, "baseline_VAR")
  cat("\n── Relative to baseline VAR (post08_p3) — excluding COVID ──────────\n")
  for (vname in TARGET_VARS) {
    cat(sprintf("\n  %s:\n", TARGET_LABELS[vname]))
    rel_base %>%
      filter(variable == vname) %>%
      select(model, h, rel_rmse) %>%
      pivot_wider(names_from = h, values_from = rel_rmse,
                  names_prefix = "h=") %>%
      arrange(`h=1`) %>%
      print(n = Inf, width = Inf)
  }
}

# ── 4. Diebold-Mariano tests ──────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 4: DIEBOLD-MARIANO TESTS vs AR(1)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("H0: equal predictive accuracy | negative DM stat = extended model better\n\n")

dm_test_pair <- function(err_ext, err_ar1, h_val, excl_covid = TRUE) {
  if (excl_covid) {
    err_ext <- err_ext %>% filter(!date %in% COVID_EVAL_EXCL)
    err_ar1 <- err_ar1 %>% filter(!date %in% COVID_EVAL_EXCL)
  }
  joined <- inner_join(
    err_ext %>% filter(h == h_val) %>% select(date, e_ext = error),
    err_ar1 %>% filter(h == h_val) %>% select(date, e_ar1 = error),
    by = "date"
  ) %>% filter(!is.na(e_ext), !is.na(e_ar1))
  
  if (nrow(joined) < 10) return(tibble(dm_stat=NA, p_value=NA, n=nrow(joined)))
  
  d <- joined$e_ext^2 - joined$e_ar1^2
  tryCatch({
    fit    <- lm(d ~ 1)
    dm_hac <- coeftest(fit, vcov = NeweyWest(fit, lag = max(1, h_val - 1)))
    tibble(dm_stat = round(dm_hac[1,"t value"], 3),
           p_value = round(dm_hac[1,"Pr(>|t|)"], 4),
           n       = nrow(joined))
  }, error = function(e) tibble(dm_stat=NA, p_value=NA, n=nrow(joined)))
}

dm_rows <- list()
for (mid in names(ext_forecast_errors)) {
  for (vname in TARGET_VARS) {
    err_ext <- ext_forecast_errors[[mid]] %>% filter(variable == vname)
    err_ar1 <- forecast_errors[[paste0("AR1_", vname)]]
    for (h_val in HORIZONS) {
      res <- dm_test_pair(err_ext, err_ar1, h_val)
      dm_rows[[length(dm_rows) + 1]] <- tibble(
        model    = mid,
        variable = vname,
        h        = h_val
      ) %>% bind_cols(res) %>%
        mutate(
          sig = case_when(
            is.na(p_value) ~ "n/a",
            p_value <= 0.01 ~ "***",
            p_value <= 0.05 ~ "**",
            p_value <= 0.10 ~ "*",
            TRUE ~ "—"
          ),
          direction = case_when(
            is.na(dm_stat) ~ "n/a",
            dm_stat < 0    ~ "ext better",
            dm_stat > 0    ~ "AR1 better",
            TRUE           ~ "—"
          )
        )
    }
  }
}

dm_tbl <- bind_rows(dm_rows)

for (vname in TARGET_VARS) {
  cat(sprintf("── %s ─────────────────────────────────────────\n",
              TARGET_LABELS[vname]))
  dm_tbl %>%
    filter(variable == vname, !is.na(dm_stat)) %>%
    select(model, h, dm_stat, p_value, sig, direction) %>%
    arrange(h, p_value) %>%
    print(n = Inf, width = Inf)
  cat("\n")
}

# ── 5. RMSE-by-horizon plots ──────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 5: RMSE PLOTS BY HORIZON\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Colour scheme: AR1 black, baseline VAR dark grey, extended by set
model_colours <- c(
  "AR1"          = "black",
  "baseline_VAR" = "gray40"
)
# Add extended model colours — cycle through a palette
ext_models  <- unique(rmse_excl$model[!rmse_excl$model %in% c("AR1","baseline_VAR")])
ext_palette <- scales::hue_pal()(length(ext_models))
names(ext_palette) <- ext_models
model_colours <- c(model_colours, ext_palette)

for (vname in TARGET_VARS) {
  plot_df <- rmse_excl %>%
    filter(variable == vname, h %in% HORIZONS) %>%
    mutate(
      h     = factor(h, levels = HORIZONS),
      is_bm = model %in% c("AR1", "baseline_VAR"),
      lwd   = ifelse(is_bm, 1.2, 0.7),
      alpha = ifelse(is_bm, 1.0, 0.75)
    )
  
  p <- ggplot(plot_df,
              aes(x = h, y = rmse, colour = model,
                  group = model, linewidth = lwd, alpha = alpha)) +
    geom_line() +
    geom_point(size = 2) +
    scale_colour_manual(values = model_colours, na.value = "grey70", name = NULL) +
    scale_linewidth_identity() +
    scale_alpha_identity() +
    labs(
      title    = sprintf("%s — RMSE by forecast horizon", TARGET_LABELS[vname]),
      subtitle = "Excluding COVID 2020Q1–2021Q4 | Expanding window 2010 Q1–2025",
      x = "Forecast horizon (quarters)", y = "RMSE"
    ) +
    theme(legend.position = "right",
          legend.text     = element_text(size = 7))
  
  print(p)
  
  # # Save (uncomment)
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("06_rmse_extended_%s.png", vname)),
  #        p, width = 12, height = 5, dpi = 150)
}

# ── 6. Best model per target × horizon ───────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 6: BEST EXTENDED MODEL PER TARGET VARIABLE × HORIZON\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("(Excluding COVID, among extended models only — not including AR1 or baseline)\n\n")

best_by_target_horizon <- rmse_excl %>%
  filter(!model %in% c("AR1", "baseline_VAR")) %>%
  group_by(variable, h) %>%
  slice_min(rmse, n = 1) %>%
  ungroup() %>%
  left_join(rel_ar1 %>% select(variable, model, h, rel_rmse_vs_ar1 = rel_rmse),
            by = c("variable", "model", "h")) %>%
  arrange(variable, h)

print(best_by_target_horizon %>%
        select(variable, h, model, rmse, rel_rmse_vs_ar1),
      n = Inf, width = Inf)

cat("\n── Summary: does any extended model beat AR(1)? ─────────────────────\n\n")
best_by_target_horizon %>%
  mutate(beats_ar1 = rel_rmse_vs_ar1 < 100) %>%
  group_by(variable) %>%
  summarise(
    n_horizons_beats_ar1 = sum(beats_ar1, na.rm = TRUE),
    best_rel_rmse        = round(min(rel_rmse_vs_ar1, na.rm = TRUE), 1),
    at_horizon           = h[which.min(rel_rmse_vs_ar1)],
    best_model           = model[which.min(rel_rmse_vs_ar1)],
    .groups = "drop"
  ) %>%
  print(width = Inf)

# # Save outputs (uncomment)
# write.csv(rmse_incl,
#           file.path(here("output","tables"), "06_rmse_extended_incl_covid.csv"),
#           row.names = FALSE)
# write.csv(rmse_excl,
#           file.path(here("output","tables"), "06_rmse_extended_excl_covid.csv"),
#           row.names = FALSE)
# write.csv(dm_tbl,
#           file.path(here("output","tables"), "06_dm_extended.csv"),
#           row.names = FALSE)

cat("\nEvaluation complete.\n")
cat(sprintf("Objects: rmse_incl, rmse_excl, rel_ar1, dm_tbl, best_by_target_horizon\n"))
cat("Next: 07_var_extended_diagnostics.R\n")