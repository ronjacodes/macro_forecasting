# ============================================================================
# scripts/models/bvar/02_bvar_baseline_evaluation.R
#
# PURPOSE:
#   Pseudo out-of-sample evaluation of all BVAR baseline models from
#   01_bvar_baseline.R. Compares three priors against each other and
#   against AR(1) and frequentist VAR benchmarks.
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/var/01_var_baseline.R")
#   source("scripts/models/var/02_var_baseline_evaluation.R") — forecast_errors
#                                                               (AR1, VAR_post08_p3)
#   source("scripts/models/bvar/01_bvar_baseline.R")          — bvar_results,
#                                                               bvar_summary,
#                                                               best_per_prior,
#                                                               var_adj, PRIORS
#
# EVALUATION DESIGN:
#   Scheme    : Expanding window (same as all previous evaluations)
#   Start     : 2010 Q1
#   Horizons  : h = 1, 2, 4, 8
#   Targets   : gdp_g, cpi_g, bond_dif (all three)
#   COVID     : reported with and without 2020Q1–2021Q4 in error calc
#   Benchmarks: AR(1) per variable | frequentist VAR post08_p3
#
# MODELS EVALUATED:
#   Best model per prior (from best_per_prior table in script 01):
#   - bvar_minnesota_<sample>_p<lag>
#   - bvar_normalwishart_<sample>_p<lag>
#   - bvar_dummyobs_<sample>_p<lag>
#   + All models on the same sample × lag grid for full comparison
#
# NOTE ON COVID ADJUSTMENT IN EVALUATION:
#   The BVAR is estimated on COVID-adjusted series (var_adj). The forecast
#   errors are computed against the ORIGINAL (unadjusted) actuals in
#   var_input. This is the correct approach — we want to evaluate how well
#   the model forecasts the actual observed values, not the adjusted ones.
#
# STRUCTURE:
#   Section 0 — Parameters and benchmarks
#   Section 1 — BVAR expanding-window forecast loop
#   Section 2 — RMSE tables: all three variables, with/without COVID
#   Section 3 — Relative RMSE vs AR(1) and vs frequentist VAR
#   Section 4 — Diebold-Mariano tests: best BVAR per prior vs AR(1)
#   Section 5 — RMSE plots: all models + benchmarks by horizon
#   Section 6 — Best model per prior × target × horizon summary
#
# KEY FINDINGS (excluding COVID, expanding window 2010–2025):
#
#   ── GDP ──────────────────────────────────────────────────────────────────
#   BVAR baseline does NOT beat AR(1) at any horizon
#   h=1: all three BVARs ~258–261% of AR(1) RMSE — much worse
#   h=2: 118–121% | h=4: 113–115% | h=8: 101–102%
#   The h=1 failure is structural: BVAR on the full sample uses p=4 lags
#   but the COVID-adjusted training data still forces large initial forecasts
#   → Full sample p=4 BVAR is not a good GDP nowcasting tool
#   → DM tests: AR(1) significantly better at h=2 (**) for all priors
#   → vs freq VAR: BVAR slightly better at h=8 (~91% of freq VAR)
#
#   ── CPI ──────────────────────────────────────────────────────────────────
#   BVAR baseline near AR(1) for CPI, slightly worse at h=1,2
#   h=1: 105–106% of AR(1) | h=2: 106–107% (worse)
#   h=4: 99.2–99.8% (marginally better — not significant)
#   h=8: 98.8–99.2% (marginally better — best result for BVAR baseline)
#   DM tests: AR(1) significantly better at h=2 (**) — BVAR struggles short-term
#   BVAR beats freq VAR at h=1,2 (89–90% of freq VAR — avoids h=2 spike)
#   → Key finding: freq VAR had a severe h=2 CPI spike (0.487 RMSE vs 0.320
#     for AR1); BVAR completely avoids this due to shrinkage ✓
#   → This is the clearest benefit of Bayesian shrinkage shown in this project
#
#   ── Bond ─────────────────────────────────────────────────────────────────
#   BVAR baseline worse than AR(1) at h=1,2,4; marginally better at h=8
#   h=1: 113–115% | h=2: 104–105% | h=4: 101–102% | h=8: 99.3–100%
#   BVAR beats freq VAR at h=1,2 (95–96% of freq VAR) — freq VAR very poor
#   → Bond BVAR does NOT achieve the large improvements seen in extended VAR
#   → Extended frequentist VAR (with EA ESI) still dominates for bond
#
#   ── Prior comparison ─────────────────────────────────────────────────────
#   Three priors give nearly identical RMSE — differences < 0.5%
#   Minnesota and dummy obs tied at 5 wins each across target×horizon cells
#   Normal-Wishart: 2 wins (CPI h=2 and GDP h=1 — by tiny margins)
#   → Prior choice does not matter for out-of-sample performance here
#   → The gain from dummy obs (soc/sur) is not visible in the 3-variable
#     baseline — may matter more in larger extended models
#
#   ── Bottom line ──────────────────────────────────────────────────────────
#   BVAR BASELINE does not beat AR(1) for any target at any horizon
#   Key advantage over freq VAR: avoids the h=2 CPI spike (shrinkage helps)
#   Main lesson: BVAR baseline needs extended variables to improve forecasts
#   → Extended BVAR (script 04) likely needed to show real gains
#   → Stationarity of variables is fine; COVID adjustment approach B works
#
# Outputs (commented — uncomment to save):
#   output/tables/bvar_02_rmse_excl_covid.csv
#   output/tables/bvar_02_dm_tests.csv
#   output/figures/bvar_02_rmse_<variable>.png
# ============================================================================

library(dplyr)
library(sandwich)

TARGET_VARS   <- c("gdp_g", "cpi_g", "bond_dif")
TARGET_LABELS <- c(gdp_g    = "GDP growth (QoQ %)",
                   cpi_g    = "CPI inflation (QoQ %)",
                   bond_dif = "Bond yield change (QoQ pp)")

# ── 0. Parameters and benchmarks ─────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 0: PARAMETERS AND BENCHMARKS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

if (!exists("EVAL_START"))    EVAL_START    <- as.Date("2010-01-01")
if (!exists("HORIZONS"))      HORIZONS      <- c(1, 2, 4, 8)
if (!exists("COVID_EVAL_EXCL"))
  COVID_EVAL_EXCL <- seq(as.Date("2020-01-01"),
                         as.Date("2021-10-01"), by = "quarter")

cat(sprintf("Evaluation: %s – %s | h=%s\n",
            format(EVAL_START), format(max(var_input$date)),
            paste(HORIZONS, collapse = ",")))
cat(sprintf("COVID quarters excluded: %d\n\n", length(COVID_EVAL_EXCL)))

# Verify benchmarks
ar1_ok   <- all(paste0("AR1_", TARGET_VARS) %in% names(forecast_errors))
var_ok   <- "VAR_post08_p3" %in% names(forecast_errors)
cat(sprintf("AR(1) benchmarks: %s\n", if (ar1_ok) "available ✓" else "MISSING"))
cat(sprintf("Freq VAR post08_p3: %s\n\n",
            if (var_ok) "available ✓" else "MISSING — run 02_var_baseline_evaluation.R"))

# ── 1. BVAR expanding-window forecast loop ────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 1: BVAR EXPANDING-WINDOW FORECASTS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Evaluating best model per prior (from best_per_prior table)\n\n")

# Helper: build COVID-adjusted training data for a given sample + window
make_bvar_train <- function(start_date, up_to_date) {
  var_adj %>%
    filter(date >= start_date, date < up_to_date) %>%
    select(gdp_g    = gdp_g_adj,
           cpi_g    = cpi_g_adj,
           bond_dif = bond_dif_adj) %>%
    as.matrix()
}

# BVAR expanding-window function
bvar_expanding <- function(prior_name, p_var, sample_name,
                           eval_start, horizons) {
  s    <- samples[[sample_name]]
  rows <- list()
  
  all_idx  <- which(var_input$date >= s$start_date)
  eval_idx <- which(var_input$date >= eval_start)
  eval_idx <- intersect(eval_idx, all_idx)
  
  exog_fc_h <- matrix(0, nrow = max(horizons), ncol = ncol(exog_base),
                      dimnames = list(NULL, exog_cols_base))
  
  for (i in eval_idx) {
    for (h in horizons) {
      target_i <- i + h - 1
      if (target_i > nrow(var_input)) next
      
      Y_train <- make_bvar_train(s$start_date, var_input$date[i])
      n_train <- nrow(Y_train)
      if (n_train < p_var * 3 + 20) next
      
      mod <- tryCatch(
        bvar(data    = Y_train, lags   = p_var,
             n_draw  = 2000, n_burn = 1000,  # fewer draws for speed
             n_thin  = 1,
             priors  = PRIORS[[prior_name]],
             mh      = bv_metropolis(scale_hess = 0.01,
                                     adjust_acc  = TRUE,
                                     acc_lower   = 0.25,
                                     acc_upper   = 0.45),
             verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(mod)) next
      
      fc <- tryCatch(
        predict(mod, horizon = max(horizons)),
        error = function(e) NULL
      )
      if (is.null(fc)) next
      
      # Store errors for all three targets
      for (vi in seq_along(TARGET_VARS)) {
        vname  <- TARGET_VARS[vi]
        actual <- var_input[[vname]][target_i]   # ACTUAL unadjusted value
        fcst   <- mean(fc$fcast[, h, vi])        # posterior mean at horizon h
        
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

# Run evaluation
bvar_fc_errors <- list()

n_eval <- nrow(best_per_prior)
cat(sprintf("Evaluating %d models (best per prior)...\n\n", n_eval))

for (i in seq_len(n_eval)) {
  row        <- best_per_prior[i, ]
  prior_name <- row$prior
  mid        <- sprintf("bvar_%s_%s_p%d", prior_name, row$sample, row$p)
  
  cat(sprintf("  [%d/%d] %s ...\n", i, n_eval, mid))
  
  fc_rows <- tryCatch(
    bvar_expanding(prior_name  = prior_name,
                   p_var       = row$p,
                   sample_name = row$sample,
                   eval_start  = EVAL_START,
                   horizons    = HORIZONS),
    error = function(e) { message("  ERROR: ", e$message); NULL }
  )
  
  if (!is.null(fc_rows) && nrow(fc_rows) > 0) {
    bvar_fc_errors[[mid]] <- fc_rows %>% mutate(model = mid, prior = prior_name)
    cat(sprintf("    → %d error rows\n", nrow(fc_rows)))
  } else {
    cat("    → SKIPPED\n")
  }
}

cat(sprintf("\nEvaluation complete: %d models\n\n",
            length(bvar_fc_errors)))

# ── 2. RMSE tables ────────────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 2: RMSE TABLES\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Combine with AR1 and frequentist VAR benchmarks
all_bvar_fc <- bind_rows(
  # AR1
  bind_rows(lapply(TARGET_VARS, function(v)
    forecast_errors[[paste0("AR1_", v)]] %>% mutate(model = "AR1", prior = "benchmark")
  )),
  # Frequentist VAR
  if (var_ok)
    forecast_errors[["VAR_post08_p3"]] %>%
    mutate(model = "VAR_post08_p3", prior = "benchmark"),
  # BVAR models
  bind_rows(bvar_fc_errors)
)

compute_rmse <- function(fc_df, excl_covid = FALSE) {
  if (excl_covid) fc_df <- fc_df %>% filter(!date %in% COVID_EVAL_EXCL)
  fc_df %>%
    group_by(model, prior, variable, h) %>%
    summarise(n    = n(),
              rmse = sqrt(mean(error^2, na.rm = TRUE)),
              mae  = mean(abs(error),   na.rm = TRUE),
              .groups = "drop")
}

rmse_incl <- compute_rmse(all_bvar_fc, excl_covid = FALSE)
rmse_excl <- compute_rmse(all_bvar_fc, excl_covid = TRUE)

# Print RMSE tables — excluding COVID (primary)
cat("── RMSE excluding COVID (primary) ───────────────────────────────────\n")
for (vname in TARGET_VARS) {
  cat(sprintf("\n  %s:\n", TARGET_LABELS[vname]))
  rmse_excl %>%
    filter(variable == vname) %>%
    select(model, h, rmse) %>%
    pivot_wider(names_from = h, values_from = rmse, names_prefix = "h=") %>%
    mutate(across(starts_with("h="), ~ round(.x, 4))) %>%
    arrange(model) %>%
    print(n = Inf, width = Inf)
}

# ── 3. Relative RMSE ──────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 3: RELATIVE RMSE\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Values < 100 = beats the benchmark | Primary: excluding COVID\n\n")

make_rel <- function(rmse_tbl, bm_model) {
  bind_rows(lapply(TARGET_VARS, function(vname) {
    bm <- rmse_tbl %>%
      filter(model == bm_model, variable == vname) %>%
      select(h, bm_rmse = rmse)
    rmse_tbl %>%
      filter(variable == vname, !model %in% c("AR1","VAR_post08_p3")) %>%
      left_join(bm, by = "h") %>%
      mutate(rel_rmse  = round(rmse / bm_rmse * 100, 1),
             variable  = vname,
             benchmark = bm_model) %>%
      select(variable, benchmark, model, prior, h, rel_rmse)
  }))
}

rel_ar1 <- make_rel(rmse_excl, "AR1")
cat("── vs AR(1) ─────────────────────────────────────────────────────────\n")
for (vname in TARGET_VARS) {
  cat(sprintf("\n  %s:\n", TARGET_LABELS[vname]))
  rel_ar1 %>%
    filter(variable == vname) %>%
    select(model, h, rel_rmse) %>%
    pivot_wider(names_from = h, values_from = rel_rmse, names_prefix = "h=") %>%
    arrange(`h=1`) %>%
    print(n = Inf, width = Inf)
}

if (var_ok) {
  rel_var <- make_rel(rmse_excl, "VAR_post08_p3")
  cat("\n── vs frequentist VAR (post08_p3) ───────────────────────────────────\n")
  for (vname in TARGET_VARS) {
    cat(sprintf("\n  %s:\n", TARGET_LABELS[vname]))
    rel_var %>%
      filter(variable == vname) %>%
      select(model, h, rel_rmse) %>%
      pivot_wider(names_from = h, values_from = rel_rmse, names_prefix = "h=") %>%
      arrange(`h=1`) %>%
      print(n = Inf, width = Inf)
  }
}

# ── 4. Diebold-Mariano tests ──────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 4: DIEBOLD-MARIANO TESTS vs AR(1)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

dm_bvar <- function(err_bvar, err_ar1, h_val, excl_covid = TRUE) {
  if (excl_covid) {
    err_bvar <- err_bvar %>% filter(!date %in% COVID_EVAL_EXCL)
    err_ar1  <- err_ar1  %>% filter(!date %in% COVID_EVAL_EXCL)
  }
  joined <- inner_join(
    err_bvar %>% filter(h == h_val) %>% select(date, e_bvar = error),
    err_ar1  %>% filter(h == h_val) %>% select(date, e_ar1  = error),
    by = "date"
  ) %>% filter(!is.na(e_bvar), !is.na(e_ar1))
  
  if (nrow(joined) < 10)
    return(tibble(dm_stat = NA, p_value = NA, n = nrow(joined)))
  
  d <- joined$e_bvar^2 - joined$e_ar1^2
  tryCatch({
    fit    <- lm(d ~ 1)
    dm_hac <- coeftest(fit, vcov = NeweyWest(fit, lag = max(1, h_val - 1)))
    tibble(dm_stat = round(dm_hac[1, "t value"],    3),
           p_value = round(dm_hac[1, "Pr(>|t|)"],   4),
           n       = nrow(joined))
  }, error = function(e) tibble(dm_stat = NA, p_value = NA, n = nrow(joined)))
}

dm_rows <- list()
for (mid in names(bvar_fc_errors)) {
  prior_name <- bvar_fc_errors[[mid]]$prior[1]
  for (vname in TARGET_VARS) {
    err_bvar <- bvar_fc_errors[[mid]] %>% filter(variable == vname)
    err_ar1  <- forecast_errors[[paste0("AR1_", vname)]]
    for (h_val in HORIZONS) {
      res <- dm_bvar(err_bvar, err_ar1, h_val)
      dm_rows[[length(dm_rows) + 1]] <- tibble(
        model     = mid,
        prior     = prior_name,
        variable  = vname,
        h         = h_val
      ) %>%
        bind_cols(res) %>%
        mutate(
          sig = case_when(
            is.na(p_value)  ~ "n/a",
            p_value <= 0.01 ~ "***",
            p_value <= 0.05 ~ "**",
            p_value <= 0.10 ~ "*",
            TRUE            ~ "—"
          ),
          direction = case_when(
            is.na(dm_stat) ~ "n/a",
            dm_stat < 0    ~ "BVAR better",
            dm_stat > 0    ~ "AR1 better",
            TRUE           ~ "—"
          )
        )
    }
  }
}

bvar_dm_tbl <- bind_rows(dm_rows)

for (vname in TARGET_VARS) {
  cat(sprintf("── %s ────────────────────────────────────────\n",
              TARGET_LABELS[vname]))
  bvar_dm_tbl %>%
    filter(variable == vname, !is.na(dm_stat)) %>%
    select(model, prior, h, dm_stat, p_value, sig, direction) %>%
    arrange(h, p_value) %>%
    print(n = Inf, width = Inf)
  cat("\n")
}

# ── 5. RMSE plots ─────────────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 5: RMSE PLOTS BY HORIZON\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Colour scheme: benchmarks bold, BVAR priors in distinct colours
model_colours <- c(
  "AR1"          = "black",
  "VAR_post08_p3"= "gray40"
)
prior_colours <- c(
  minnesota     = COL_GDP,
  normalwishart = COL_CPI,
  dummyobs      = COL_BOND
)
# Map model IDs to prior colours
for (mid in names(bvar_fc_errors)) {
  pn <- bvar_fc_errors[[mid]]$prior[1]
  model_colours[mid] <- prior_colours[pn]
}

for (vname in TARGET_VARS) {
  plot_df <- rmse_excl %>%
    filter(variable == vname, h %in% HORIZONS) %>%
    mutate(
      h     = factor(h, levels = HORIZONS),
      is_bm = model %in% c("AR1","VAR_post08_p3"),
      lwd   = ifelse(is_bm, 1.4, 0.9),
      alpha = ifelse(is_bm, 1.0, 0.85),
      ltype = ifelse(is_bm, "solid", "dashed")
    )
  
  p <- ggplot(plot_df,
              aes(x = h, y = rmse, colour = model, group = model,
                  linewidth = lwd, alpha = alpha, linetype = ltype)) +
    geom_line() +
    geom_point(size = 2.5) +
    scale_colour_manual(values = model_colours, na.value = "grey70",
                        name = NULL) +
    scale_linewidth_identity() +
    scale_alpha_identity() +
    scale_linetype_identity() +
    labs(
      title    = sprintf("%s — BVAR RMSE by horizon", TARGET_LABELS[vname]),
      subtitle = paste0("Excl. COVID 2020Q1–2021Q4 | Expanding window 2010–2025\n",
                        "Blue=Minnesota | Red=Normal-Wishart | Green=Dummy obs | ",
                        "Black=AR(1) | Grey=Freq VAR"),
      x = "Forecast horizon (quarters)", y = "RMSE"
    ) +
    theme(legend.position = "right",
          legend.text     = element_text(size = 7))
  
  print(p)
}

# ── 6. Summary: best BVAR per target × horizon ───────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 6: BEST BVAR PER TARGET × HORIZON\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

bvar_best_by_target <- rmse_excl %>%
  filter(!model %in% c("AR1","VAR_post08_p3")) %>%
  group_by(variable, h) %>%
  slice_min(rmse, n = 1) %>%
  ungroup() %>%
  left_join(rel_ar1 %>% select(variable, model, h,
                               rel_vs_ar1 = rel_rmse),
            by = c("variable","model","h")) %>%
  arrange(variable, h)

print(bvar_best_by_target %>%
        select(variable, h, model, prior, rmse, rel_vs_ar1),
      n = Inf, width = Inf)

cat("\n── Does any BVAR beat AR(1)? ─────────────────────────────────────────\n\n")
bvar_best_by_target %>%
  mutate(beats_ar1 = rel_vs_ar1 < 100) %>%
  group_by(variable) %>%
  summarise(
    n_horizons_beats = sum(beats_ar1, na.rm = TRUE),
    best_rel         = round(min(rel_vs_ar1, na.rm = TRUE), 1),
    at_horizon       = h[which.min(rel_vs_ar1)],
    best_model       = model[which.min(rel_vs_ar1)],
    best_prior       = prior[which.min(rel_vs_ar1)],
    .groups = "drop"
  ) %>%
  print(width = Inf)

cat("\n── Which prior wins most often? ─────────────────────────────────────\n\n")
bvar_best_by_target %>%
  count(prior, name = "n_wins") %>%
  arrange(desc(n_wins)) %>%
  print()

# # Save outputs (uncomment)
# write.csv(rmse_excl,
#   file.path(here("output","tables"), "bvar_02_rmse_excl_covid.csv"),
#   row.names = FALSE)
# write.csv(bvar_dm_tbl,
#   file.path(here("output","tables"), "bvar_02_dm_tests.csv"),
#   row.names = FALSE)

cat("\nEvaluation complete.\n")
cat("Objects: rmse_incl, rmse_excl, rel_ar1, bvar_dm_tbl, bvar_best_by_target\n")
cat("Next: bvar/03_bvar_baseline_diagnostics.R\n")