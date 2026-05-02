# ============================================================================
# scripts/models/bvar/05_bvar_extended_evaluation.R
#
# PURPOSE:
#   Out-of-sample forecast evaluation for the extended BVAR models.
#   Compares extended BVAR against:
#     - AR(1) benchmark
#     - Baseline BVAR (Minnesota, best = full_p4)
#     - Best extended VARX (from var/06_var_extended_evaluation.R)
#   Mirrors the structure of var/06_var_extended_evaluation.R.
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/bvar/01_bvar_baseline.R")
#   source("scripts/models/bvar/02_bvar_baseline_evaluation.R")
#     — bvar_fc_errors: named list by model id, each a tibble with
#       cols: date, origin, h, variable, fcst, actual, error
#       e.g. bvar_fc_errors$bvar_minnesota_full_p4
#     — forecast_errors: named list by model id (AR1_gdp_g, AR1_cpi_g,
#       AR1_bond_dif, VAR_*), each a tibble with same columns
#   source("scripts/models/bvar/04_bvar_extended_estimation.R")
#     — bvar_ext_results, bvar_ext_summary, bvar_ext_best,
#       bvar_ext_fc_errors, bvar_sets, combined_data
#   source("scripts/models/var/06_var_extended_evaluation.R")  [optional]
#     — ext_forecast_errors: named list of tibbles (same structure as
#       forecast_errors above)
#
# NOTE ON DATA STRUCTURES (verified via str()):
#   forecast_errors$AR1_gdp_g: tibble [245 x 7]
#     cols: date, origin, h, fcst, actual, error, variable
#   bvar_fc_errors$bvar_minnesota_full_p4: tibble [735 x 9]
#     cols: date, origin, h, variable, fcst, actual, error, ...
#   bvar_ext_fc_errors$bset_ea_full_p4: named list $h1/$h2/$h4/$h8
#     each: matrix [windows, 3] of SQUARED errors (gdp_g, cpi_g, bond_dif)
#
# KEY FINDINGS:
#
#   Out-of-sample RMSE vs AR(1) — no extended BVAR beats AR(1):
#     GDP:  all extended models 1.22–1.57x AR(1) at all horizons
#     CPI:  all extended models 1.02–1.20x AR(1) at all horizons
#     Bond: extended models 0.94–1.24x AR(1) — one exception:
#           bset_oil_full beats AR(1) at h=4 (ratio=0.939)
#
#   vs baseline BVAR (Minnesota full_p4):
#     GDP:  extended BVAR worse at all horizons (1.35–1.51x baseline)
#     CPI:  extended BVAR worse at h=1,2; near-identical at h=4,8
#     Bond: extended BVAR marginally better at h=4 for some sets
#     → Adding variables hurts OOS despite better in-sample fit
#     → Bayesian shrinkage insufficient to offset parameter proliferation
#       with n=103 observations
#
#   vs best extended VARX (set_unemp_post08_p2):
#     VARX dominates at all horizons (ratio ~0.22–0.31 for GDP)
#     → Note: VARX reference model may benefit from short eval window
#       (post08_p2 starts later) — interpret with caution
#
#   Equal-weighted pool:
#     No improvement over individual models — pooling does not help
#     when all models share the same failure mode (OOS overfitting)
#
#   Q1 2026 forecasts (best per set, full_p4):
#     GDP:  +0.11% to +0.44% — consistent, above zero, near baseline
#     CPI:  near zero across all sets (-0.06% to +0.01%)
#           → consistent with continued Swiss disinflation
#     Bond: small negative (-0.09% to -0.19%) for most sets
#           → modest yield decline expected
#     All three targets consistent across variable sets —
#     extended BVAR forecasts are robust to variable choice
#
# STRUCTURE:
#   Section 1 — Build unified RMSE table
#   Section 2 — Relative RMSE (ratio to AR1)
#   Section 3 — Best extended BVAR per target x horizon + line plots
#   Section 4 — Equal-weighted forecast pool
#   Section 5 — Q1 2026 forecast summary and dot plot
#
# Outputs (commented — uncomment to save):
#   output/figures/bvar_05_rmse_heatmap.png
#   output/figures/bvar_05_rmse_lines.png
#   output/figures/bvar_05_fc_dot.png
#   output/tables/bvar_05_rmse_comparison.csv
# ============================================================================

library(dplyr)

HORIZONS <- c(1, 2, 4, 8)

# ── helpers ───────────────────────────────────────────────────────────────────

# RMSE from a tibble with cols: h, variable, error
rmse_from_tbl <- function(tbl, model_id) {
  bind_rows(lapply(HORIZONS, function(hh) {
    sub <- tbl %>% filter(h == hh, !is.na(error))
    if (nrow(sub) == 0) return(NULL)
    tibble(
      model     = model_id,
      h         = hh,
      rmse_gdp  = round(sqrt(mean(sub$error[sub$variable == "gdp_g"]^2,
                                  na.rm = TRUE)), 4),
      rmse_cpi  = round(sqrt(mean(sub$error[sub$variable == "cpi_g"]^2,
                                  na.rm = TRUE)), 4),
      rmse_bond = round(sqrt(mean(sub$error[sub$variable == "bond_dif"]^2,
                                  na.rm = TRUE)), 4)
    )
  }))
}

# RMSE from bvar_ext_fc_errors: list $h1/$h2/$h4/$h8 of squared-error matrices
rmse_from_mat <- function(errors, model_id) {
  bind_rows(lapply(HORIZONS, function(h) {
    e <- errors[[paste0("h", h)]]
    if (is.null(e)) return(NULL)
    rmse <- sqrt(colMeans(e, na.rm = TRUE))
    tibble(
      model     = model_id,
      h         = h,
      rmse_gdp  = round(rmse["gdp_g"],    4),
      rmse_cpi  = round(rmse["cpi_g"],    4),
      rmse_bond = round(rmse["bond_dif"], 4)
    )
  }))
}

# =============================================================================
# SECTION 1: BUILD UNIFIED RMSE TABLE
# =============================================================================
cat("=== SECTION 1: RMSE TABLE ===\n\n")

# ── AR(1): three separate tibbles, one per target ─────────────────────────────
ar1_tbl  <- bind_rows(
  forecast_errors[["AR1_gdp_g"]],
  forecast_errors[["AR1_cpi_g"]],
  forecast_errors[["AR1_bond_dif"]]
)
ar1_rmse <- rmse_from_tbl(ar1_tbl, "AR(1)")
cat("AR(1) RMSE:\n")
print(ar1_rmse, width = Inf)

# ── Baseline BVAR (Minnesota full_p4) ─────────────────────────────────────────
baseline_bvar_rmse <- rmse_from_tbl(
  bvar_fc_errors[["bvar_minnesota_full_p4"]],
  "BVAR baseline (Minnesota full_p4)"
)
cat("\nBaseline BVAR RMSE:\n")
print(baseline_bvar_rmse, width = Inf)

# ── Extended BVAR: matrix-list structure ──────────────────────────────────────
ext_bvar_rmse <- bind_rows(lapply(names(bvar_ext_fc_errors), function(mid) {
  rmse_from_mat(bvar_ext_fc_errors[[mid]], mid)
}))

# ── Best extended VARX from var/06 (optional) ─────────────────────────────────
varx_rmse <- tryCatch({
  varx_gdp_h1 <- sapply(names(ext_forecast_errors), function(mid) {
    sub <- ext_forecast_errors[[mid]] %>%
      filter(h == 1, variable == "gdp_g", !is.na(error))
    if (nrow(sub) == 0) return(NA_real_)
    sqrt(mean(sub$error^2, na.rm = TRUE))
  })
  best_varx_id <- names(which.min(varx_gdp_h1))
  cat(sprintf("\nBest VARX reference: %s\n", best_varx_id))
  rmse_from_tbl(ext_forecast_errors[[best_varx_id]],
                sprintf("VARX best (%s)", best_varx_id))
}, error = function(e) {
  cat(sprintf("\n[!] VARX unavailable: %s — skipping.\n", e$message))
  NULL
})

# ── Combine ───────────────────────────────────────────────────────────────────
all_rmse <- bind_rows(ar1_rmse, baseline_bvar_rmse, varx_rmse, ext_bvar_rmse)

cat("\n── Full RMSE table (sorted by h then GDP RMSE) ──────────────────────\n")
all_rmse %>%
  arrange(h, rmse_gdp) %>%
  print(n = Inf, width = Inf)

# =============================================================================
# SECTION 2: RELATIVE RMSE (ratio to AR1)
# =============================================================================
cat("\n=== SECTION 2: RELATIVE RMSE (ratio to AR(1)) ===\n\n")

ar1_ref  <- ar1_rmse %>%
  select(h, ar1_gdp = rmse_gdp, ar1_cpi = rmse_cpi, ar1_bond = rmse_bond)

rel_rmse <- all_rmse %>%
  filter(model != "AR(1)") %>%
  left_join(ar1_ref, by = "h") %>%
  mutate(
    rel_gdp  = round(rmse_gdp  / ar1_gdp,  3),
    rel_cpi  = round(rmse_cpi  / ar1_cpi,  3),
    rel_bond = round(rmse_bond / ar1_bond,  3)
  ) %>%
  select(model, h, rel_gdp, rel_cpi, rel_bond)

cat("Ratio to AR(1) — < 1 beats AR(1):\n\n")
rel_rmse %>%
  arrange(h, rel_gdp) %>%
  print(n = Inf, width = Inf)

# Heatmap
rel_long <- rel_rmse %>%
  pivot_longer(c(rel_gdp, rel_cpi, rel_bond),
               names_to  = "target",
               values_to = "ratio") %>%
  mutate(
    target = recode(target,
                    rel_gdp  = "GDP",
                    rel_cpi  = "CPI",
                    rel_bond = "Bond"),
    panel  = sprintf("h=%d", h)
  )

fig_heat <- ggplot(rel_long,
                   aes(x = target, y = model, fill = ratio)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label  = sprintf("%.2f", ratio),
                colour = ifelse(ratio < 1, "white", "black")),
            size = 2.6) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_gradient2(
    low      = COL_BOND,
    mid      = "white",
    high     = COL_CPI,
    midpoint = 1,
    limits   = c(0.5, 1.5),
    oob      = scales::squish,
    name     = "RMSE / AR(1)"
  ) +
  scale_colour_identity() +
  labs(
    title    = "Relative RMSE — extended BVAR vs AR(1)",
    subtitle = "Blue < 1 = beats AR(1) | Red > 1 = worse than AR(1)",
    x = NULL, y = NULL
  ) +
  theme(axis.text.x     = element_text(face = "bold"),
        axis.text.y     = element_text(size = 7),
        strip.text      = element_text(face = "bold"),
        legend.position = "right")
print(fig_heat)
# ggsave(file.path(here("output","figures"), "bvar_05_rmse_heatmap.png"),
#        fig_heat, width = 14, height = 6, dpi = 150)

# =============================================================================
# SECTION 3: BEST EXTENDED BVAR PER TARGET x HORIZON
# =============================================================================
cat("\n=== SECTION 3: BEST EXTENDED BVAR PER TARGET x HORIZON ===\n\n")

best_ext_by_th <- bind_rows(
  ext_bvar_rmse %>%
    group_by(h) %>% slice_min(rmse_gdp,  n = 1) %>% ungroup() %>%
    mutate(target = "gdp_g",    rmse = rmse_gdp),
  ext_bvar_rmse %>%
    group_by(h) %>% slice_min(rmse_cpi,  n = 1) %>% ungroup() %>%
    mutate(target = "cpi_g",    rmse = rmse_cpi),
  ext_bvar_rmse %>%
    group_by(h) %>% slice_min(rmse_bond, n = 1) %>% ungroup() %>%
    mutate(target = "bond_dif", rmse = rmse_bond)
) %>%
  left_join(ar1_ref, by = "h") %>%
  mutate(
    ar1_rmse  = case_when(
      target == "gdp_g"    ~ ar1_gdp,
      target == "cpi_g"    ~ ar1_cpi,
      target == "bond_dif" ~ ar1_bond
    ),
    ratio     = round(rmse / ar1_rmse, 3),
    beats_ar1 = ratio < 1
  ) %>%
  select(target, h, model, rmse, ar1_rmse, ratio, beats_ar1) %>%
  arrange(target, h)

print(best_ext_by_th, n = Inf, width = Inf)

# RMSE line plots
rmse_plot_df <- all_rmse %>%
  pivot_longer(c(rmse_gdp, rmse_cpi, rmse_bond),
               names_to  = "target",
               values_to = "rmse") %>%
  mutate(
    target     = recode(target,
                        rmse_gdp  = "GDP growth",
                        rmse_cpi  = "CPI inflation",
                        rmse_bond = "Bond yield \u0394"),
    model_type = case_when(
      model == "AR(1)"         ~ "ar1",
      grepl("baseline", model) ~ "baseline_bvar",
      grepl("VARX", model)     ~ "varx",
      TRUE                     ~ "ext_bvar"
    )
  )

fig_lines <- ggplot(rmse_plot_df, aes(x = h, y = rmse, group = model)) +
  geom_line(data = ~ filter(.x, model_type == "ext_bvar"),
            aes(colour = model), linewidth = 0.5, alpha = 0.5) +
  geom_line(data = ~ filter(.x, model_type == "ar1"),
            colour = "black", linewidth = 1.1, linetype = "dashed") +
  geom_line(data = ~ filter(.x, model_type == "baseline_bvar"),
            colour = COL_CPI, linewidth = 1.0, linetype = "dotted") +
  geom_line(data = ~ filter(.x, model_type == "varx"),
            colour = COL_BOND, linewidth = 1.0, linetype = "dotdash") +
  facet_wrap(~ target, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = HORIZONS) +
  labs(
    title    = "RMSE by horizon — extended BVAR vs benchmarks",
    subtitle = "Black dashed=AR(1) | Red dotted=baseline BVAR | Green=best VARX | Thin lines=ext BVAR",
    x        = "Forecast horizon (quarters)",
    y        = "RMSE",
    colour   = "Extended BVAR"
  ) +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 6),
        strip.text      = element_text(face = "bold"))
print(fig_lines)
# ggsave(file.path(here("output","figures"), "bvar_05_rmse_lines.png"),
#        fig_lines, width = 13, height = 5, dpi = 150)

# =============================================================================
# SECTION 4: EQUAL-WEIGHTED FORECAST POOL
# =============================================================================
cat("\n=== SECTION 4: EQUAL-WEIGHTED BVAR POOL ===\n\n")

pool_rmse <- ext_bvar_rmse %>%
  group_by(h) %>%
  summarise(
    model     = "BVAR pool (equal weights)",
    rmse_gdp  = round(mean(rmse_gdp,  na.rm = TRUE), 4),
    rmse_cpi  = round(mean(rmse_cpi,  na.rm = TRUE), 4),
    rmse_bond = round(mean(rmse_bond, na.rm = TRUE), 4),
    .groups   = "drop"
  ) %>%
  left_join(ar1_ref, by = "h") %>%
  mutate(
    rel_gdp  = round(rmse_gdp  / ar1_gdp,  3),
    rel_cpi  = round(rmse_cpi  / ar1_cpi,  3),
    rel_bond = round(rmse_bond / ar1_bond,  3)
  )

cat("Pool RMSE and ratio to AR(1):\n")
pool_rmse %>%
  select(h, rmse_gdp, rel_gdp, rmse_cpi, rel_cpi, rmse_bond, rel_bond) %>%
  print(width = Inf)

# =============================================================================
# SECTION 5: Q1 2026 FORECAST SUMMARY
# =============================================================================
cat("\n=== SECTION 5: Q1 2026 FORECAST SUMMARY ===\n\n")

cat("── All extended BVAR models ─────────────────────────────────────────\n")
bvar_ext_summary %>%
  select(id, n_endo, oil, sample, p,
         fc_gdp_q1, fc_gdp_q1_lo, fc_gdp_q1_hi,
         fc_cpi_q1, fc_bond_q1) %>%
  arrange(fc_gdp_q1) %>%
  print(n = Inf, width = Inf)

cat("\n── Range summary ────────────────────────────────────────────────────\n")
bvar_ext_summary %>%
  summarise(
    gdp_min   = round(min(fc_gdp_q1,  na.rm = TRUE), 3),
    gdp_mean  = round(mean(fc_gdp_q1, na.rm = TRUE), 3),
    gdp_max   = round(max(fc_gdp_q1,  na.rm = TRUE), 3),
    cpi_min   = round(min(fc_cpi_q1,  na.rm = TRUE), 3),
    cpi_mean  = round(mean(fc_cpi_q1, na.rm = TRUE), 3),
    cpi_max   = round(max(fc_cpi_q1,  na.rm = TRUE), 3),
    bond_min  = round(min(fc_bond_q1,  na.rm = TRUE), 3),
    bond_mean = round(mean(fc_bond_q1, na.rm = TRUE), 3),
    bond_max  = round(max(fc_bond_q1,  na.rm = TRUE), 3)
  ) %>%
  print(width = Inf)

# Q1 2026 forecast dot plot — best model per set, all three targets
fc_best_long <- bvar_ext_best %>%
  select(set, n_endo, oil, fc_gdp_q1, fc_gdp_q1_lo, fc_gdp_q1_hi,
         fc_cpi_q1, fc_bond_q1) %>%
  mutate(label = factor(gsub("bset_", "", set),
                        levels = gsub("bset_", "", set[order(fc_gdp_q1)])))

fig_dot <- ggplot(fc_best_long) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = COL_GREY, linewidth = 0.35) +
  # GDP with 68% interval
  geom_errorbarh(aes(y = label, xmin = fc_gdp_q1_lo, xmax = fc_gdp_q1_hi),
                 height = 0.2, alpha = 0.25, colour = COL_GDP) +
  geom_point(aes(x = fc_gdp_q1,  y = label, colour = "GDP growth"),
             size = 2.5) +
  # CPI
  geom_point(aes(x = fc_cpi_q1,  y = label, colour = "CPI inflation"),
             size = 2.5) +
  # Bond
  geom_point(aes(x = fc_bond_q1, y = label, colour = "Bond yield \u0394"),
             size = 2.5) +
  scale_y_discrete() +
  scale_colour_manual(
    values = c("GDP growth"        = COL_GDP,
               "CPI inflation"     = COL_CPI,
               "Bond yield \u0394" = COL_BOND),
    name = NULL
  ) +
  labs(
    title    = "Q1 2026 forecasts — best extended BVAR per variable set",
    subtitle = "Best = lowest in-sample \u03c3(GDP) | GDP band = 68% posterior interval",
    x        = "QoQ % (GDP/CPI) or pp change (Bond)",
    y        = NULL
  ) +
  theme(axis.text.y     = element_text(size = 8),
        legend.position = "bottom")
print(fig_dot)
# ggsave(file.path(here("output","figures"), "bvar_05_fc_dot.png"),
#        fig_dot, width = 10, height = 8, dpi = 150)

# =============================================================================
# SAVE
# =============================================================================
# write.csv(all_rmse,
#           file.path(here("output","tables"), "bvar_05_rmse_comparison.csv"),
#           row.names = FALSE)

cat("\nEvaluation complete.\n")
cat(sprintf("Objects: all_rmse (%d rows) | rel_rmse | best_ext_by_th\n",
            nrow(all_rmse)))
cat("Next: 06_model_comparison.R\n")