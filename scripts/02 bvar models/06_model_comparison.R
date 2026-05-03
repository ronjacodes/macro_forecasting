# ============================================================================
# scripts/models/bvar/06_model_comparison.R
#
# PURPOSE:
#   Final model comparison across all approaches. Puts the best model
#   from each family side by side in a single RMSE table and forecast plot.
#
# COMPARISON DESIGN — apples to apples:
#
#   Block A — 3-variable baseline (same variables, different estimation):
#     AR(1)          : univariate benchmark, no COVID dummies
#     VAR baseline   : best 3-var VAR = post08_p3 (excl. COVID from RMSE)
#     BVAR baseline  : best 3-var BVAR = Minnesota full_p4 (incl. COVID)
#
#   Block B — extended models (adding variables, different estimation):
#     VARX extended  : best per target x horizon (excl. COVID from RMSE)
#                      ⚠ small n warning: h=8 bond n=3, h=8 CPI n=6
#     BVAR extended  : best per target x horizon (incl. COVID)
#
#   NOTE ON COVID TREATMENT:
#     VAR evaluation (scripts 02/06) excludes COVID quarters from RMSE.
#     BVAR evaluation (scripts 02b/04/05) includes COVID quarters.
#     This makes direct RMSE comparison imperfect — we flag this clearly
#     and focus on rankings rather than absolute RMSE values.
#     The AR(1) RMSE is computed both ways for reference.
#
#   NOTE ON VARX SMALL SAMPLES:
#     best_by_target_horizon uses very short windows for some models
#     (n=3 for bond h=8, n=6 for CPI h=8). These results are not
#     reliable and are flagged with ⚠ in the comparison table.
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/var/02_var_baseline_evaluation.R")
#     — forecast_errors, EVAL_START, HORIZONS, COVID_EVAL_EXCL
#   source("scripts/models/var/06_var_extended_evaluation.R")
#     — ext_forecast_errors, best_by_target_horizon
#   source("scripts/models/bvar/02_bvar_baseline_evaluation.R")
#     — bvar_fc_errors
#   source("scripts/models/bvar/04_bvar_extended_estimation.R")
#     — bvar_ext_fc_errors, bvar_ext_best, bvar_ext_summary
#   source("scripts/models/bvar/05_bvar_extended_evaluation.R")
#     — ext_bvar_rmse (optional, recomputed here if not available)
#
# KEY FINDINGS (fill in after running):
#   3-var block:
#     GDP  h=1: AR(1) wins | BVAR baseline close | VAR baseline worst
#     CPI  h=4: VAR baseline beats AR(1) | BVAR baseline similar
#     Bond h=8: VAR baseline beats AR(1) | BVAR baseline similar
#   Extended block:
#     CPI/Bond: VARX extended clearly wins (but small n caveat)
#     GDP: no model reliably beats AR(1)
#   Q1 2026 consensus: GDP ~+0.3-0.4%, CPI ~0%, Bond ~-0.1pp
#
# STRUCTURE:
#   Section 1 — 3-variable baseline RMSE comparison
#   Section 2 — Extended model RMSE comparison
#   Section 3 — Summary heatmap: relative RMSE vs AR(1)
#   Section 4 — Q1 2026 forecast comparison: all families, all targets
#
# Outputs (commented — uncomment to save):
#   output/figures/06_rmse_heatmap.png
#   output/figures/06_forecast_comparison.png
#   output/tables/06_rmse_comparison.csv
# ============================================================================

library(dplyr)

HORIZONS <- c(1, 2, 4, 8)

# ── helpers ───────────────────────────────────────────────────────────────────
rmse_from_tbl <- function(tbl, model_id, excl_covid = FALSE) {
  if (excl_covid) tbl <- tbl %>% filter(!date %in% COVID_EVAL_EXCL)
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
                                  na.rm = TRUE)), 4),
      n_gdp     = sum(sub$variable == "gdp_g"),
      n_cpi     = sum(sub$variable == "cpi_g"),
      n_bond    = sum(sub$variable == "bond_dif")
    )
  }))
}

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
      rmse_bond = round(rmse["bond_dif"], 4),
      n_gdp     = sum(!is.na(e[, "gdp_g"])),
      n_cpi     = sum(!is.na(e[, "cpi_g"])),
      n_bond    = sum(!is.na(e[, "bond_dif"]))
    )
  }))
}

# =============================================================================
# SECTION 1: 3-VARIABLE BASELINE COMPARISON
# =============================================================================
cat("=== SECTION 1: 3-VARIABLE BASELINE COMPARISON ===\n\n")
cat("Same 3 variables (GDP, CPI, Bond) — different estimation methods\n\n")

# AR(1) — both with and without COVID for reference
ar1_tbl <- bind_rows(
  forecast_errors[["AR1_gdp_g"]],
  forecast_errors[["AR1_cpi_g"]],
  forecast_errors[["AR1_bond_dif"]]
)
ar1_incl <- rmse_from_tbl(ar1_tbl, "AR(1) [incl COVID]", excl_covid = FALSE)
ar1_excl <- rmse_from_tbl(ar1_tbl, "AR(1) [excl COVID]", excl_covid = TRUE)

# VAR baseline (post08_p3) — excl COVID (consistent with var/02)
var_base <- rmse_from_tbl(
  forecast_errors[["VAR_post08_p3"]],
  "VAR baseline post08_p3", excl_covid = TRUE
)

# BVAR baseline (Minnesota full_p4) — incl COVID (consistent with bvar/02)
bvar_base <- rmse_from_tbl(
  bvar_fc_errors[["bvar_minnesota_full_p4"]],
  "BVAR baseline Minnesota full_p4", excl_covid = FALSE
)

baseline_rmse <- bind_rows(ar1_excl, ar1_incl, var_base, bvar_base)

cat("── RMSE table ───────────────────────────────────────────────────────\n")
baseline_rmse %>%
  select(model, h, rmse_gdp, rmse_cpi, rmse_bond) %>%
  arrange(h, rmse_gdp) %>%
  print(n = Inf, width = Inf)

# Relative to AR1 excl COVID
ar1_ref_excl <- ar1_excl %>%
  select(h, ar1_gdp = rmse_gdp, ar1_cpi = rmse_cpi, ar1_bond = rmse_bond)

cat("\n── Relative RMSE vs AR(1) excl COVID ───────────────────────────────\n")
baseline_rmse %>%
  filter(!grepl("AR\\(1\\)", model)) %>%
  left_join(ar1_ref_excl, by = "h") %>%
  mutate(
    rel_gdp  = round(rmse_gdp  / ar1_gdp,  3),
    rel_cpi  = round(rmse_cpi  / ar1_cpi,  3),
    rel_bond = round(rmse_bond / ar1_bond,  3)
  ) %>%
  select(model, h, rel_gdp, rel_cpi, rel_bond) %>%
  arrange(h) %>%
  print(n = Inf, width = Inf)

# =============================================================================
# SECTION 2: EXTENDED MODEL COMPARISON
# =============================================================================
cat("\n=== SECTION 2: EXTENDED MODEL COMPARISON ===\n\n")
cat("Adding variables — VARX (frequentist) vs BVAR extended (Bayesian)\n")
cat("NOTE: VARX uses excl-COVID RMSE; BVAR uses incl-COVID RMSE\n")
cat("NOTE: VARX small-n models flagged (n < 20 unreliable)\n\n")

# VARX best per target x horizon — from best_by_target_horizon
# Pick one representative model: most frequent in best_by_target_horizon
best_varx_id <- best_by_target_horizon %>%
  count(model, sort = TRUE) %>%
  slice(1) %>%
  pull(model)
cat(sprintf("Representative VARX model: %s\n", best_varx_id))

varx_best <- rmse_from_tbl(
  ext_forecast_errors[[best_varx_id]],
  sprintf("VARX best (%s)", best_varx_id),
  excl_covid = TRUE
)

# Also show the per-target best from best_by_target_horizon directly
cat("\nVARX best per target x horizon (from var/06):\n")
best_by_target_horizon %>%
  select(variable, h, model, rmse, n, rel_rmse_vs_ar1) %>%
  mutate(n_flag = ifelse(n < 20, " ⚠", "")) %>%
  arrange(variable, h) %>%
  print(n = Inf, width = Inf)

# BVAR extended best — best per set overall
bvar_ext_gdp_h1 <- sapply(names(bvar_ext_fc_errors), function(mid) {
  e <- bvar_ext_fc_errors[[mid]][["h1"]]
  if (is.null(e)) return(NA_real_)
  sqrt(mean(e[, "gdp_g"]^2, na.rm = TRUE))
})
best_bvar_ext_id <- names(which.min(bvar_ext_gdp_h1))
cat(sprintf("\nRepresentative extended BVAR model: %s\n\n", best_bvar_ext_id))

bvar_ext_best_rmse <- rmse_from_mat(
  bvar_ext_fc_errors[[best_bvar_ext_id]],
  sprintf("BVAR extended (%s)", best_bvar_ext_id)
)

ext_rmse <- bind_rows(ar1_excl, varx_best, bvar_ext_best_rmse)

cat("── Extended RMSE table ──────────────────────────────────────────────\n")
ext_rmse %>%
  select(model, h, rmse_gdp, rmse_cpi, rmse_bond) %>%
  arrange(h, rmse_gdp) %>%
  print(n = Inf, width = Inf)

# =============================================================================
# SECTION 3: SUMMARY HEATMAP
# =============================================================================
cat("\n=== SECTION 3: SUMMARY HEATMAP ===\n\n")

# Build clean comparison: use AR1 excl as reference throughout
# Use excl-COVID RMSE for VAR family, incl for BVAR family
# Label clearly

model_labels <- c(
  "VAR baseline"         = "VAR baseline\n(post08_p3, excl COVID)",
  "BVAR baseline"        = "BVAR baseline\n(Minnesota full_p4)",
  "VARX extended"        = sprintf("VARX extended\n(%s)", best_varx_id),
  "BVAR extended"        = sprintf("BVAR extended\n(%s)", best_bvar_ext_id)
)

all_comp <- bind_rows(
  var_base    %>% mutate(family = "VAR baseline"),
  bvar_base   %>% mutate(family = "BVAR baseline"),
  varx_best   %>% mutate(family = "VARX extended"),
  bvar_ext_best_rmse %>% mutate(family = "BVAR extended")
) %>%
  left_join(ar1_ref_excl, by = "h") %>%
  mutate(
    rel_gdp  = round(rmse_gdp  / ar1_gdp,  3),
    rel_cpi  = round(rmse_cpi  / ar1_cpi,  3),
    rel_bond = round(rmse_bond / ar1_bond,  3),
    family   = factor(family,
                      levels = c("VAR baseline", "BVAR baseline",
                                 "VARX extended", "BVAR extended"))
  )

rel_long <- all_comp %>%
  select(family, h, rel_gdp, rel_cpi, rel_bond) %>%
  pivot_longer(c(rel_gdp, rel_cpi, rel_bond),
               names_to  = "target",
               values_to = "ratio") %>%
  mutate(
    target = recode(target,
                    rel_gdp  = "GDP",
                    rel_cpi  = "CPI",
                    rel_bond = "Bond"),
    panel  = sprintf("h=%d", h),
    label  = sprintf("%.2f", ratio)
  )

fig_heat <- ggplot(rel_long,
                   aes(x = target, y = family, fill = ratio)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label  = label,
                colour = ifelse(ratio < 1, "white", "black")),
            size = 3, fontface = "bold") +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_gradient2(
    low      = COL_BOND,
    mid      = "white",
    high     = COL_CPI,
    midpoint = 1,
    limits   = c(0.3, 1.5),
    oob      = scales::squish,
    name     = "RMSE /\nAR(1)"
  ) +
  scale_colour_identity() +
  scale_y_discrete(limits = rev(levels(rel_long$family))) +
  labs(
    title    = "Relative RMSE vs AR(1) — all model families",
    subtitle = paste("Blue < 1 = beats AR(1) | AR(1) uses excl-COVID RMSE as reference",
                     "| VAR/VARX: excl-COVID | BVAR: incl-COVID"),
    x = NULL, y = NULL
  ) +
  theme(axis.text.x     = element_text(face = "bold", size = 10),
        axis.text.y     = element_text(size = 8),
        strip.text      = element_text(face = "bold", size = 10),
        legend.position = "right",
        plot.title      = element_text(face = "bold", size = 12),
        plot.subtitle   = element_text(size = 7, color = "gray40"))
print(fig_heat)
# ggsave(file.path(here("output","figures"), "06_rmse_heatmap.png"),
#        fig_heat, width = 12, height = 4, dpi = 150)

# =============================================================================
# SECTION 4: Q1 2026 FORECAST COMPARISON
# =============================================================================
cat("\n=== SECTION 4: Q1 2026 FORECAST COMPARISON ===\n\n")

# ── AR(1): last h=1 forecast from the expanding window ────────────────────────
ar1_fc <- bind_rows(
  forecast_errors[["AR1_gdp_g"]],
  forecast_errors[["AR1_cpi_g"]],
  forecast_errors[["AR1_bond_dif"]]
) %>%
  filter(h == 1) %>%
  group_by(variable) %>%
  slice_max(origin, n = 1) %>%
  ungroup() %>%
  select(variable, fcst) %>%
  pivot_wider(names_from = variable, values_from = fcst) %>%
  mutate(family = "AR(1)",
         fc_gdp = round(gdp_g, 3),
         fc_cpi = round(cpi_g, 3),
         fc_bond = round(bond_dif, 3)) %>%
  select(family, fc_gdp, fc_cpi, fc_bond)

# ── VAR baseline (post08_p3) ──────────────────────────────────────────────────
var_fc <- tryCatch({
  res <- results[["post08"]][["p3"]]
  exog_fc <- matrix(0, nrow = 1, ncol = ncol(exog_base),
                    dimnames = list(NULL, colnames(exog_base)))
  fc  <- predict(res$model, n.ahead = 1, ci = 0.95, dumvar = exog_fc)
  tibble(family  = "VAR baseline",
         fc_gdp  = round(fc$fcst$gdp_g[1,    "fcst"], 3),
         fc_cpi  = round(fc$fcst$cpi_g[1,    "fcst"], 3),
         fc_bond = round(fc$fcst$bond_dif[1, "fcst"], 3))
}, error = function(e) {
  cat(sprintf("  [!] VAR forecast: %s\n", e$message))
  tibble(family = "VAR baseline",
         fc_gdp = NA_real_, fc_cpi = NA_real_, fc_bond = NA_real_)
})

# ── VARX best (representative model) ─────────────────────────────────────────
varx_fc <- tryCatch({
  tbl <- ext_forecast_errors[[best_varx_id]] %>%
    filter(h == 1) %>%
    group_by(variable) %>%
    slice_max(origin, n = 1) %>%
    ungroup() %>%
    select(variable, fcst) %>%
    pivot_wider(names_from = variable, values_from = fcst)
  tibble(family  = sprintf("VARX extended\n(%s)", best_varx_id),
         fc_gdp  = round(tbl$gdp_g,    3),
         fc_cpi  = round(tbl$cpi_g,    3),
         fc_bond = round(tbl$bond_dif, 3))
}, error = function(e) {
  cat(sprintf("  [!] VARX forecast: %s\n", e$message))
  tibble(family = "VARX extended",
         fc_gdp = NA_real_, fc_cpi = NA_real_, fc_bond = NA_real_)
})

# ── BVAR baseline (Minnesota full_p4) ────────────────────────────────────────
bvar_base_fc <- tryCatch({
  res <- bvar_results[["minnesota"]][["full"]][["p4"]]
  fc  <- predict(res$model, horizon = 1, conf_bands = c(0.16, 0.84))
  tibble(family  = "BVAR baseline",
         fc_gdp  = round(median(fc$fcast[, 1, 1]), 3),
         fc_cpi  = round(median(fc$fcast[, 1, 2]), 3),
         fc_bond = round(median(fc$fcast[, 1, 3]), 3))
}, error = function(e) {
  cat(sprintf("  [!] BVAR baseline forecast: %s\n", e$message))
  tibble(family = "BVAR baseline",
         fc_gdp = NA_real_, fc_cpi = NA_real_, fc_bond = NA_real_)
})

# ── BVAR extended best ────────────────────────────────────────────────────────
bvar_ext_fc <- tryCatch({
  row <- bvar_ext_summary %>%
    filter(id == best_bvar_ext_id) %>% slice(1)
  tibble(family  = sprintf("BVAR extended\n(%s)", best_bvar_ext_id),
         fc_gdp  = row$fc_gdp_q1,
         fc_cpi  = row$fc_cpi_q1,
         fc_bond = row$fc_bond_q1)
}, error = function(e) {
  cat(sprintf("  [!] BVAR extended forecast: %s\n", e$message))
  tibble(family = "BVAR extended",
         fc_gdp = NA_real_, fc_cpi = NA_real_, fc_bond = NA_real_)
})

# ── Combine and print ─────────────────────────────────────────────────────────
fc_comparison <- bind_rows(ar1_fc, var_fc, varx_fc, bvar_base_fc, bvar_ext_fc)

cat("Q1 2026 point forecasts:\n\n")
fc_comparison %>%
  mutate(across(c(fc_gdp, fc_cpi, fc_bond),
                ~ sprintf("%+.3f", .x))) %>%
  print(width = Inf)

# ── Forecast comparison plot ──────────────────────────────────────────────────
family_order <- c("AR(1)", "VAR baseline",
                  sprintf("VARX extended\n(%s)", best_varx_id),
                  "BVAR baseline",
                  sprintf("BVAR extended\n(%s)", best_bvar_ext_id))

family_cols <- c(
  "AR(1)"        = "black",
  "VAR baseline" = "#7b3294"
)
family_cols[sprintf("VARX extended\n(%s)", best_varx_id)] <- "#e66101"
family_cols["BVAR baseline"] <- COL_GDP
family_cols[sprintf("BVAR extended\n(%s)", best_bvar_ext_id)] <- COL_CPI

fc_long <- fc_comparison %>%
  mutate(family = factor(family, levels = family_order)) %>%
  pivot_longer(c(fc_gdp, fc_cpi, fc_bond),
               names_to  = "target",
               values_to = "fcst") %>%
  mutate(
    target = recode(target,
                    fc_gdp  = "GDP growth (QoQ %)",
                    fc_cpi  = "CPI inflation (QoQ %)",
                    fc_bond = "Bond yield \u0394 (pp)"),
    target = factor(target,
                    levels = c("GDP growth (QoQ %)",
                               "CPI inflation (QoQ %)",
                               "Bond yield \u0394 (pp)"))
  )

fig_fc <- ggplot(fc_long,
                 aes(x = fcst, y = family, colour = family)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = COL_GREY, linewidth = 0.4) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%+.2f", fcst)),
            hjust = -0.3, size = 3, fontface = "bold") +
  facet_wrap(~ target, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = family_cols, guide = "none") +
  scale_y_discrete(limits = rev(family_order)) +
  labs(
    title    = "Q1 2026 forecasts — all model families",
    subtitle = "GDP growth | CPI inflation | Bond yield change",
    x = NULL, y = NULL
  ) +
  theme(strip.text         = element_text(face = "bold", size = 10),
        axis.text.y        = element_text(size = 8),
        axis.text.x        = element_blank(),
        axis.ticks.x       = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title         = element_text(face = "bold", size = 12))
print(fig_fc)
# ggsave(file.path(here("output","figures"), "06_forecast_comparison.png"),
#        fig_fc, width = 12, height = 4, dpi = 150)

# =============================================================================
# SAVE
# =============================================================================
# write.csv(bind_rows(baseline_rmse, ext_rmse),
#           file.path(here("output","tables"), "06_rmse_comparison.csv"),
#           row.names = FALSE)

cat("\nModel comparison complete.\n")
cat("Key objects: all_comp | fc_comparison | fig_heat | fig_fc\n")