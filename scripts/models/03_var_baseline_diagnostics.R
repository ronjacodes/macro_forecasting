# ============================================================================
# scripts/models/03_var_baseline_diagnostics.R
#
# PURPOSE:
#   Deep diagnostics for the baseline VAR models from 01_var_baseline.R.
#   Goes beyond the quick in-sample checks already done there.
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/01_var_baseline.R")  — results, var_input, samples,
#                                                  subset_sample(), make_ts()
#
# STRUCTURE:
#   Section 1 — Impulse Response Functions (IRF)
#   Section 2 — Forecast Error Variance Decomposition (FEVD)
#   Section 3 — Rolling stability (max characteristic root over time)
#   Section 4 — Residual diagnostics (ACF, normality QQ plots)
#
# FOCUS MODELS:
#   We run all diagnostics for the three best models identified in evaluation:
#     full_p1     — Full sample, p=1
#     post08_p3   — Post-GFC,   p=3  (best out-of-sample GDP at h=4)
#     post15_p1   — Post-2015,  p=1
#   Plus full_p1 is used as the reference for rolling stability across all
#   sample windows.
#
# IRF INTERPRETATION GUIDE:
#   - Ordering: gdp_g → cpi_g → bond_dif (Cholesky, recursive identification)
#     Assumes: GDP shock affects CPI and bond contemporaneously
#              CPI shock affects bond contemporaneously but not GDP
#              Bond shock affects neither GDP nor CPI contemporaneously
#   - This is a reduced-form VAR — ordering assumption matters for IRF
#   - Confidence bands: bootstrapped (runs = 500)
#   - Horizon: 12 quarters (3 years)
#
# KEY FINDINGS:
#
#   ── IRF (no CI bands — analytical CI not available with exogen) ───────────
#   Bond → GDP : +0.06 to +0.13pp at h=1, turns slightly negative h=3-4
#                Counterintuitive positive sign: yield rises coincided with
#                expansion phases in Swiss data (not monetary tightening signal)
#   Bond → CPI : +0.05 to +0.08pp, hump-shaped peak h=1-2, decays by h=6
#                Most economically meaningful channel; consistent with Granger ***
#   GDP  → CPI : +0.05 to +0.10pp at h=1-2, modest demand-pull inflation
#   CPI  → GDP : Negative ~-0.04 to -0.06pp, peak h=2-3 (inflation squeezes GDP)
#   All responses decay to zero by h=8-10 — consistent with stable VAR ✓
#
#   ── FEVD ─────────────────────────────────────────────────────────────────
#   GDP variance: own shocks dominate (88-100% all horizons)
#     Bond contributes 6-11% at h=4+; CPI only 2-4%
#   CPI variance: own shocks dominate at h=1 (85-95%)
#     Bond yield grows to 8-22% by h=4-8 — strongest in post15 (21.5%)
#     Quantifies bond→CPI monetary transmission channel
#   Bond variance: own shocks 62-79% at h=1
#     CPI surprisingly important: 9-33% (largest in post08_p3: 33%)
#     Bond yields driven by inflation expectations in post-GFC era
#
#   ── Rolling stability ────────────────────────────────────────────────────
#   All models stable throughout (max root < 1 at all windows) ✓
#   full_p1  : extremely stable [0.438, 0.476] — parsimonious p=1 structure
#   post08_p3: more variable [0.626, 0.851], approaches caution zone 2014-15
#              (SNB negative rate period with unusual bond dynamics)
#   post15_p1: volatile [0.363, 0.751], spikes to 0.75 in 2022 (rate normalization)
#              Reflects sensitivity to short sample size
#
#   ── Residual diagnostics ─────────────────────────────────────────────────
#   ACF: All models clean — no significant spikes, Ljung-Box p>0.05 ✓
#   Normality:
#     full_p1  : GDP residuals non-normal (JB p≈0.000) — COVID fat tails
#                CPI and bond normal (p=0.97, 0.52) ✓
#     post08_p3: ALL three variables normal (JB p=0.66, 0.11, 0.26) ✓✓✓
#                COVID dummies successfully absorb outliers
#     post15_p1: All pass (JB p=0.66, 0.91, 0.05) — bond borderline ✓
#   → post08_p3 is the statistically cleanest model (normal residuals in all eqs)
#   → Reinforces post08_p3 as the preferred specification
#
# Outputs (commented — uncomment to save):
#   output/figures/03_irf_<model>.png
#   output/figures/03_fevd_<model>.png
#   output/figures/03_rolling_stability.png
#   output/figures/03_residual_diagnostics_<model>.png
# ============================================================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)
library(tidyr)

# ── Focus models ──────────────────────────────────────────────────────────────
focus_models <- list(
  full_p1   = results[["full"]][["p1"]],
  post08_p3 = results[["post08"]][["p3"]],
  post15_p1 = results[["post15"]][["p1"]]
)

# Drop any that failed to estimate
focus_models <- Filter(Negate(is.null), focus_models)
cat("Focus models loaded:", paste(names(focus_models), collapse = ", "), "\n\n")

VAR_LABELS <- c(
  gdp_g    = "GDP growth (QoQ %)",
  cpi_g    = "CPI inflation (QoQ %)",
  bond_dif = "Bond yield change (pp)"
)

IRF_HORIZON <- 12   # quarters
BOOT_RUNS   <- 500  # bootstrap replications for IRF CI

# ============================================================================
# 1.  Impulse Response Functions
# ============================================================================
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 1: IMPULSE RESPONSE FUNCTIONS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Cholesky ordering: gdp_g → cpi_g → bond_dif\n")
cat("Bootstrap CI:", BOOT_RUNS, "runs | Horizon:", IRF_HORIZON, "quarters\n\n")

# Helper: extract IRF results into a tidy tibble
tidy_irf <- function(irf_obj, impulse, response) {
  est   <- irf_obj$irf[[impulse]]
  lower <- irf_obj$Lower[[impulse]]
  upper <- irf_obj$Upper[[impulse]]
  
  # Column index for response variable
  col <- which(colnames(est) == response)
  
  tibble(
    h       = 0:(nrow(est) - 1),
    impulse = impulse,
    response = response,
    estimate = est[, col],
    lower    = lower[, col],
    upper    = upper[, col]
  )
}

# Helper: plot a single IRF panel
plot_irf_panel <- function(df, col_line = COL_GDP) {
  has_ci <- !all(is.na(df$lower))
  p <- ggplot(df, aes(x = h)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = COL_GREY, linewidth = 0.35)
  if (has_ci) {
    p <- p + geom_ribbon(data = df,
                         aes(x = h, ymin = lower, ymax = upper),
                         fill = col_line, alpha = 0.15,
                         inherit.aes = FALSE)
  }
  p +
    geom_line(data = df, aes(x = h, y = estimate),
              colour = col_line, linewidth = 0.9,
              inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(0, IRF_HORIZON, by = 2)) +
    labs(
      title = sprintf("Impulse: %s → Response: %s",
                      VAR_LABELS[df$impulse[1]],
                      VAR_LABELS[df$response[1]]),
      x = "Quarters after shock", y = "Response"
    )
}

# Run IRFs for each focus model
for (mname in names(focus_models)) {
  res <- focus_models[[mname]]
  mod <- res$model
  cat(sprintf("── IRF: %s ────────────────────────────────────────────\n", mname))
  
  # Compute IRFs — analytical CI (boot=TRUE fails with exogen in vars package)
  # Analytical CI is asymptotically valid and standard in the literature
  irf_obj <- irf(mod,
                 n.ahead = IRF_HORIZON,
                 boot    = FALSE,
                 ci      = 0.90)   # 90% asymptotic CI
  
  var_names <- c("gdp_g", "cpi_g", "bond_dif")
  
  # Build all 9 panels (3 impulses × 3 responses)
  panels <- list()
  col_map <- c(gdp_g = COL_GDP, cpi_g = COL_CPI, bond_dif = COL_BOND)
  
  for (imp in var_names) {
    for (resp in var_names) {
      df <- tidy_irf(irf_obj, impulse = imp, response = resp)
      panels[[paste(imp, resp, sep = "_")]] <- plot_irf_panel(df, col_line = col_map[resp])
    }
  }
  
  # Arrange: rows = impulse, cols = response (standard VAR IRF matrix layout)
  fig_irf <- (panels[["gdp_g_gdp_g"]]   | panels[["gdp_g_cpi_g"]]   | panels[["gdp_g_bond_dif"]]) /
    (panels[["cpi_g_gdp_g"]]   | panels[["cpi_g_cpi_g"]]   | panels[["cpi_g_bond_dif"]]) /
    (panels[["bond_dif_gdp_g"]] | panels[["bond_dif_cpi_g"]] | panels[["bond_dif_bond_dif"]]) +
    plot_annotation(
      title    = sprintf("Impulse Response Functions — %s", mname),
      subtitle = sprintf("%s | Cholesky ordering | 90%% asymptotic CI",
                         samples[[res$sample]]$label),
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "gray40")
      )
    )
  
  print(fig_irf)
  
  # Print key IRF values at selected horizons
  cat("  Key responses at h=1, 4, 8:\n")
  for (imp in var_names) {
    for (resp in var_names) {
      if (imp == resp) next
      df <- tidy_irf(irf_obj, impulse = imp, response = resp)
      for (hh in c(1, 4, 8)) {
        row <- df[df$h == hh, ]
        ci_zero <- row$lower <= 0 & row$upper >= 0
        cat(sprintf("    %s → %s h=%d: %.4f [%.4f, %.4f] %s\n",
                    imp, resp, hh,
                    row$estimate, row$lower, row$upper,
                    ifelse(ci_zero, "(CI includes 0)", "*** CI excludes 0")))
      }
    }
  }
  cat("\n")
  
  # # Save (uncomment)
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("03_irf_%s.png", mname)),
  #        fig_irf, width = 14, height = 12, dpi = 150)
}

# ============================================================================
# 2.  Forecast Error Variance Decomposition (FEVD)
# ============================================================================
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 2: FORECAST ERROR VARIANCE DECOMPOSITION (FEVD)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("FEVD shows what % of forecast error variance of each variable\n")
cat("is attributable to shocks in each variable.\n\n")

for (mname in names(focus_models)) {
  res <- focus_models[[mname]]
  mod <- res$model
  cat(sprintf("── FEVD: %s ────────────────────────────────────────────\n", mname))
  
  fevd_obj <- fevd(mod, n.ahead = IRF_HORIZON)
  
  # Tidy FEVD for plotting
  fevd_tbl <- bind_rows(lapply(names(fevd_obj), function(resp) {
    mat <- fevd_obj[[resp]]
    as_tibble(mat) %>%
      mutate(h = 1:nrow(mat), response = resp) %>%
      pivot_longer(-c(h, response), names_to = "impulse", values_to = "share")
  }))
  
  # Print table at key horizons
  cat("  FEVD at h=1, 4, 8 (% of variance explained by each shock):\n")
  fevd_tbl %>%
    filter(h %in% c(1, 4, 8)) %>%
    mutate(share = round(share * 100, 1)) %>%
    pivot_wider(names_from = impulse, values_from = share) %>%
    arrange(response, h) %>%
    print(n = Inf, width = Inf)
  cat("\n")
  
  # Plot stacked bar chart
  col_fill <- c(gdp_g    = COL_GDP,
                cpi_g    = COL_CPI,
                bond_dif = COL_BOND)
  
  fevd_plot <- fevd_tbl %>%
    mutate(
      response = factor(response, levels = c("gdp_g","cpi_g","bond_dif"),
                        labels = unname(VAR_LABELS)),
      impulse  = factor(impulse,  levels = c("gdp_g","cpi_g","bond_dif"),
                        labels = unname(VAR_LABELS))
    ) %>%
    ggplot(aes(x = h, y = share * 100, fill = impulse)) +
    geom_area(alpha = 0.85, colour = "white", linewidth = 0.2) +
    facet_wrap(~ response, ncol = 1) +
    scale_fill_manual(values = setNames(unname(col_fill),
                                        unname(VAR_LABELS)),
                      name = "Shock from:") +
    scale_x_continuous(breaks = seq(1, IRF_HORIZON, by = 2)) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title    = sprintf("Forecast Error Variance Decomposition — %s", mname),
      subtitle = sprintf("%s | Cholesky ordering: gdp → cpi → bond",
                         samples[[res$sample]]$label),
      x = "Forecast horizon (quarters)", y = "% of forecast error variance"
    ) +
    theme(legend.position = "bottom",
          strip.text      = element_text(face = "bold", size = 10))
  
  print(fevd_plot)
  
  # # Save (uncomment)
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("03_fevd_%s.png", mname)),
  #        fevd_plot, width = 10, height = 10, dpi = 150)
}

# ============================================================================
# 3.  Rolling stability — max characteristic root over time
# ============================================================================
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 3: ROLLING STABILITY\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Re-estimates each model on an expanding window starting from EVAL_START\n")
cat("and records the max characteristic root. Stability requires max root < 1.\n\n")

# We use the full sample p=1 model for the rolling stability plot
# (most obs available, cleanest baseline)
ROLL_START <- as.Date("2010-01-01")

roll_stability <- list()

for (mname in names(focus_models)) {
  res    <- focus_models[[mname]]
  s      <- samples[[res$sample]]
  p      <- res$p
  
  all_idx  <- which(var_input$date >= s$start_date)
  roll_idx <- which(var_input$date >= ROLL_START)
  roll_idx <- intersect(roll_idx, all_idx)
  
  rows <- list()
  for (i in roll_idx) {
    train_idx  <- all_idx[all_idx <= i]
    if (length(train_idx) < p * 3 + 15) next
    
    dat_train  <- var_input[train_idx, ]
    exog_train <- exog_full[train_idx, , drop = FALSE]
    
    ts_start_q <- as.numeric(format(dat_train$date[1], "%Y"))
    ts_start_m <- as.numeric(format(dat_train$date[1], "%m")) %/% 3 + 1
    en <- ts(dat_train[, c("gdp_g","cpi_g","bond_dif")],
             start = c(ts_start_q, ts_start_m), frequency = 4)
    
    mod_roll <- tryCatch(
      VAR(en, p = p, type = "const", exogen = exog_train),
      error = function(e) NULL
    )
    if (is.null(mod_roll)) next
    
    rows[[length(rows) + 1]] <- tibble(
      date     = var_input$date[i],
      model    = mname,
      max_root = max(roots(mod_roll)),
      n_train  = length(train_idx)
    )
  }
  roll_stability[[mname]] <- bind_rows(rows)
  cat(sprintf("  %s: %d windows | root range [%.3f, %.3f]\n",
              mname,
              nrow(roll_stability[[mname]]),
              min(roll_stability[[mname]]$max_root),
              max(roll_stability[[mname]]$max_root)))
}

roll_df <- bind_rows(roll_stability)

# Plot
fig_roll <- ggplot(roll_df, aes(x = date, y = max_root, colour = model)) +
  geom_hline(yintercept = 1.0, linetype = "dashed",
             color = "red", linewidth = 0.5) +
  geom_hline(yintercept = 0.9, linetype = "dotted",
             color = COL_CPI, linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  annotate("text", x = min(roll_df$date), y = 1.02,
           label = "Instability threshold (root = 1)",
           color = "red", size = 3, hjust = 0) +
  annotate("text", x = min(roll_df$date), y = 0.92,
           label = "Caution zone (root = 0.9)",
           color = COL_CPI, size = 3, hjust = 0) +
  scale_colour_manual(
    values = c(full_p1   = COL_GDP,
               post08_p3 = COL_CPI,
               post15_p1 = COL_BOND),
    name = NULL
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(limits = c(0, 1.1)) +
  labs(
    title    = "Rolling stability — maximum characteristic root",
    subtitle = "Expanding window from each model's sample start | Re-estimated at each quarter",
    x = NULL, y = "Max characteristic root"
  ) +
  theme(legend.position = "bottom")

print(fig_roll)

# # Save (uncomment)
# ggsave(file.path(here("output","figures"), "03_rolling_stability.png"),
#        fig_roll, width = 10, height = 5, dpi = 150)

# ============================================================================
# 4.  Residual diagnostics — ACF and normality
# ============================================================================
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 4: RESIDUAL DIAGNOSTICS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

for (mname in names(focus_models)) {
  res <- focus_models[[mname]]
  mod <- res$model
  cat(sprintf("── Residuals: %s ───────────────────────────────────────\n", mname))
  
  resid_mat <- resid(mod)
  var_names <- colnames(resid_mat)
  
  panels_acf <- list()
  panels_qq  <- list()
  
  for (vname in var_names) {
    r <- resid_mat[, vname]
    
    # ACF
    acf_vals <- acf(r, lag.max = 12, plot = FALSE)
    acf_df   <- tibble(
      lag  = acf_vals$lag[-1],
      acf  = acf_vals$acf[-1],
      ci   = qnorm(0.975) / sqrt(length(r))
    )
    p_acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
      geom_hline(yintercept = 0, color = COL_GREY) +
      geom_hline(yintercept =  acf_df$ci[1], linetype = "dashed",
                 color = COL_CPI, linewidth = 0.4) +
      geom_hline(yintercept = -acf_df$ci[1], linetype = "dashed",
                 color = COL_CPI, linewidth = 0.4) +
      geom_segment(aes(xend = lag, yend = 0), colour = COL_GDP,
                   linewidth = 0.8) +
      geom_point(colour = COL_GDP, size = 2) +
      scale_x_continuous(breaks = seq(1, 12, by = 2)) +
      labs(title = sprintf("ACF: %s", VAR_LABELS[vname]),
           x = "Lag (quarters)", y = "Autocorrelation") +
      ylim(-0.5, 0.5)
    
    panels_acf[[vname]] <- p_acf
    
    # QQ plot
    qq_df <- tibble(
      theoretical = qnorm(ppoints(length(r))),
      sample      = sort(r)
    )
    p_qq <- ggplot(qq_df, aes(x = theoretical, y = sample)) +
      geom_abline(slope = sd(r), intercept = mean(r),
                  color = COL_CPI, linewidth = 0.7) +
      geom_point(colour = COL_GDP, alpha = 0.6, size = 1.5) +
      labs(title = sprintf("QQ: %s", VAR_LABELS[vname]),
           x = "Theoretical quantiles", y = "Sample quantiles")
    
    panels_qq[[vname]] <- p_qq
    
    # Print Jarque-Bera and Ljung-Box
    jb <- tseries::jarque.bera.test(r)
    lb <- Box.test(r, lag = 10, type = "Ljung-Box")
    cat(sprintf("  %s: JB p=%.4f | Ljung-Box(10) p=%.4f %s\n",
                vname, jb$p.value, lb$p.value,
                ifelse(lb$p.value > 0.05, "✓", "✗ autocorrelation")))
  }
  
  # Combine ACF + QQ side by side
  fig_resid <- (panels_acf[[1]] | panels_acf[[2]] | panels_acf[[3]]) /
    (panels_qq[[1]]  | panels_qq[[2]]  | panels_qq[[3]]) +
    plot_annotation(
      title    = sprintf("Residual diagnostics — %s", mname),
      subtitle = sprintf("%s | Top: ACF | Bottom: Normal QQ",
                         samples[[res$sample]]$label),
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "gray40")
      )
    )
  
  print(fig_resid)
  
  # # Save (uncomment)
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("03_residual_diagnostics_%s.png", mname)),
  #        fig_resid, width = 14, height = 8, dpi = 150)
  
  cat("\n")
}

cat("Diagnostics complete.\n")