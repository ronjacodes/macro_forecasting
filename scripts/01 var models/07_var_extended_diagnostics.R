# ============================================================================
# scripts/models/07_var_extended_diagnostics.R
#
# PURPOSE:
#   Deep diagnostics for the best extended VARX models identified in
#   06_var_extended_evaluation.R. Focus is on understanding WHAT the added
#   variables contribute — i.e. the new transmission channels — rather than
#   simply replicating the baseline diagnostics.
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/01_var_baseline.R")
#   source("scripts/models/05_var_extended_estimation.R")  — ext_results,
#                                                             var_sets, ext_data,
#                                                             exog_base,
#                                                             exog_with_oil
#   source("scripts/models/06_var_extended_evaluation.R")  — best_by_target_horizon
#
# FOCUS MODELS (from 06 evaluation — best per target):
#   GDP  : set_unemp_post08_p2    (only extended model beating AR1 for GDP)
#   CPI  : set_ea_post15_p1       (best h=4, good overall)
#          set_price_post08_p1    (best h=8 CPI)
#          set_world_post15_p1    (best h=2 CPI)
#   Bond : set_medium_full_p2     (best h=1 bond)
#          set_ea_post15_p1       (best h=2,4 bond)
#          set_oil_fx_post15_p1   (best h=8 bond)
#
#   We run full diagnostics on 3 representative models, one per target:
#     set_unemp_post08_p2    — GDP representative
#     set_ea_post15_p1       — CPI + bond representative (best both)
#     set_medium_full_p2     — bond h=1 representative (full sample)
#
# KEY QUESTIONS this script answers:
#   1. Do the new IRF channels make economic sense?
#      (EA ESI → GDP? Import prices → CPI? Unemployment → GDP?)
#   2. How much of CPI/bond variance is explained by the new variables?
#      (FEVD: does EA ESI contribute more than bond did in baseline?)
#   3. Are the extended models still stable despite more parameters?
#   4. Do residuals improve vs baseline (normality, serial correlation)?
#
# STRUCTURE:
#   Section 1 — Select focus models from ext_results
#   Section 2 — Impulse Response Functions (new variables → all targets)
#   Section 3 — Forecast Error Variance Decomposition
#   Section 4 — Rolling stability
#   Section 5 — Residual diagnostics
#
# CHOLESKY ORDERING (extended):
#   For models including EA ESI / world demand:
#     External (ea_esi / msci_wld) → gdp_g → cpi_g → bond_dif [→ fx_eur → imp_prc → unemp]
#   Rationale: global/external shocks affect Switzerland contemporaneously;
#              GDP affects CPI and bond within-quarter;
#              financial/FX variables ordered last (most endogenous)
#
# KEY FINDINGS:
#
#   ── IRF — transmission channels ─────────────────────────────────────────
#   unemp → gdp : +0.049 at h=1, decays smoothly to ~0 by h=8
#                 Sign counterintuitive — rising unemployment associated with
#                 slightly higher GDP in post-GFC Swiss data (labour hoarding /
#                 lagged adjustment). Small magnitude, no CI available.
#   unemp → bond: +0.049 at h=1, then near zero — unemployment rise slightly
#                 raises bond yields (risk-off? or inflation expectations?)
#   ea_esi → gdp: +0.027 h=1, peaks +0.033 at h=2, decays by h=8
#                 Positive and persistent — EA sentiment leads Swiss GDP ✓
#   ea_esi → cpi: +0.010 h=1, rising to +0.024 at h=2, slow decay
#                 EA demand pulls Swiss inflation upward with a lag ✓
#   ea_esi → bond:+0.069 h=1 (set_ea), largest response — EA sentiment
#                 strongly drives Swiss bond yield changes ✓
#                 Confirms: SNB follows ECB / EA financial conditions
#   In set_medium: ea_esi → gdp peaks +0.112 at h=2 (stronger, full sample)
#                 fx_eur → cpi: +0.042 h=1 (CHF weakening → higher import prices ✓)
#                 kof_bar → gdp: +0.026 h=4 (domestic leading indicator ✓)
#
#   ── FEVD — variance contributions ───────────────────────────────────────
#   set_unemp: unemployment explains only 3.6–3.9% of GDP variance at h=4-8
#              — small but consistent contribution
#   set_ea   : ea_esi explains 1.7–2.1% of GDP, 1.2–2.5% of CPI,
#              2.3–3.2% of bond variance at h=4-8
#              — modest in FEVD but forecasting improvement is large
#              → FEVD understates forecasting value (common in VAR literature)
#   set_medium: ea_esi explains 10.6–11.1% of GDP variance at h=4-8!
#               (largest new-variable contribution observed)
#               fx_eur explains 1.9% of CPI variance (imported inflation ✓)
#               CHF/EUR and KOF barometer explain <1% of bond variance
#
#   ── Rolling stability ────────────────────────────────────────────────────
#   All three models stable throughout (max root < 1) ✓
#   set_ea_post15:   root [0.751, 0.796] — most stable, small p=1 model
#   set_medium_full: root [0.846, 0.863] — stable despite 6 variables
#   set_unemp_post08:root [0.802, 0.926] — approaches caution zone (0.9)
#                    at some post-2022 windows — worth monitoring
#
#   ── Residual diagnostics ─────────────────────────────────────────────────
#   set_unemp_post08_p2: ALL pass JB normality and Ljung-Box ✓✓✓
#   set_ea_post15_p1   : ALL pass JB normality and Ljung-Box ✓✓✓
#   set_medium_full_p2 : GDP non-normal (JB p=0.0001 ✗) — COVID fat tails
#                        from full sample; CPI and bond normal ✓
#   → Post-GFC and post-2015 models have cleaner residuals than full-sample
#   → set_unemp and set_ea both achieve normality across all three equations
#
# Outputs (commented — uncomment to save):
#   output/figures/07_irf_<model>.png
#   output/figures/07_fevd_<model>.png
#   output/figures/07_stability_extended.png
#   output/figures/07_residuals_<model>.png
# ============================================================================

library(dplyr)

IRF_HORIZON <- 12
VAR_LABELS_EXT <- c(
  gdp_g    = "GDP growth",
  cpi_g    = "CPI inflation",
  bond_dif = "Bond yield Δ",
  fx_eur   = "CHF/EUR",
  ea_esi   = "EA Sentiment",
  msci_wld = "MSCI World",
  imp_prc  = "Import prices",
  kof_bar  = "KOF barometer",
  unemp    = "Unemployment"
)

# ── 1. Select focus models ────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 1: FOCUS MODELS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Pull models from ext_results
focus_ext <- list(
  set_unemp_post08_p2  = ext_results[["set_unemp"]][["post08"]][["p2"]],
  set_ea_post15_p1     = ext_results[["set_ea"]][["post15"]][["p1"]],
  set_medium_full_p2   = ext_results[["set_medium"]][["full"]][["p2"]]
)

# Drop any that are NULL
focus_ext <- Filter(Negate(is.null), focus_ext)

cat("Focus models:\n")
for (nm in names(focus_ext)) {
  res <- focus_ext[[nm]]
  cat(sprintf("  %-25s  %d-var  sample=%-8s  p=%d\n",
              nm, res$n_endo, res$sample, res$p))
  cat(sprintf("    endo: %s\n", paste(res$endo, collapse = ", ")))
}

# ── Helper: tidy IRF extraction (same as baseline diagnostics) ────────────────
tidy_irf <- function(irf_obj, impulse, response) {
  est   <- irf_obj$irf[[impulse]]
  lower <- irf_obj$Lower[[impulse]]
  upper <- irf_obj$Upper[[impulse]]
  col   <- which(colnames(est) == response)
  df    <- tibble(
    h        = 0:(nrow(est) - 1),
    impulse  = impulse,
    response = response,
    estimate = est[, col]
  )
  if (!is.null(lower) && !is.null(upper)) {
    df$lower <- lower[, col]
    df$upper <- upper[, col]
  } else {
    df$lower <- NA_real_
    df$upper <- NA_real_
  }
  df
}

# ── Helper: plot IRF panel ────────────────────────────────────────────────────
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
      title = sprintf("%s → %s",
                      VAR_LABELS_EXT[df$impulse[1]],
                      VAR_LABELS_EXT[df$response[1]]),
      x = "Quarters", y = "Response"
    )
}

# Colour map for variables
col_map <- c(
  gdp_g    = COL_GDP,
  cpi_g    = COL_CPI,
  bond_dif = COL_BOND,
  fx_eur   = "#7b3294",
  ea_esi   = "#e66101",
  msci_wld = "#d01c8b",
  imp_prc  = "#8c510a",
  kof_bar  = "#4393c3",
  unemp    = "#4dac26"
)

# ── 2. Impulse Response Functions ─────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 2: IMPULSE RESPONSE FUNCTIONS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Key focus: responses OF new variables TO target vars, and vice versa\n")
cat("Analytical CI (boot=TRUE fails with exogen in vars package)\n\n")

for (mname in names(focus_ext)) {
  res <- focus_ext[[mname]]
  mod <- res$model
  var_names <- res$endo
  
  cat(sprintf("── IRF: %s ─────────────────────────────────────────\n", mname))
  
  irf_obj <- tryCatch(
    irf(mod, n.ahead = IRF_HORIZON, boot = FALSE, ci = 0.90),
    error = function(e) { message("IRF error: ", e$message); NULL }
  )
  if (is.null(irf_obj)) next
  
  # Build all panels (k × k grid)
  panels <- list()
  for (imp in var_names) {
    for (resp in var_names) {
      df <- tryCatch(tidy_irf(irf_obj, imp, resp), error = function(e) NULL)
      if (is.null(df)) next
      panels[[paste(imp, resp, sep = "_")]] <-
        plot_irf_panel(df, col_line = col_map[resp])
    }
  }
  
  k <- length(var_names)
  fig_irf <- wrap_plots(panels, nrow = k, ncol = k,
                        byrow = TRUE) +
    plot_annotation(
      title    = sprintf("IRF — %s", mname),
      subtitle = sprintf("%s | p=%d | %d-var | 90%% asymptotic CI",
                         samples[[res$sample]]$label, res$p, k),
      theme = theme(
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "gray40")
      )
    )
  print(fig_irf)
  
  # Print key NEW-variable IRF values
  new_vars <- var_names[!var_names %in% BASE_COLS]
  target_vars <- c("gdp_g", "cpi_g", "bond_dif")
  
  cat("  Key responses — new variables → targets:\n")
  for (imp in new_vars) {
    for (resp in target_vars) {
      df <- tryCatch(tidy_irf(irf_obj, imp, resp), error = function(e) NULL)
      if (is.null(df)) next
      for (hh in c(1, 2, 4, 8)) {
        row <- df[df$h == hh, ]
        ci_zero <- is.na(row$lower) || (row$lower <= 0 & row$upper >= 0)
        cat(sprintf("    %s → %s h=%2d: %+.4f [%s] %s\n",
                    imp, resp, hh, row$estimate,
                    if (is.na(row$lower)) "no CI"
                    else sprintf("%.3f, %.3f", row$lower, row$upper),
                    if (!is.na(row$lower) && !ci_zero) "*** CI excl 0" else ""))
      }
    }
  }
  cat("\n")
  
  # # Save
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("07_irf_%s.png", mname)),
  #        fig_irf, width = 3.5 * k, height = 3 * k, dpi = 150)
}

# ── 3. Forecast Error Variance Decomposition ──────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 3: FORECAST ERROR VARIANCE DECOMPOSITION (FEVD)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("KEY QUESTION: how much variance in GDP/CPI/bond is explained by\n")
cat("the new variables vs. the baseline variables?\n\n")

for (mname in names(focus_ext)) {
  res <- focus_ext[[mname]]
  mod <- res$model
  var_names <- res$endo
  
  cat(sprintf("── FEVD: %s ─────────────────────────────────────────\n", mname))
  
  fevd_obj <- tryCatch(
    fevd(mod, n.ahead = IRF_HORIZON),
    error = function(e) { message("FEVD error: ", e$message); NULL }
  )
  if (is.null(fevd_obj)) next
  
  # Tidy FEVD
  fevd_tbl <- bind_rows(lapply(names(fevd_obj), function(resp) {
    mat <- fevd_obj[[resp]]
    as_tibble(mat) %>%
      mutate(h = 1:nrow(mat), response = resp) %>%
      pivot_longer(-c(h, response),
                   names_to = "impulse", values_to = "share")
  }))
  
  # Print at key horizons — focus on the three target variables
  cat("  FEVD at h=1,4,8 for target variables (% variance):\n")
  fevd_tbl %>%
    filter(response %in% BASE_COLS, h %in% c(1, 4, 8)) %>%
    mutate(share = round(share * 100, 1)) %>%
    pivot_wider(names_from = impulse, values_from = share) %>%
    arrange(response, h) %>%
    print(n = Inf, width = Inf)
  
  # Highlight new variable contributions
  new_vars <- var_names[!var_names %in% BASE_COLS]
  if (length(new_vars) > 0) {
    cat("\n  Contribution of NEW variables to target variance (h=4, h=8):\n")
    fevd_tbl %>%
      filter(response %in% BASE_COLS,
             impulse %in% new_vars,
             h %in% c(4, 8)) %>%
      mutate(share = round(share * 100, 1)) %>%
      arrange(response, h, impulse) %>%
      print(n = Inf, width = Inf)
  }
  
  # Stacked area plot — one panel per target variable
  col_fill <- col_map[var_names]
  names(col_fill) <- VAR_LABELS_EXT[var_names]
  
  fevd_plot <- fevd_tbl %>%
    filter(response %in% BASE_COLS) %>%
    mutate(
      response = factor(response, levels = BASE_COLS,
                        labels = c("GDP growth","CPI inflation","Bond yield Δ")),
      impulse  = VAR_LABELS_EXT[impulse]
    ) %>%
    ggplot(aes(x = h, y = share * 100, fill = impulse)) +
    geom_area(alpha = 0.85, colour = "white", linewidth = 0.2) +
    facet_wrap(~ response, ncol = 1) +
    scale_fill_manual(values = col_fill, name = "Shock from:") +
    scale_x_continuous(breaks = seq(1, IRF_HORIZON, by = 2)) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title    = sprintf("FEVD — %s", mname),
      subtitle = sprintf("%s | Cholesky: external → gdp → cpi → bond → financial",
                         samples[[res$sample]]$label),
      x = "Forecast horizon (quarters)",
      y = "% of forecast error variance"
    ) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  
  print(fevd_plot)
  
  # # Save
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("07_fevd_%s.png", mname)),
  #        fevd_plot, width = 10, height = 10, dpi = 150)
  cat("\n")
}

# ── 4. Rolling stability ───────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 4: ROLLING STABILITY — MAX CHARACTERISTIC ROOT\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

ROLL_START <- as.Date("2010-01-01")
roll_ext   <- list()

for (mname in names(focus_ext)) {
  res      <- focus_ext[[mname]]
  s        <- samples[[res$sample]]
  p        <- res$p
  endo_cols <- res$endo
  k         <- length(endo_cols)
  
  # Select exog
  vs_def   <- Filter(function(v) v$id == res$set, var_sets)[[1]]
  exog_mat <- if (vs_def$oil) exog_with_oil else exog_base
  
  all_idx  <- which(var_input$date >= s$start_date)
  roll_idx <- which(var_input$date >= ROLL_START)
  roll_idx <- intersect(roll_idx, all_idx)
  
  rows <- list()
  for (i in roll_idx) {
    extra_cols <- endo_cols[!endo_cols %in% BASE_COLS]
    
    train_dat <- var_input %>%
      select(date, all_of(BASE_COLS)) %>%
      { if (length(extra_cols) > 0)
        left_join(., ext_data %>% select(date, all_of(extra_cols)), by = "date")
        else . } %>%
      filter(date >= s$start_date, date <= var_input$date[i]) %>%
      filter(complete.cases(.[, endo_cols]))
    
    n_train  <- nrow(train_dat)
    n_params <- k * (k * p + 1) + ncol(exog_mat)
    if (n_train < n_params + 5) next
    
    exog_train <- exog_mat[var_input$date %in% train_dat$date, , drop = FALSE]
    
    ts_y <- as.numeric(format(train_dat$date[1], "%Y"))
    ts_q <- as.numeric(format(train_dat$date[1], "%m")) %/% 3 + 1
    en   <- ts(train_dat[, endo_cols],
               start = c(ts_y, ts_q), frequency = 4)
    
    mod_roll <- tryCatch(
      VAR(en, p = p, type = "const", exogen = exog_train),
      error = function(e) NULL
    )
    if (is.null(mod_roll)) next
    
    rows[[length(rows) + 1]] <- tibble(
      date     = var_input$date[i],
      model    = mname,
      max_root = max(roots(mod_roll)),
      n_train  = n_train
    )
  }
  
  roll_ext[[mname]] <- bind_rows(rows)
  if (nrow(roll_ext[[mname]]) > 0) {
    cat(sprintf("  %s: %d windows | root [%.3f, %.3f]\n",
                mname,
                nrow(roll_ext[[mname]]),
                min(roll_ext[[mname]]$max_root),
                max(roll_ext[[mname]]$max_root)))
  }
}

roll_df <- bind_rows(roll_ext)

if (nrow(roll_df) > 0) {
  fig_roll <- ggplot(roll_df, aes(x = date, y = max_root, colour = model)) +
    geom_hline(yintercept = 1.0, linetype = "dashed",
               color = "red", linewidth = 0.5) +
    geom_hline(yintercept = 0.9, linetype = "dotted",
               color = COL_CPI, linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    annotate("text", x = min(roll_df$date), y = 1.02,
             label = "Instability threshold", color = "red",
             size = 3, hjust = 0) +
    annotate("text", x = min(roll_df$date), y = 0.92,
             label = "Caution zone (0.9)", color = COL_CPI,
             size = 3, hjust = 0) +
    scale_colour_manual(
      values = c(set_unemp_post08_p2 = COL_GDP,
                 set_ea_post15_p1    = COL_CPI,
                 set_medium_full_p2  = COL_BOND),
      name = NULL
    ) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(limits = c(0, 1.1)) +
    labs(
      title    = "Rolling stability — extended VARX models",
      subtitle = "Expanding window | Max characteristic root per quarter",
      x = NULL, y = "Max characteristic root"
    ) +
    theme(legend.position = "bottom")
  
  print(fig_roll)
  
  # # Save
  # ggsave(file.path(here("output","figures"), "07_stability_extended.png"),
  #        fig_roll, width = 10, height = 5, dpi = 150)
}

# ── 5. Residual diagnostics ────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 5: RESIDUAL DIAGNOSTICS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Focus: do extended models have cleaner residuals than baseline?\n")
cat("Key comparison: JB normality and Ljung-Box serial correlation\n\n")

for (mname in names(focus_ext)) {
  res <- focus_ext[[mname]]
  mod <- res$model
  var_names <- res$endo
  
  cat(sprintf("── Residuals: %s ─────────────────────────────────────\n", mname))
  
  resid_mat <- resid(mod)
  # Only diagnose the three target variables
  diag_vars <- intersect(BASE_COLS, colnames(resid_mat))
  
  panels_acf <- list()
  panels_qq  <- list()
  
  for (vname in diag_vars) {
    r <- resid_mat[, vname]
    
    # ACF
    acf_vals <- acf(r, lag.max = 12, plot = FALSE)
    acf_df   <- tibble(
      lag = acf_vals$lag[-1],
      acf = acf_vals$acf[-1],
      ci  = qnorm(0.975) / sqrt(length(r))
    )
    p_acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
      geom_hline(yintercept = 0, color = COL_GREY) +
      geom_hline(yintercept =  acf_df$ci[1], linetype = "dashed",
                 color = COL_CPI, linewidth = 0.4) +
      geom_hline(yintercept = -acf_df$ci[1], linetype = "dashed",
                 color = COL_CPI, linewidth = 0.4) +
      geom_segment(aes(xend = lag, yend = 0),
                   colour = col_map[vname], linewidth = 0.8) +
      geom_point(colour = col_map[vname], size = 2) +
      scale_x_continuous(breaks = seq(1, 12, by = 2)) +
      labs(title = sprintf("ACF: %s", VAR_LABELS_EXT[vname]),
           x = "Lag", y = "ACF") +
      ylim(-0.5, 0.5)
    panels_acf[[vname]] <- p_acf
    
    # QQ
    qq_df <- tibble(
      theoretical = qnorm(ppoints(length(r))),
      sample      = sort(r)
    )
    p_qq <- ggplot(qq_df, aes(x = theoretical, y = sample)) +
      geom_abline(slope = sd(r), intercept = mean(r),
                  color = COL_CPI, linewidth = 0.7) +
      geom_point(colour = col_map[vname], alpha = 0.6, size = 1.5) +
      labs(title = sprintf("QQ: %s", VAR_LABELS_EXT[vname]),
           x = "Theoretical", y = "Sample")
    panels_qq[[vname]] <- p_qq
    
    # Stats
    jb <- tryCatch(tseries::jarque.bera.test(r), error = function(e) NULL)
    lb <- Box.test(r, lag = 10, type = "Ljung-Box")
    cat(sprintf("  %s: JB p=%.4f %s | Ljung-Box(10) p=%.4f %s\n",
                vname,
                if (!is.null(jb)) jb$p.value else NA,
                if (!is.null(jb) && jb$p.value > 0.05) "✓" else "✗ non-normal",
                lb$p.value,
                if (lb$p.value > 0.05) "✓" else "✗ autocorrelation"))
  }
  
  fig_resid <- (panels_acf[[1]] | panels_acf[[2]] | panels_acf[[3]]) /
    (panels_qq[[1]]  | panels_qq[[2]]  | panels_qq[[3]]) +
    plot_annotation(
      title    = sprintf("Residual diagnostics — %s", mname),
      subtitle = sprintf("%s | Top: ACF residuals | Bottom: Normal QQ",
                         samples[[res$sample]]$label),
      theme = theme(
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "gray40")
      )
    )
  
  print(fig_resid)
  
  # # Save
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("07_residuals_%s.png", mname)),
  #        fig_resid, width = 14, height = 8, dpi = 150)
  cat("\n")
}

cat("Diagnostics complete.\n")
cat("Next: 08_bvar_baseline.R\n")