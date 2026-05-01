# ============================================================================
# scripts/models/bvar/01_bvar_baseline.R
#
# PURPOSE:
#   Estimate Bayesian VAR (BVAR) baseline models comparing THREE prior
#   specifications for the same 3-variable system: GDP growth, CPI inflation,
#   bond yield change. Priors are compared on the same model grid.
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/var/01_var_baseline.R")  — var_input, samples,
#                                                      exog_full, covid_cols
#   install.packages("BVAR")
#
# PACKAGE:
#   BVAR (Kuschnig & Vashold 2021) — Minnesota, Normal-Wishart, dummy obs
#   priors; Metropolis-Hastings sampler; predict/irf/fevd methods.
#   Citation: Kuschnig N, Vashold L (2021). Journal of Statistical Software.
#
# THREE PRIORS COMPARED:
#
#   1. MINNESOTA PRIOR (Litterman 1986) — baseline
#      Shrinks coefficients toward: own lag 1 = 1 (random walk), all else = 0
#      Lambda (tightness) optimised hierarchically via marginal likelihood
#      Most common prior in applied macro forecasting
#
#   2. NORMAL-WISHART CONJUGATE PRIOR
#      Fully conjugate extension of Minnesota: puts a joint prior on
#      coefficients AND residual covariance matrix Sigma jointly
#      (standard Minnesota conditions on Sigma — treats it as fixed)
#      Advantage: analytic posterior → faster sampling, proper uncertainty
#      on Sigma propagated to forecasts and IRFs
#      Reference: Kadiyala & Karlsson (1997), Karlsson (2013)
#
#   3. DUMMY OBSERVATIONS PRIOR (Banbura, Giannone & Reichlin 2010)
#      Implements Minnesota-type shrinkage by augmenting the data with
#      artificial "dummy" observations rather than specifying prior analytically
#      Two components:
#        - Sum-of-coefficients (SOC): shrinks toward variable-specific means
#          → captures persistence / near-unit-root behaviour
#        - Single-unit-root (SUR): shrinks toward common trend
#          → allows for cointegration-like behaviour without imposing it
#      Most common prior at central banks (ECB, SNB, Fed)
#      Particularly appropriate for Swiss macro given persistent dynamics
#
# COVID TREATMENT:
#   Approach B (same as before): OLS demean each series on COVID dummies
#   before BVAR estimation. Clean and comparable across all three priors.
#
# MODEL GRID:
#   Priors  : minnesota, normalwishart, dummyobs
#   Lag orders: p = 1, 2, 4
#   Samples : full, post08, post15
#   → 3 priors × 3 lags × 3 samples = up to 27 models
#     (many post15 skipped for df — especially p=4)
#
# MCMC SETTINGS:
#   n_draw = 10000, n_burn = 5000, n_thin = 1
#
# STRUCTURE:
#   Section 0 — Setup, COVID adjustment
#   Section 1 — Define prior configurations
#   Section 2 — Estimate all models
#   Section 3 — Summary table and prior comparison
#   Section 4 — Forecasts: fan charts for best model per prior
#   Section 5 — IRF with posterior bands: best model per prior
#
# KEY FINDINGS:
#
#   ── In-sample fit (sigma GDP) ────────────────────────────────────────────
#   All priors: best model is full sample p=4 (most lags, long history)
#   Dummy obs (full p4): σ=1.233  ← lowest across all 27 models
#   Minnesota  (full p4): σ=1.253
#   Norm-Wishart(full p4): σ=1.250
#   → Dummy obs prior gives best in-sample fit, consistent with stronger
#     shrinkage from extra soc/sur hyperparameters
#
#   ── Shrinkage hyperparameters ────────────────────────────────────────────
#   Lambda: 0.46–0.74 across all models (moderate shrinkage)
#   Dummy obs soc: ~1.4–1.9 (posterior > mode of 1 → moderate persistence)
#   Dummy obs sur: ~0.6–1.1 (posterior ≈ mode → mild common trend)
#   Minnesota and NW give near-identical lambda (prior choice barely matters
#   for tightness — the data dominates)
#
#   ── Forecast agreement (Q1 2026) ─────────────────────────────────────────
#   GDP:  +0.38% to +0.47% across best models — priors agree ✓
#   CPI:  +0.00% to +0.01% — all priors agree on near-zero CPI ✓
#   Bond: -0.06pp to -0.08pp — all priors agree on slight decline ✓
#   → Prior choice affects magnitude but NOT direction for any variable
#
#   ── IRF: all three priors give consistent transmission ───────────────────
#   GDP → CPI : positive, 68% band excludes 0 ✓ (all priors)
#   Bond → CPI: positive, hump-shaped (consistent across priors) ✓
#   GDP → Bond: persistent positive (bond rises after GDP shock) ✓
#   CPI → GDP : initially negative then positive — counterintuitive but
#               within posterior bands
#   → Prior robustness confirmed: IRF shapes are very similar across all
#     three priors, giving confidence in the economic interpretation
#
# Outputs (commented — uncomment to save):
#   output/figures/bvar_01_forecast_<id>.png
#   output/figures/bvar_01_irf_<prior>.png
#   output/tables/bvar_01_summary.csv
# ============================================================================

# ── 0. Setup and COVID adjustment ─────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 0: SETUP AND COVID ADJUSTMENT\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

if (!requireNamespace("BVAR", quietly = TRUE)) install.packages("BVAR")
library(BVAR)
library(dplyr)

cat("BVAR package version:", as.character(packageVersion("BVAR")), "\n\n")

TARGET_VARS <- c("gdp_g", "cpi_g", "bond_dif")

# COVID adjustment: OLS demean on COVID dummies, preserve mean level
covid_adjust <- function(y, exog_mat) {
  n   <- length(y)
  df  <- as.data.frame(exog_mat[seq_len(n), , drop = FALSE])
  df$y <- y
  fit  <- lm(y ~ . - y, data = df)
  fv   <- fitted(fit)
  y - fv + mean(fv)
}

var_adj <- var_input
for (vname in TARGET_VARS) {
  adj <- covid_adjust(var_input[[vname]], exog_full)
  var_adj[[paste0(vname, "_adj")]] <- adj
  cat(sprintf("  %s adjusted | max COVID effect removed: %.4f pp\n",
              vname, max(abs(var_input[[vname]] - adj))))
}
cat("\n")

N_DRAW <- 10000
N_BURN <- 5000
N_THIN <- 1
N_AHEAD <- 8
set.seed(2026)
cat(sprintf("MCMC: %d draws | %d burn-in | seed 2026\n\n", N_DRAW, N_BURN))

# ── 1. Prior configurations ───────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 1: PRIOR CONFIGURATIONS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# ── Prior 1: Minnesota (hierarchical lambda) ──────────────────────────────────
prior_minnesota <- bv_priors(
  hyper = "auto",
  mn    = bv_minnesota(
    lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 1e-4, max = 5),
    alpha  = bv_alpha(mode = 2),
    psi    = bv_psi()
  )
)
cat("Prior 1 — Minnesota (Litterman 1986):\n")
cat("  lambda ~ half-normal(mode=0.2) optimised via marginal likelihood\n")
cat("  alpha = 2 (quadratic lag decay), psi = OLS residual variance\n\n")

# ── Prior 2: Normal-Wishart conjugate ─────────────────────────────────────────
# Fully conjugate: joint prior on (B, Sigma)
# Uses same Minnesota-style coefficient shrinkage but also puts an
# inverse-Wishart prior on Sigma → proper uncertainty propagation
prior_normalwishart <- bv_priors(
  hyper = "auto",
  mn    = bv_minnesota(
    lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 1e-4, max = 5),
    alpha  = bv_alpha(mode = 2),
    psi    = bv_psi(
      scale  = 0.004,   # prior scale for Sigma — calibrated to QoQ data
      shape  = 0.004    # degrees of freedom for inv-Wishart
    )
  )
)
cat("Prior 2 — Normal-Wishart conjugate (Kadiyala & Karlsson 1997):\n")
cat("  Joint prior on (B, Sigma) via inverse-Wishart on Sigma\n")
cat("  psi: scale=0.004, shape=0.004 (diffuse but proper)\n")
cat("  Key difference vs Minnesota: Sigma uncertainty propagated to\n")
cat("  forecasts and IRFs → wider, more honest posterior bands\n\n")

# ── Prior 3: Dummy observations (Banbura, Giannone & Reichlin 2010) ───────────
# Implemented by augmenting the data matrix with artificial observations
# Two hyperparameters:
#   soc (sum-of-coefficients mu): shrinkage toward variable-specific RW
#   sur (single-unit-root delta): shrinkage toward common trend/level
prior_dummyobs <- bv_priors(
  hyper = "auto",
  mn    = bv_minnesota(
    lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 1e-4, max = 5),
    alpha  = bv_alpha(mode = 2),
    psi    = bv_psi()
  ),
  soc   = bv_soc(mode = 1, sd = 1, min = 1e-4, max = 50),
  sur   = bv_sur(mode = 1, sd = 1, min = 1e-4, max = 50)
)
cat("Prior 3 — Dummy observations (Banbura, Giannone & Reichlin 2010):\n")
cat("  Minnesota base + sum-of-coefficients (soc) + single-unit-root (sur)\n")
cat("  soc ~ half-normal(mode=1): shrinks toward variable-specific means\n")
cat("  sur ~ half-normal(mode=1): shrinks toward common stochastic trend\n")
cat("  All three hyperparameters (lambda, soc, sur) optimised jointly\n")
cat("  Standard prior at ECB, SNB and other central banks\n\n")

# Collect priors in named list for loop
PRIORS <- list(
  minnesota    = prior_minnesota,
  normalwishart = prior_normalwishart,
  dummyobs     = prior_dummyobs
)

PRIOR_LABELS <- c(
  minnesota     = "Minnesota (Litterman 1986)",
  normalwishart = "Normal-Wishart conjugate",
  dummyobs      = "Dummy observations (BGR 2010)"
)

# ── 2. Estimate all models ────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 2: ESTIMATE ALL MODELS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("3 priors × 3 lag orders × 3 samples = up to 27 models\n\n")

ml_opt <- bv_metropolis(
  scale_hess = 0.01,
  adjust_acc = TRUE,
  acc_lower  = 0.25,
  acc_upper  = 0.45
)

LAG_ORDERS <- c(1, 2, 4)

bvar_results  <- list()   # bvar_results[[prior]][[sample]][[paste0("p",p)]]
bvar_summary_rows <- list()

for (prior_name in names(PRIORS)) {
  bvar_results[[prior_name]] <- list()
  cat(sprintf("╔══════════════════════════════════════════════════════════╗\n"))
  cat(sprintf("║  Prior: %-50s║\n", PRIOR_LABELS[prior_name]))
  cat(sprintf("╚══════════════════════════════════════════════════════════╝\n"))
  
  for (sname in names(samples)) {
    s <- samples[[sname]]
    bvar_results[[prior_name]][[sname]] <- list()
    
    dat <- var_adj %>%
      filter(date >= s$start_date) %>%
      select(date,
             gdp_g    = gdp_g_adj,
             cpi_g    = cpi_g_adj,
             bond_dif = bond_dif_adj) %>%
      arrange(date)
    n <- nrow(dat)
    Y <- as.matrix(dat[, c("gdp_g", "cpi_g", "bond_dif")])
    
    cat(sprintf("\n  [%s | n=%d]\n", sname, n))
    
    for (p in LAG_ORDERS) {
      pkey <- paste0("p", p)
      id   <- sprintf("bvar_%s_%s_p%d", prior_name, sname, p)
      
      if (n < p * 3 + 20) {
        cat(sprintf("    p=%d SKIPPED (n=%d insufficient)\n", p, n))
        bvar_results[[prior_name]][[sname]][[pkey]] <- NULL
        next
      }
      
      cat(sprintf("    p=%d estimating... ", p))
      
      mod <- tryCatch(
        bvar(data   = Y, lags   = p,
             n_draw = N_DRAW, n_burn = N_BURN, n_thin = N_THIN,
             priors = PRIORS[[prior_name]], mh = ml_opt,
             verbose = FALSE),
        error = function(e) { message("ERROR: ", e$message); NULL }
      )
      
      if (is.null(mod)) {
        bvar_results[[prior_name]][[sname]][[pkey]] <- NULL
        next
      }
      
      fc <- tryCatch(
        predict(mod, horizon = N_AHEAD,
                conf_bands = c(0.05, 0.16, 0.84, 0.95)),
        error = function(e) NULL
      )
      
      # Posterior mean Q1 2026 forecasts
      fc_gdp  <- if (!is.null(fc)) mean(fc$fcast[, 1, 1]) else NA_real_
      fc_cpi  <- if (!is.null(fc)) mean(fc$fcast[, 1, 2]) else NA_real_
      fc_bond <- if (!is.null(fc)) mean(fc$fcast[, 1, 3]) else NA_real_
      
      # Lambda posterior
      lambda_post <- tryCatch({
        h <- mod$hyper
        if ("lambda" %in% colnames(h)) mean(h[, "lambda"])
        else mean(h[, 1])
      }, error = function(e) NA_real_)
      
      # SOC and SUR posteriors (dummy obs prior only)
      soc_post <- tryCatch(mean(mod$hyper[, "soc"]), error = function(e) NA_real_)
      sur_post <- tryCatch(mean(mod$hyper[, "sur"]), error = function(e) NA_real_)
      
      # Posterior mean sigma GDP
      sigma_gdp <- tryCatch(
        mean(sqrt(mod$sigma[, 1, 1])),
        error = function(e) NA_real_
      )
      
      bvar_results[[prior_name]][[sname]][[pkey]] <- list(
        id          = id,
        prior       = prior_name,
        sample      = sname,
        p           = p,
        n_obs       = n,
        model       = mod,
        fc          = fc,
        lambda_post = lambda_post,
        soc_post    = soc_post,
        sur_post    = sur_post,
        sigma_gdp   = sigma_gdp
      )
      
      bvar_summary_rows[[length(bvar_summary_rows) + 1]] <- tibble(
        id          = id,
        prior       = prior_name,
        sample      = sname,
        p           = p,
        n_obs       = n,
        lambda      = round(lambda_post, 3),
        soc         = round(soc_post, 3),
        sur         = round(sur_post, 3),
        sigma_gdp   = round(sigma_gdp, 4),
        fc_gdp_q1   = round(fc_gdp,   3),
        fc_cpi_q1   = round(fc_cpi,   3),
        fc_bond_q1  = round(fc_bond,  3)
      )
      
      # Compact print
      hyper_str <- sprintf("lambda=%.3f", lambda_post)
      if (!is.na(soc_post)) hyper_str <- paste0(hyper_str,
                                                sprintf(" soc=%.3f sur=%.3f",
                                                        soc_post, sur_post))
      cat(sprintf("done | %s | σ=%.4f | GDP=%+.3f%%\n",
                  hyper_str, sigma_gdp, fc_gdp))
    }
  }
  cat("\n")
}

bvar_summary <- bind_rows(bvar_summary_rows)
n_estimated  <- nrow(bvar_summary)
cat(sprintf("Estimation complete: %d models\n\n", n_estimated))

# ── 3. Summary table and prior comparison ─────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 3: SUMMARY AND PRIOR COMPARISON\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

cat("── Full summary table ───────────────────────────────────────────────\n")
print(bvar_summary, n = Inf, width = Inf)

cat("\n── Best model per prior (lowest GDP sigma) ──────────────────────────\n")
best_per_prior <- bvar_summary %>%
  group_by(prior) %>%
  slice_min(sigma_gdp, n = 1) %>%
  ungroup() %>%
  select(prior, sample, p, lambda, soc, sur,
         sigma_gdp, fc_gdp_q1, fc_cpi_q1, fc_bond_q1)
print(best_per_prior, width = Inf)

cat("\n── Prior comparison: how different are the forecasts? ───────────────\n")
cat("Same sample and lag, different prior:\n\n")
bvar_summary %>%
  filter(sample == "post08", p == 1) %>%
  select(prior, lambda, soc, sur, sigma_gdp,
         fc_gdp_q1, fc_cpi_q1, fc_bond_q1) %>%
  print(width = Inf)

cat("\n── Lambda by prior: does conjugacy / dummy obs change shrinkage? ────\n")
bvar_summary %>%
  group_by(prior) %>%
  summarise(
    mean_lambda = round(mean(lambda, na.rm = TRUE), 3),
    sd_lambda   = round(sd(lambda,   na.rm = TRUE), 3),
    mean_sigma  = round(mean(sigma_gdp, na.rm = TRUE), 4),
    mean_gdp_q1 = round(mean(fc_gdp_q1, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print(width = Inf)

# ── 4. Fan chart plots: best model per prior ──────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 4: FORECAST FAN CHARTS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Helper: fan chart for one variable from a BVAR result object
plot_bvar_fan <- function(res, var_idx, var_label, col,
                          nowcast_ref = NULL, nowcast_label = NULL) {
  if (is.null(res$fc)) return(NULL)
  s  <- samples[[res$sample]]
  fc <- res$fc
  
  hist_col <- c("gdp_g_adj", "cpi_g_adj", "bond_dif_adj")[var_idx]
  hist_dat <- var_adj %>%
    filter(date >= s$start_date) %>%
    tail(16) %>%
    select(date, value = !!sym(hist_col))
  
  last_date <- max(hist_dat$date)
  fc_dates  <- seq(last_date + months(3), by = "quarter", length.out = N_AHEAD)
  
  fc_draws  <- fc$fcast[, , var_idx]
  fc_q      <- apply(fc_draws, 2, quantile,
                     probs = c(0.05, 0.16, 0.50, 0.84, 0.95))
  
  fc_df <- tibble(date = fc_dates,
                  q05 = fc_q[1,], q16 = fc_q[2,],
                  med = fc_q[3,], q84 = fc_q[4,], q95 = fc_q[5,])
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = COL_GREY, linewidth = 0.3) +
    geom_ribbon(data = fc_df, aes(x = date, ymin = q05, ymax = q95),
                fill = col, alpha = 0.10, inherit.aes = FALSE) +
    geom_ribbon(data = fc_df, aes(x = date, ymin = q16, ymax = q84),
                fill = col, alpha = 0.20, inherit.aes = FALSE) +
    geom_line(data = hist_dat, aes(x = date, y = value),
              colour = col, linewidth = 0.8) +
    geom_point(data = hist_dat, aes(x = date, y = value),
               colour = col, size = 1.5) +
    geom_line(data = fc_df, aes(x = date, y = med),
              colour = col, linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = fc_df$date[1], y = fc_df$med[1] + 0.08,
             label = sprintf("%+.2f", fc_df$med[1]),
             size = 3.2, colour = col, fontface = "bold") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(title = var_label, x = NULL, y = NULL)
  
  if (!is.null(nowcast_ref)) {
    p <- p +
      geom_hline(yintercept = nowcast_ref, linetype = "dotted",
                 color = "gray50", linewidth = 0.4) +
      annotate("text", x = min(fc_df$date), y = nowcast_ref + 0.08,
               label = nowcast_label, size = 2.5, color = "gray40", hjust = 0)
  }
  p
}

# Plot best model per prior (3 panels each)
for (prior_name in names(PRIORS)) {
  best_row <- best_per_prior %>% filter(prior == prior_name)
  res <- bvar_results[[prior_name]][[best_row$sample]][[paste0("p", best_row$p)]]
  if (is.null(res)) next
  
  hyper_str <- sprintf("λ=%.3f", res$lambda_post)
  if (!is.na(res$soc_post))
    hyper_str <- paste0(hyper_str,
                        sprintf(" soc=%.3f sur=%.3f",
                                res$soc_post, res$sur_post))
  
  p_gdp  <- plot_bvar_fan(res, 1, "GDP growth (QoQ %)",
                          COL_GDP, 0.30, "Nowcasting Lab: +0.30%")
  p_cpi  <- plot_bvar_fan(res, 2, "CPI inflation (QoQ %)",       COL_CPI)
  p_bond <- plot_bvar_fan(res, 3, "Bond yield change (QoQ pp)",  COL_BOND)
  
  fig <- wrap_plots(Filter(Negate(is.null), list(p_gdp, p_cpi, p_bond)),
                    ncol = 3) +
    plot_annotation(
      title    = sprintf("BVAR forecast — %s (%s)",
                         res$id, PRIOR_LABELS[prior_name]),
      subtitle = sprintf("%s | p=%d | %s | 68%% & 90%% posterior CI",
                         samples[[res$sample]]$label, res$p, hyper_str),
      theme = theme(plot.title    = element_text(size = 11, face = "bold"),
                    plot.subtitle = element_text(size = 8,  color = "gray40"))
    )
  print(fig)
  
  # # Save (uncomment)
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("bvar_01_forecast_%s.png", prior_name)),
  #        fig, width = 14, height = 5, dpi = 150)
}

# ── 5. IRF with posterior bands: best model per prior ─────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 5: IRF WITH POSTERIOR UNCERTAINTY BANDS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Key advantage of BVAR: proper posterior bands on IRFs\n")
cat("(frequentist boot=TRUE failed with exogen — solved here)\n\n")

IRF_HORIZON <- 12
VAR_LABELS_B <- c("GDP growth", "CPI inflation", "Bond yield Δ")

for (prior_name in names(PRIORS)) {
  best_row <- best_per_prior %>% filter(prior == prior_name)
  res      <- bvar_results[[prior_name]][[best_row$sample]][[paste0("p", best_row$p)]]
  if (is.null(res)) next
  
  cat(sprintf("── IRF: %s ──────────────────────────────────────\n", res$id))
  
  bvar_irf <- tryCatch(
    irf(res$model, horizon = IRF_HORIZON, conf_bands = c(0.16, 0.84)),
    error = function(e) { message("IRF error: ", e$message); NULL }
  )
  if (is.null(bvar_irf)) next
  
  # Built-in BVAR plot (clean 3×3 grid with posterior bands)
  plot(bvar_irf, area = TRUE,
       col  = c(COL_GDP, COL_CPI, COL_BOND),
       mar  = c(2, 2, 2, 0.5))
  title(main = sprintf("BVAR IRF — %s | %s | 68%% posterior",
                       res$id, PRIOR_LABELS[prior_name]),
        outer = TRUE, line = -1, cex.main = 0.9)
  
  # Print key cross-variable responses
  tryCatch({
    irf_med <- apply(bvar_irf$irf, c(1, 2, 3), median)
    vn <- c("GDP", "CPI", "Bond")
    cat("  Posterior median IRF at h=1,2,4,8:\n")
    for (imp in 1:3) {
      for (resp in 1:3) {
        if (imp == resp) next
        v <- irf_med[c(1,2,4,8), resp, imp]
        cat(sprintf("    %s → %s: %+.4f  %+.4f  %+.4f  %+.4f\n",
                    vn[imp], vn[resp], v[1], v[2], v[3], v[4]))
      }
    }
  }, error = function(e) message("Could not extract IRF: ", e$message))
  cat("\n")
}

cat("Script complete.\n")
cat(sprintf("Objects: bvar_results (%d models), bvar_summary, best_per_prior\n",
            n_estimated))
cat("Next: bvar/02_bvar_baseline_evaluation.R\n")