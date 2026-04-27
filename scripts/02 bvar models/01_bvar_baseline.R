# ============================================================================
# scripts/models/bvar/01_bvar_baseline.R
#
# PURPOSE:
#   Estimate Bayesian VAR (BVAR) baseline models with Minnesota prior for the
#   same 3-variable system as the frequentist baseline: GDP growth, CPI
#   inflation, and bond yield change. Compare to frequentist VAR and AR(1).
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/var/01_var_baseline.R")  — var_input, samples,
#                                                      exog_full, covid_cols
#   install.packages("BVAR")
#
# PACKAGE:
#   BVAR (Kuschnig & Vashold 2021) — Minnesota prior, Metropolis-Hastings
#   sampler, built-in predict/irf/fevd methods.
#   Reference: Kuschnig N, Vashold L (2021). "BVAR: Bayesian Vector
#   Autoregressions with Hierarchical Prior Selection in R."
#   Journal of Statistical Software.
#
# MINNESOTA PRIOR:
#   The Minnesota prior (Litterman 1986) shrinks VAR coefficients toward:
#     - Own first lag: coefficient = 1 (random walk prior per variable)
#     - Other lags / cross-variable lags: coefficient = 0 (shrinkage)
#   Hyperparameters:
#     lambda (overall tightness): how strongly to shrink toward prior
#       → small lambda = strong shrinkage → closer to random walk
#       → large lambda = weak shrinkage → closer to OLS VAR
#     alpha  (lag decay): how fast shrinkage increases with lag length
#     psi    (cross-variable shrinkage): relative tightness on cross-variable
#   We use hierarchical prior on lambda (optimize over the marginal likelihood)
#   which lets the data choose the degree of shrinkage automatically.
#
# MODEL GRID:
#   Same structure as frequentist baseline:
#   Lag orders: p = 1, 2, 4
#   Samples   : full (Q2 2000–Q4 2025), post08 (Q1 2009–Q4 2025),
#               post15 (Q1 2015–Q4 2025)
#   Exogenous : 8 COVID dummies (note: BVAR package handles exog differently —
#               see Section 0 for approach)
#
# COVID TREATMENT IN BVAR:
#   The BVAR package does not support exogenous variables natively in the
#   same way as vars::VAR(). We handle COVID via two approaches:
#   Approach A: Include COVID dummies as additional endogenous variables
#               (not ideal — adds unnecessary equations)
#   Approach B: Demean/adjust the data for COVID outliers before estimation
#               (cleaner — remove COVID effect from the series directly)
#   Approach C: Use dummy observations prior (Banbura et al. 2010)
#               (most principled — treats COVID as prior information)
#   → We use Approach B: adjust GDP/CPI/bond for COVID quarters using
#     the residuals from the frequentist VAR COVID dummies, then estimate
#     BVAR on the adjusted series.
#
# MCMC SETTINGS:
#   n_draw   = 10000  (posterior draws)
#   n_burn   = 5000   (burn-in draws, discarded)
#   n_thin   = 1      (thinning — keep every draw)
#   These settings give reliable posterior estimates for a 3-variable system.
#   Increase n_draw to 25000 for final paper results if time permits.
#
# STRUCTURE:
#   Section 0 — Install/load BVAR, COVID adjustment
#   Section 1 — Estimate BVAR models (p × sample grid)
#   Section 2 — Posterior summaries and coefficient inspection
#   Section 3 — Forecasts with fan charts (Q1 2026 nowcast)
#   Section 4 — Compare to frequentist VAR and AR(1) (in-sample)
#   Section 5 — IRF with posterior uncertainty bands
#
# KEY FINDINGS:
#
#   ── Lambda (shrinkage hyperparameter) ────────────────────────────────────
#   All models: moderate to weak shrinkage (lambda 0.46–0.72)
#   → Data prefer less shrinkage than the prior mode (0.2)
#   → Suggests Swiss macro dynamics are sufficiently informative;
#     pure random walk prior is too restrictive
#   full:   p1=0.636, p2=0.515, p4=0.461  (shrinkage increases with lags ✓)
#   post08: p1=0.717, p2=NA,    p4=0.531
#   post15: p1=0.675, p2=0.642, p4=skipped
#
#   ── Q1 2026 forecasts (posterior medians) ────────────────────────────────
#   GDP  range: +0.28% to +0.54% across models (all above Nowcasting +0.30%)
#   CPI  range: -0.05% to +0.04% (near-zero — consistent with extended VAR)
#   Bond range: -0.13pp to +0.01pp (slight yield decline expected)
#   Note: post15 models show wider GDP uncertainty (small sample)
#
#   ── BVAR vs frequentist VAR (in-sample sigma) ────────────────────────────
#   (To be filled after Section 4 runs with results available)
#   Key comparison: bvar_post08_p1 vs VAR post08_p3
#   BVAR shrinkage trades off in-sample fit for better out-of-sample stability
#
#   ── IRF with posterior bands (bvar_post08_p1) ────────────────────────────
#   GDP  → CPI : positive, peak h=1-2, 68% band excludes 0 ✓
#   Bond → CPI : +0.04 to +0.08pp, hump-shaped, band mostly positive ✓
#   CPI  → Bond: large positive response ~0.3pp at h=1, decays — strong ✓
#   GDP  → Bond: negative contemporaneous response (counterintuitive —
#               consistent with frequentist finding; safe-haven channel)
#   Key advantage vs frequentist: proper posterior bands on all IRFs ✓
#   (frequentist IRF bands failed with exogen — BVAR solves this)
#
# Outputs (commented — uncomment to save):
#   output/figures/bvar_01_forecast_<model>.png
#   output/figures/bvar_01_irf_<model>.png
#   output/tables/bvar_01_summary.csv
# ============================================================================

# ── 0. Setup ──────────────────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 0: SETUP — BVAR PACKAGE AND COVID ADJUSTMENT\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

if (!requireNamespace("BVAR", quietly = TRUE)) {
  install.packages("BVAR")
}
library(BVAR)
library(dplyr)

cat("BVAR package version:", as.character(packageVersion("BVAR")), "\n\n")

# ── COVID adjustment (Approach B) ─────────────────────────────────────────────
# Use the coefficients of the COVID dummies from the frequentist VAR to
# remove the COVID effect from the raw series. This gives us a
# "COVID-adjusted" series that can be fed directly into the BVAR.
#
# Step 1: Fit a simple OLS regression of each target on COVID dummies
# Step 2: Subtract fitted COVID effect from the raw series
# Step 3: Use adjusted series for BVAR estimation

cat("── COVID adjustment ──────────────────────────────────────────────────\n")
cat("Removing COVID dummy effects from series before BVAR estimation\n\n")

TARGET_VARS <- c("gdp_g", "cpi_g", "bond_dif")

covid_adjust <- function(y, exog_mat) {
  # Simple OLS: y ~ COVID dummies
  # Align lengths — exog_mat may have different number of rows than y
  n    <- length(y)
  df   <- as.data.frame(exog_mat[1:n, , drop = FALSE])
  df$y <- y
  fit  <- lm(y ~ . - y, data = df)
  fitted_vals <- predict(fit)
  # Subtract COVID effect, preserve mean level
  y - fitted_vals + mean(fitted_vals)
}

# Build COVID-adjusted var_input
var_adj <- var_input

for (vname in TARGET_VARS) {
  adj <- covid_adjust(var_input[[vname]], exog_full)
  var_adj[[paste0(vname, "_adj")]] <- adj
  cat(sprintf("  %s: COVID adjustment applied | max change: %.4f pp\n",
              vname,
              max(abs(var_input[[vname]] - adj))))
}

cat("\nAdjusted columns created:", paste(paste0(TARGET_VARS, "_adj"), collapse=", "), "\n\n")

# ── MCMC settings ─────────────────────────────────────────────────────────────
N_DRAW <- 10000
N_BURN <- 5000
N_THIN <- 1
set.seed(2026)

cat(sprintf("MCMC: %d draws, %d burn-in, thinning=%d\n",
            N_DRAW, N_BURN, N_THIN))
cat("Note: this may take several minutes per model\n\n")

# ── 1. Estimate BVAR models ───────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 1: ESTIMATE BVAR MODELS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Minnesota prior setup:
# - lambda: hierarchical (data-driven optimal tightness)
# - alpha:  2 (standard quadratic lag decay)
# - psi:    diagonal of residual variance (standard Minnesota)
mn_prior <- bv_priors(
  hyper  = "auto",
  mn     = bv_minnesota(
    lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 1e-4, max = 5),
    alpha  = bv_alpha(mode = 2),
    psi    = bv_psi()
  )
)

# Marginal likelihood for hyperparameter optimization
ml_opt <- bv_metropolis(
  scale_hess = 0.01,
  adjust_acc = TRUE,
  acc_lower  = 0.25,
  acc_upper  = 0.45
)

LAG_ORDERS <- c(1, 2, 4)
N_AHEAD    <- 8

bvar_results <- list()
bvar_summary_rows <- list()

for (sname in names(samples)) {
  s      <- samples[[sname]]
  bvar_results[[sname]] <- list()
  
  # Extract COVID-adjusted data for this sample
  dat <- var_adj %>%
    filter(date >= s$start_date) %>%
    select(date,
           gdp_g    = gdp_g_adj,
           cpi_g    = cpi_g_adj,
           bond_dif = bond_dif_adj) %>%
    arrange(date)
  
  n <- nrow(dat)
  cat(sprintf("── Sample: %s [n=%d, %s – %s] ─────────────────────────────\n",
              sname, n,
              format(min(dat$date)), format(max(dat$date))))
  
  for (p in LAG_ORDERS) {
    id <- sprintf("bvar_%s_p%d", sname, p)
    
    # Min obs check
    if (n < p * 3 + 20) {
      cat(sprintf("  [p=%d] SKIPPED — insufficient obs (n=%d)\n", p, n))
      bvar_results[[sname]][[paste0("p", p)]] <- NULL
      next
    }
    
    cat(sprintf("  [p=%d] Estimating BVAR... ", p))
    
    # Convert to matrix for BVAR package
    Y <- as.matrix(dat[, c("gdp_g","cpi_g","bond_dif")])
    
    mod <- tryCatch({
      bvar(
        data     = Y,
        lags     = p,
        n_draw   = N_DRAW,
        n_burn   = N_BURN,
        n_thin   = N_THIN,
        priors   = mn_prior,
        mh       = ml_opt,
        verbose  = FALSE
      )
    }, error = function(e) {
      message("ERROR: ", e$message)
      NULL
    })
    
    if (is.null(mod)) {
      bvar_results[[sname]][[paste0("p", p)]] <- NULL
      next
    }
    
    # Forecast
    fc <- tryCatch(
      predict(mod, horizon = N_AHEAD, conf_bands = c(0.16, 0.84, 0.05, 0.95)),
      error = function(e) NULL
    )
    
    # Extract posterior mean and CI for Q1 2026 (h=1)
    fc_gdp_mean <- if (!is.null(fc)) mean(fc$fcast[, 1, 1]) else NA
    fc_cpi_mean <- if (!is.null(fc)) mean(fc$fcast[, 1, 2]) else NA
    fc_bond_mean<- if (!is.null(fc)) mean(fc$fcast[, 1, 3]) else NA
    
    # Lambda posterior (optimal shrinkage)
    lambda_post <- tryCatch({
      # Column name varies by BVAR version — try both
      h <- mod$hyper
      if ("lambda" %in% colnames(h)) mean(h[, "lambda"])
      else if (ncol(h) >= 1) mean(h[, 1])
      else NA_real_
    }, error = function(e) NA_real_)
    
    # In-sample fit: posterior mean residual SD for GDP
    sigma_gdp <- tryCatch({
      fitted_mean <- apply(mod$beta, c(2,3), mean)  # mean coefficients
      # approximate: use posterior mean of sigma
      mean(sqrt(mod$sigma[, 1, 1]))
    }, error = function(e) NA_real_)
    
    bvar_results[[sname]][[paste0("p", p)]] <- list(
      id          = id,
      sample      = sname,
      p           = p,
      n_obs       = n,
      model       = mod,
      fc          = fc,
      lambda_post = lambda_post,
      sigma_gdp   = sigma_gdp
    )
    
    bvar_summary_rows[[length(bvar_summary_rows) + 1]] <- tibble(
      id           = id,
      sample       = sname,
      p            = p,
      n_obs        = n,
      lambda_post  = round(lambda_post, 4),
      sigma_gdp    = round(sigma_gdp, 4),
      fc_gdp_q1    = round(fc_gdp_mean, 3),
      fc_cpi_q1    = round(fc_cpi_mean, 3),
      fc_bond_q1   = round(fc_bond_mean, 3)
    )
    
    cat(sprintf("done | lambda=%.3f | σ(GDP)=%.4f | Q1 GDP=%+.3f%%\n",
                lambda_post, sigma_gdp, fc_gdp_mean))
  }
  cat("\n")
}

bvar_summary <- bind_rows(bvar_summary_rows)

# ── 2. Posterior summaries ────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 2: POSTERIOR SUMMARIES\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

cat("── Model summary table ──────────────────────────────────────────────\n")
print(bvar_summary, n = Inf, width = Inf)

cat("\n── Lambda interpretation ─────────────────────────────────────────────\n")
cat("lambda = posterior mean of Minnesota tightness hyperparameter\n")
cat("  Small lambda (<0.1) : strong shrinkage → data support random walk prior\n")
cat("  Large lambda (>1.0) : weak shrinkage → data support unrestricted VAR\n\n")

for (sname in names(samples)) {
  for (p in LAG_ORDERS) {
    res <- bvar_results[[sname]][[paste0("p", p)]]
    if (is.null(res)) next
    
    cat(sprintf("  %s p=%d: lambda=%.3f → %s\n",
                sname, p, res$lambda_post,
                ifelse(res$lambda_post < 0.2, "strong shrinkage",
                       ifelse(res$lambda_post < 0.5, "moderate shrinkage",
                              "weak shrinkage"))))
  }
}

# ── 3. Forecasts with fan charts ──────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 3: FORECASTS WITH FAN CHARTS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Helper: plot fan chart for one model, one variable
plot_bvar_fan <- function(res, var_idx, var_name, var_label, col,
                          nowcast_ref = NULL, nowcast_label = NULL) {
  if (is.null(res$fc)) return(NULL)
  
  s    <- samples[[res$sample]]
  fc   <- res$fc
  
  # Historical data (last 16 quarters)
  hist_col <- switch(var_name,
                     gdp_g    = "gdp_g_adj",
                     cpi_g    = "cpi_g_adj",
                     bond_dif = "bond_dif_adj")
  
  hist_dat <- var_adj %>%
    filter(date >= s$start_date) %>%
    tail(16) %>%
    select(date, value = !!sym(hist_col))
  
  last_date <- max(hist_dat$date)
  fc_dates  <- seq(last_date + months(3),
                   by = "quarter", length.out = N_AHEAD)
  
  # Extract quantiles from BVAR forecast object
  # fc$fcast: array [n_draw, horizon, variable]
  fc_draws <- fc$fcast[, , var_idx]  # [n_draw, horizon]
  fc_q     <- apply(fc_draws, 2, quantile,
                    probs = c(0.05, 0.16, 0.50, 0.84, 0.95))
  
  fc_df <- tibble(
    date  = fc_dates,
    q05   = fc_q[1, ],
    q16   = fc_q[2, ],
    med   = fc_q[3, ],
    q84   = fc_q[4, ],
    q95   = fc_q[5, ]
  )
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = COL_GREY, linewidth = 0.3) +
    # 90% CI (outer band)
    geom_ribbon(data = fc_df, aes(x = date, ymin = q05, ymax = q95),
                fill = col, alpha = 0.10, inherit.aes = FALSE) +
    # 68% CI (inner band — BoE style)
    geom_ribbon(data = fc_df, aes(x = date, ymin = q16, ymax = q84),
                fill = col, alpha = 0.20, inherit.aes = FALSE) +
    # Historical
    geom_line(data = hist_dat, aes(x = date, y = value),
              colour = col, linewidth = 0.8) +
    geom_point(data = hist_dat, aes(x = date, y = value),
               colour = col, size = 1.5) +
    # Median forecast
    geom_line(data = fc_df, aes(x = date, y = med),
              colour = col, linewidth = 0.8, linetype = "dashed") +
    # Point estimate label
    annotate("text", x = fc_df$date[1], y = fc_df$med[1] + 0.08,
             label = sprintf("%+.2f", fc_df$med[1]),
             size = 3.2, colour = col, fontface = "bold") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(title = var_label, x = NULL, y = NULL)
  
  # Reference line (e.g. Nowcasting Lab)
  if (!is.null(nowcast_ref)) {
    p <- p +
      geom_hline(yintercept = nowcast_ref, linetype = "dotted",
                 color = "gray50", linewidth = 0.4) +
      annotate("text", x = min(fc_df$date), y = nowcast_ref + 0.08,
               label = nowcast_label, size = 2.5, color = "gray40", hjust = 0)
  }
  p
}

# Plot 3-panel forecast for each model
for (sname in names(samples)) {
  for (p in LAG_ORDERS) {
    res <- bvar_results[[sname]][[paste0("p", p)]]
    if (is.null(res) || is.null(res$fc)) next
    
    p_gdp  <- plot_bvar_fan(res, 1, "gdp_g",    "GDP growth (QoQ %)",
                            COL_GDP, 0.30, "Nowcasting Lab: +0.30%")
    p_cpi  <- plot_bvar_fan(res, 2, "cpi_g",    "CPI inflation (QoQ %)",  COL_CPI)
    p_bond <- plot_bvar_fan(res, 3, "bond_dif", "Bond yield change (QoQ pp)", COL_BOND)
    
    panels <- Filter(Negate(is.null), list(p_gdp, p_cpi, p_bond))
    if (length(panels) == 0) next
    
    fig <- wrap_plots(panels, ncol = 3) +
      plot_annotation(
        title    = sprintf("BVAR forecast — %s", res$id),
        subtitle = sprintf("%s | p=%d | Minnesota prior (λ=%.3f) | 68%% & 90%% posterior CI",
                           samples[[sname]]$label, p, res$lambda_post),
        theme = theme(
          plot.title    = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9,  color = "gray40")
        )
      )
    print(fig)
    
    # # Save (uncomment)
    # ggsave(file.path(here("output","figures"),
    #                  sprintf("bvar_01_forecast_%s.png", res$id)),
    #        fig, width = 14, height = 5, dpi = 150)
  }
}

# ── 4. Compare BVAR vs frequentist VAR ───────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 4: BVAR vs FREQUENTIST VAR — IN-SAMPLE COMPARISON\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Comparing posterior mean residual SD vs frequentist sigma\n\n")

# Pull frequentist sigmas from baseline results (from 01_var_baseline.R)
freq_sigma <- tryCatch({
  bind_rows(lapply(names(samples), function(sname) {
    lapply(c("p1","p2","p4"), function(pk) {
      res <- results[[sname]][[pk]]
      if (is.null(res)) return(NULL)
      tibble(
        sample     = sname,
        p          = res$p,
        freq_sigma = round(res$diag$sigma_gdp, 4)
      )
    }) %>% bind_rows()
  }))
}, error = function(e) NULL)

comp_tbl <- bvar_summary %>%
  select(sample, p, lambda_post, bvar_sigma = sigma_gdp,
         fc_gdp_q1, fc_cpi_q1, fc_bond_q1)

if (!is.null(freq_sigma)) {
  comp_tbl <- comp_tbl %>%
    left_join(freq_sigma, by = c("sample","p")) %>%
    mutate(sigma_improvement = round((freq_sigma - bvar_sigma) / freq_sigma * 100, 1))
  
  cat("── Sigma comparison (GDP equation) ──────────────────────────────────\n")
  cat("Positive improvement % = BVAR has lower residual SD than frequentist\n\n")
  print(comp_tbl %>%
          select(sample, p, lambda_post, freq_sigma, bvar_sigma, sigma_improvement),
        n = Inf, width = Inf)
} else {
  cat("Frequentist results not available — run 01_var_baseline.R first\n")
  print(comp_tbl, n = Inf, width = Inf)
}

# ── 5. IRF with posterior uncertainty ─────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 5: IMPULSE RESPONSE FUNCTIONS WITH POSTERIOR BANDS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Key advantage of BVAR: proper posterior uncertainty bands on IRFs\n")
cat("(vs asymptotic CI from frequentist VAR which failed with exogen)\n\n")

VAR_LABELS_B <- c("GDP growth", "CPI inflation", "Bond yield Δ")
IRF_HORIZON  <- 12

# Run IRF for best model (post08_p1 or post08_p3 equivalent)
best_sname <- "post08"
best_p     <- 1
best_res   <- bvar_results[[best_sname]][[paste0("p", best_p)]]

if (!is.null(best_res)) {
  cat(sprintf("IRF for: %s\n\n", best_res$id))
  
  bvar_irf <- tryCatch(
    irf(best_res$model, horizon = IRF_HORIZON,
        conf_bands = c(0.16, 0.84)),
    error = function(e) { message("IRF error: ", e$message); NULL }
  )
  
  if (!is.null(bvar_irf)) {
    # Plot IRF matrix — use built-in BVAR plot method
    # rows = response, cols = impulse (BVAR convention)
    cat("Plotting IRF (using BVAR built-in plot)...\n")
    plot(bvar_irf, area = TRUE,
         col = c(COL_GDP, COL_CPI, COL_BOND),
         mar = c(2, 2, 2, 0.5))
    title(main = sprintf("BVAR IRF — %s | 68%% posterior bands", best_res$id),
          outer = TRUE, line = -1)
    
    # Also extract and print key values
    cat("\nKey IRF values at h=1,2,4,8 (posterior median):\n")
    tryCatch({
      # BVAR irf object: irf$irf[horizon, response, impulse, draw]
      irf_med <- apply(bvar_irf$irf, c(1,2,3), median)
      var_names <- c("GDP", "CPI", "Bond")
      for (imp in 1:3) {
        for (resp in 1:3) {
          if (imp == resp) next
          vals <- irf_med[c(1,2,4,8), resp, imp]
          cat(sprintf("  %s → %s: h=1:%+.4f h=2:%+.4f h=4:%+.4f h=8:%+.4f\n",
                      var_names[imp], var_names[resp],
                      vals[1], vals[2], vals[3], vals[4]))
        }
      }
    }, error = function(e) message("Could not extract IRF values: ", e$message))
  }
}

cat("\nBVAR baseline estimation complete.\n")
cat(sprintf("Models estimated: %d\n",
            sum(sapply(bvar_results, function(s)
              sum(!sapply(s, is.null))))))
cat("Objects: bvar_results, bvar_summary\n")
cat("Next: bvar/02_bvar_baseline_evaluation.R\n")