# ============================================================================
# scripts/models/bvar/04_bvar_extended_estimation.R
#
# PURPOSE:
#   Estimate extended BVAR models adding candidate predictors to the
#   3-variable baseline (GDP, CPI, bond yield).
#   Mirrors the structure of var/05_var_extended_estimation.R but uses
#   Bayesian estimation via the BVAR package (v1.0.5).
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/bvar/01_bvar_baseline.R") — bvar_results,
#     bvar_summary, best_per_prior, PRIORS, PRIOR_LABELS,
#     var_adj, samples, BASE_COLS
#   source("scripts/models/var/05_var_extended_estimation.R") — ext_data
#     (provides the aligned extended endogenous variables)
#
# DESIGN CHOICES:
#   Prior: Minnesota only (best baseline performance, fastest to estimate)
#   Lags:  p = 1, 4 (p=2 dropped — adds little over p=1 in baseline)
#   Samples: full, post08, post15
#   COVID treatment: Approach B (OLS demean, *_adj series) — same as baseline
#   Draws: 5000 burn + 5000 save (same as baseline)
#   Oil: included as EXOGENOUS in all oil-inclusive sets (same as VARX)
#     → passed via exogen= alongside COVID dummies
#
# VARIABLE SETS (endogenous, mirrors var/05_var_extended_estimation.R):
#   Focus on the sets that performed best in the frequentist VARX:
#
#   4-variable (baseline + 1):
#     bset_ea    : + ea_esi        (EA ESI — best for bond)
#     bset_world : + msci_wld      (MSCI World — global demand)
#     bset_sw    : + kof_bar       (KOF barometer — best for GDP)
#
#   5-variable (baseline + 2):
#     bset_ext   : + fx_eur + ea_esi          (core external)
#     bset_dom   : + kof_bar + imp_prc        (domestic + prices)
#
#   6-variable (baseline + 3):
#     bset_medium: + fx_eur + ea_esi + kof_bar
#
#   With oil exogenous:
#     bset_oil      : baseline only + oil exog
#     bset_oil_sw   : + kof_bar    [oil exog]
#     bset_oil_ext  : + fx_eur + ea_esi  [oil exog]
#     bset_oil_med  : + fx_eur + ea_esi + kof_bar  [oil exog]
#
#   Total: 10 sets × 2 lags × 3 samples = 60 models
#   (many skipped for df constraints)
#
# NOTE ON FORECAST DIMENSIONS (verified):
#   fc$fcast: [draw, horizon, k] where k = n endogenous variables
#   predict.bvar does NOT accept exogen= or horizon= with newdata=
#   Pass only: predict(mod, horizon = N_AHEAD, conf_bands = c(0.16, 0.84))
#   newdata= replaces conditioning history, not future exog values
#
# EVALUATION:
#   In-sample: sigma(GDP), sigma(CPI), sigma(bond) from posterior mean residuals
#   Out-of-sample: expanding-window RMSE at h=1,2,4,8 (2000 draws)
#   Forecast: Q1 2026 point forecast + 68% posterior interval
#
# KEY FINDINGS:
#
#   In-sample (sigma, full sample, p=4):
#     Best overall: bset_ea_full_p4 — σ(GDP)=1.06 vs baseline 1.25 (-15%)
#     Best for CPI/bond: bset_ext/bset_oil_ext — marginal improvement only
#     All best models: full sample, p=4 (consistent with baseline BVAR)
#     Oil exog makes negligible difference at 4-var level
#     Larger models (6-var) do NOT improve over 4-5-var
#
#   Out-of-sample RMSE (expanding window, best model per set):
#     GDP h=1: best = bset_oil_med (2.25) — all models similar, none beat AR(1)
#     GDP h=2: best = bset_world (2.06)
#     GDP h=4: best = bset_world (1.66)
#     GDP h=8: best = bset_ea (1.75)
#     CPI: bset_oil_med best at h=1,2,4 — still well above AR(1)
#     Bond: bset_oil best at h=1,4 — marginal differences across sets
#     Extended BVAR does NOT beat AR(1) for any target at any horizon
#     Extended BVAR slightly worse than baseline BVAR out-of-sample
#     (adding variables hurts OOS despite better in-sample fit — overfitting)
#
#   Q1 2026 GDP forecasts:
#     Range: +0.11% (bset_dom) to +0.68% (bset_world_full_p1)
#     Best models (full_p4): +0.32% to +0.44%
#     All consistent with baseline BVAR and Nowcasting Lab (+0.30%)
#     CPI: near zero across all models (-0.12% to +0.04%)
#     Bond: small negative change (-0.09% to -0.19%) for most models
#
# Outputs (commented — uncomment to save):
#   output/figures/bvar_04_fc_<set>.png
#   output/tables/bvar_04_summary.csv
# ============================================================================

library(dplyr)

N_AHEAD     <- 8
LAG_ORDERS  <- c(1, 4)
PRIOR_NAME  <- "minnesota"
N_DRAW      <- 10000
N_BURN      <- 5000
N_THIN      <- 1

# =============================================================================
# SECTION 0: VARIABLE SETS
# =============================================================================
cat("=== SECTION 0: VARIABLE SETS ===\n\n")

bvar_sets <- list(
  # ── No oil ────────────────────────────────────────────────────────────────
  list(id = "bset_ea",     oil = FALSE, n_endo = 4,
       label = "Baseline + EA ESI",
       extra = "ea_esi"),
  list(id = "bset_world",  oil = FALSE, n_endo = 4,
       label = "Baseline + MSCI World",
       extra = "msci_wld"),
  list(id = "bset_sw",     oil = FALSE, n_endo = 4,
       label = "Baseline + KOF barometer",
       extra = "kof_bar"),
  list(id = "bset_ext",    oil = FALSE, n_endo = 5,
       label = "Baseline + CHF/EUR + EA ESI",
       extra = c("fx_eur", "ea_esi")),
  list(id = "bset_dom",    oil = FALSE, n_endo = 5,
       label = "Baseline + KOF barometer + Import prices",
       extra = c("kof_bar", "imp_prc")),
  list(id = "bset_medium", oil = FALSE, n_endo = 6,
       label = "Baseline + CHF/EUR + EA ESI + KOF barometer",
       extra = c("fx_eur", "ea_esi", "kof_bar")),
  # ── With oil exogenous ────────────────────────────────────────────────────
  list(id = "bset_oil",     oil = TRUE,  n_endo = 3,
       label = "Baseline only + oil exog",
       extra = character(0)),
  list(id = "bset_oil_sw",  oil = TRUE,  n_endo = 4,
       label = "Baseline + KOF barometer  [oil exog]",
       extra = "kof_bar"),
  list(id = "bset_oil_ext", oil = TRUE,  n_endo = 5,
       label = "Baseline + CHF/EUR + EA ESI  [oil exog]",
       extra = c("fx_eur", "ea_esi")),
  list(id = "bset_oil_med", oil = TRUE,  n_endo = 6,
       label = "Baseline + CHF/EUR + EA ESI + KOF  [oil exog]",
       extra = c("fx_eur", "ea_esi", "kof_bar"))
)

for (vs in bvar_sets) {
  cat(sprintf("  %-14s [%d-var]  oil=%-5s  extra: %s\n",
              vs$id, vs$n_endo,
              ifelse(vs$oil, "TRUE", "FALSE"),
              if (length(vs$extra) > 0) paste(vs$extra, collapse = ", ")
              else "(none)"))
}
cat(sprintf("\nLag orders: %s | Samples: full, post08, post15\n\n",
            paste(LAG_ORDERS, collapse = ", ")))

# =============================================================================
# SECTION 1: PREPARE DATA
# =============================================================================
cat("=== SECTION 1: PREPARE DATA ===\n\n")

adj_cols <- c("gdp_g_adj", "cpi_g_adj", "bond_dif_adj")

ext_cols_available <- setdiff(names(ext_data), "date")
cat(sprintf("ext_data columns available: %s\n\n",
            paste(ext_cols_available, collapse = ", ")))

# Combined data: COVID-adjusted base series + extended variables
combined_data <- var_adj %>%
  select(date, all_of(adj_cols)) %>%
  left_join(ext_data, by = "date") %>%
  arrange(date)

cat(sprintf("combined_data: %d rows x %d cols (%s to %s)\n\n",
            nrow(combined_data),
            ncol(combined_data),
            format(min(combined_data$date)),
            format(max(combined_data$date))))

exog_cols_base     <- colnames(exog_base)
exog_cols_with_oil <- colnames(exog_with_oil)

# =============================================================================
# SECTION 2: ESTIMATE EXTENDED BVAR MODELS
# =============================================================================
cat("=== SECTION 2: ESTIMATE EXTENDED BVAR MODELS ===\n\n")
cat(sprintf("Prior: %s | Draws: %d burn + %d save | Thin: %d\n\n",
            PRIOR_NAME, N_BURN, N_DRAW - N_BURN, N_THIN))

bvar_ext_results  <- list()
bvar_ext_rows     <- list()
n_estimated       <- 0
n_skipped         <- 0

for (vs in bvar_sets) {
  bvar_ext_results[[vs$id]] <- list()
  endo_cols_raw <- c(BASE_COLS, vs$extra)
  endo_cols_adj <- c("gdp_g_adj", "cpi_g_adj", "bond_dif_adj", vs$extra)
  k             <- length(endo_cols_raw)
  
  cat(sprintf("\n── %s [%d-var]: %s ──\n", vs$id, k, vs$label))
  
  for (sname in names(samples)) {
    s <- samples[[sname]]
    bvar_ext_results[[vs$id]][[sname]] <- list()
    
    dat <- combined_data %>%
      filter(date >= s$start_date) %>%
      select(date, all_of(endo_cols_adj)) %>%
      arrange(date)
    
    complete_rows <- complete.cases(dat[, endo_cols_adj])
    dat           <- dat[complete_rows, ]
    n_obs         <- nrow(dat)
    
    if (n_obs < 20) {
      cat(sprintf("  [%s] SKIPPED — only %d complete obs\n", sname, n_obs))
      n_skipped <- n_skipped + length(LAG_ORDERS)
      next
    }
    
    Y <- as.matrix(dat[, endo_cols_adj])
    
    exog_mat <- if (vs$oil) exog_with_oil else exog_base
    exog_s   <- exog_mat[var_adj$date %in% dat$date, , drop = FALSE]
    exog_s   <- exog_s[complete_rows, , drop = FALSE]
    
    for (p in LAG_ORDERS) {
      pkey    <- paste0("p", p)
      id      <- sprintf("%s_%s_p%d", vs$id, sname, p)
      min_obs <- k * p + 1 + ncol(exog_s) + 10
      
      if (n_obs < min_obs) {
        cat(sprintf("  [%s p=%d] SKIPPED — n=%d < min_obs=%d\n",
                    sname, p, n_obs, min_obs))
        n_skipped <- n_skipped + 1
        bvar_ext_results[[vs$id]][[sname]][[pkey]] <- NULL
        next
      }
      
      cat(sprintf("  [%s p=%d] n=%d | estimating... ", sname, p, n_obs))
      
      mn_prior <- bv_minnesota(
        lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 1e-4, max = 5),
        alpha  = bv_alpha(mode = 2),
        psi    = bv_psi(),
        var    = 1
      )
      priors <- bv_priors(hyper = "auto", mn = mn_prior)
      
      mod <- tryCatch(
        bvar(
          data    = Y,
          lags    = p,
          n_draw  = N_DRAW,
          n_burn  = N_BURN,
          n_thin  = N_THIN,
          priors  = priors,
          exogen  = exog_s,
          verbose = FALSE
        ),
        error = function(e) {
          cat(sprintf("ERROR: %s\n", e$message))
          NULL
        }
      )
      
      if (is.null(mod)) {
        n_skipped <- n_skipped + 1
        bvar_ext_results[[vs$id]][[sname]][[pkey]] <- NULL
        next
      }
      
      # ── Posterior mean residuals ────────────────────────────────────────
      # mod$beta: [draws, K, M] where K = 1 + k*p (const + lags ONLY)
      # Exogenous variables are NOT in mod$beta in BVAR v1.0.5
      beta_mean <- apply(mod$beta, c(2, 3), mean)  # [K, M]
      K_beta    <- nrow(beta_mean)
      K_expect  <- 1 + k * p
      
      X_list <- lapply((p + 1):n_obs, function(t) {
        c(1, as.vector(t(Y[(t - 1):(t - p), ])))  # const + lags only
      })
      X_mat <- do.call(rbind, X_list)
      Y_mat <- Y[(p + 1):n_obs, ]
      
      if (ncol(X_mat) != K_beta) {
        cat(sprintf("\n  [!] Dimension mismatch: X cols=%d, beta rows=%d — NA sigma\n",
                    ncol(X_mat), K_beta))
        sigma_gdp <- sigma_cpi <- sigma_bond <- NA_real_
      } else {
        resids     <- Y_mat - X_mat %*% beta_mean
        sigma_gdp  <- round(sd(resids[, 1], na.rm = TRUE), 4)
        sigma_cpi  <- round(sd(resids[, 2], na.rm = TRUE), 4)
        sigma_bond <- round(sd(resids[, 3], na.rm = TRUE), 4)
      }
      
      # ── Lambda posterior ────────────────────────────────────────────────
      lambda_mean <- tryCatch({
        h <- mod$hyper
        if ("lambda" %in% colnames(h)) mean(h[, "lambda"]) else mean(h[, 1])
      }, error = function(e) NA_real_)
      
      # ── Acceptance rate ─────────────────────────────────────────────────
      acc_rate <- tryCatch(
        mod$meta$accepted / mod$meta$n_save,
        error = function(e) NA_real_
      )
      
      # ── Forecast ────────────────────────────────────────────────────────
      
      # predict.bvar uses newdata= to replace conditioning history, NOT for
      # future exogenous values. Simply pass horizon + conf_bands.
      # The model conditions on its estimated history automatically.
      fc <- tryCatch(
        predict(mod, horizon = N_AHEAD, conf_bands = c(0.16, 0.84)),
        error = function(e) {
          cat(sprintf("\n  [!] predict() error: %s\n", e$message))
          NULL
        }
      )
      
      # fc$fcast: [draw, horizon, variable] — verified: [5000, N_AHEAD, k]
      if (!is.null(fc) && length(dim(fc$fcast)) == 3) {
        fc_gdp_med  <- apply(fc$fcast[, , 1, drop = FALSE], 2, median)
        fc_gdp_q16  <- apply(fc$fcast[, , 1, drop = FALSE], 2, quantile, 0.16)
        fc_gdp_q84  <- apply(fc$fcast[, , 1, drop = FALSE], 2, quantile, 0.84)
        fc_cpi_med  <- apply(fc$fcast[, , 2, drop = FALSE], 2, median)
        fc_bond_med <- apply(fc$fcast[, , 3, drop = FALSE], 2, median)
      } else if (!is.null(fc) && length(dim(fc$fcast)) == 2) {
        # Some BVAR versions return [horizon, variable] if draws=1
        fc_gdp_med  <- fc$fcast[, 1]
        fc_gdp_q16  <- fc$fcast[, 1]
        fc_gdp_q84  <- fc$fcast[, 1]
        fc_cpi_med  <- fc$fcast[, 2]
        fc_bond_med <- fc$fcast[, 3]
      } else {
        cat(sprintf("\n  [!] predict() returned unexpected structure — NA forecasts\n"))
        fc_gdp_med  <- rep(NA_real_, N_AHEAD)
        fc_gdp_q16  <- rep(NA_real_, N_AHEAD)
        fc_gdp_q84  <- rep(NA_real_, N_AHEAD)
        fc_cpi_med  <- rep(NA_real_, N_AHEAD)
        fc_bond_med <- rep(NA_real_, N_AHEAD)
      }
      
      cat(sprintf("ok | sigma: gdp=%.4f cpi=%.4f bond=%.4f | lambda=%.3f | acc=%.1f%% | Q1 gdp=%+.3f cpi=%+.3f bond=%+.3f\n",
                  sigma_gdp, sigma_cpi, sigma_bond,
                  lambda_mean, acc_rate * 100,
                  fc_gdp_med[1], fc_cpi_med[1], fc_bond_med[1]))
      
      bvar_ext_results[[vs$id]][[sname]][[pkey]] <- list(
        id          = id,
        set         = vs$id,
        label       = vs$label,
        sample      = sname,
        p           = p,
        n_endo      = k,
        endo        = endo_cols_raw,
        oil         = vs$oil,
        model       = mod,
        n_obs       = n_obs,
        sigma_gdp   = sigma_gdp,
        sigma_cpi   = sigma_cpi,
        sigma_bond  = sigma_bond,
        lambda_mean = lambda_mean,
        acc_rate    = acc_rate,
        fc_gdp_med  = fc_gdp_med,
        fc_gdp_q16  = fc_gdp_q16,
        fc_gdp_q84  = fc_gdp_q84,
        fc_cpi_med  = fc_cpi_med,
        fc_bond_med = fc_bond_med
      )
      
      bvar_ext_rows[[length(bvar_ext_rows) + 1]] <- tibble(
        id           = id,
        set          = vs$id,
        label        = vs$label,
        n_endo       = k,
        oil          = vs$oil,
        sample       = sname,
        p            = p,
        n_obs        = n_obs,
        sigma_gdp    = sigma_gdp,
        sigma_cpi    = sigma_cpi,
        sigma_bond   = sigma_bond,
        lambda_mean  = round(lambda_mean, 3),
        acc_rate     = round(acc_rate, 3),
        fc_gdp_q1    = round(fc_gdp_med[1],  3),
        fc_gdp_q1_lo = round(fc_gdp_q16[1],  3),
        fc_gdp_q1_hi = round(fc_gdp_q84[1],  3),
        fc_cpi_q1    = round(fc_cpi_med[1],  3),
        fc_bond_q1   = round(fc_bond_med[1], 3)
      )
      
      n_estimated <- n_estimated + 1
    }
  }
}

cat(sprintf("\nEstimation complete: %d models estimated, %d skipped\n\n",
            n_estimated, n_skipped))

bvar_ext_summary <- bind_rows(bvar_ext_rows)

# =============================================================================
# SECTION 3: SUMMARY TABLES
# =============================================================================
cat("=== SECTION 3: SUMMARY TABLES ===\n\n")

cat("── Best model per set (lowest sigma_gdp) ────────────────────────────\n")
bvar_ext_best <- bvar_ext_summary %>%
  group_by(set) %>%
  slice_min(sigma_gdp, n = 1) %>%
  ungroup() %>%
  arrange(n_endo, sigma_gdp)

bvar_ext_best %>%
  select(set, n_endo, oil, sample, p, n_obs,
         sigma_gdp, sigma_cpi, sigma_bond,
         lambda_mean, fc_gdp_q1, fc_cpi_q1, fc_bond_q1) %>%
  print(n = Inf, width = Inf)

cat("\n── Average sigma by model size ──────────────────────────────────────\n")
bvar_ext_summary %>%
  group_by(n_endo) %>%
  summarise(
    n_models        = n(),
    mean_sigma_gdp  = round(mean(sigma_gdp,   na.rm = TRUE), 4),
    mean_sigma_cpi  = round(mean(sigma_cpi,   na.rm = TRUE), 4),
    mean_sigma_bond = round(mean(sigma_bond,  na.rm = TRUE), 4),
    min_sigma_gdp   = round(min(sigma_gdp,    na.rm = TRUE), 4),
    mean_lambda     = round(mean(lambda_mean, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(n_endo) %>%
  print(width = Inf)

cat("\n── Comparison: baseline BVAR vs best extended BVAR ──────────────────\n")
# bvar_summary only has sigma_gdp — pull best baseline row and add label
baseline_ref <- bvar_summary %>%
  filter(prior == "minnesota") %>%
  slice_min(sigma_gdp, n = 1) %>%
  select(id, sigma_gdp) %>%
  mutate(type = "baseline_bvar")

best_ext_ref <- bvar_ext_best %>%
  slice_min(sigma_gdp, n = 1) %>%
  select(id, sigma_gdp, sigma_cpi, sigma_bond) %>%
  mutate(type = "extended_bvar")

cat("Baseline BVAR (Minnesota, best):\n")
print(baseline_ref, width = Inf)
cat("\nBest extended BVAR:\n")
print(best_ext_ref, width = Inf)

# =============================================================================
# SECTION 4: OUT-OF-SAMPLE EVALUATION (expanding window)
# =============================================================================
cat("\n=== SECTION 4: EXPANDING-WINDOW RMSE ===\n\n")
cat("Horizons: h=1,2,4,8 | Min training: 40 obs | Draws: 2000\n\n")

EVAL_DRAWS <- 2000
MIN_TRAIN  <- 40
HORIZONS   <- c(1, 2, 4, 8)

bvar_ext_fc_errors <- list()

for (i in seq_len(nrow(bvar_ext_best))) {
  row <- bvar_ext_best[i, ]
  res <- bvar_ext_results[[row$set]][[row$sample]][[paste0("p", row$p)]]
  if (is.null(res)) next
  
  cat(sprintf("── Eval: %s ────────────────────────────────\n", res$id))
  
  vs       <- bvar_sets[[which(sapply(bvar_sets, function(v) v$id == res$set))]]
  endo_adj <- c("gdp_g_adj", "cpi_g_adj", "bond_dif_adj", vs$extra)
  p        <- res$p
  s        <- samples[[res$sample]]
  
  dat <- combined_data %>%
    filter(date >= s$start_date) %>%
    select(date, all_of(endo_adj)) %>%
    arrange(date)
  complete_rows <- complete.cases(dat[, endo_adj])
  dat           <- dat[complete_rows, ]
  n_total       <- nrow(dat)
  Y_full        <- as.matrix(dat[, endo_adj])
  
  exog_mat <- if (vs$oil) exog_with_oil else exog_base
  exog_s   <- exog_mat[var_adj$date %in% dat$date, , drop = FALSE]
  exog_s   <- exog_s[complete_rows, , drop = FALSE]
  
  errors <- vector("list", length(HORIZONS))
  names(errors) <- paste0("h", HORIZONS)
  
  for (h in HORIZONS) {
    t_start <- MIN_TRAIN + p
    t_end   <- n_total - h
    if (t_end <= t_start) next
    
    errs <- matrix(NA, nrow = t_end - t_start + 1, ncol = 3,
                   dimnames = list(NULL, c("gdp_g", "cpi_g", "bond_dif")))
    
    for (t in t_start:t_end) {
      Y_train    <- Y_full[1:t, ]
      exog_train <- exog_s[1:t, , drop = FALSE]
      
      mn_prior <- bv_minnesota(
        lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 1e-4, max = 5),
        alpha  = bv_alpha(mode = 2),
        psi    = bv_psi(),
        var    = 1
      )
      priors <- bv_priors(hyper = "auto", mn = mn_prior)
      
      mod_t <- tryCatch(
        bvar(
          data    = Y_train,
          lags    = p,
          n_draw  = EVAL_DRAWS + 1000,
          n_burn  = 1000,
          n_thin  = N_THIN,
          priors  = priors,
          exogen  = exog_train,
          verbose = FALSE
        ),
        error = function(e) NULL
      )
      if (is.null(mod_t)) next
      
      fc_t <- tryCatch(
        predict(mod_t, horizon = h),
        error = function(e) NULL
      )
      if (is.null(fc_t)) next
      
      fc_med <- apply(fc_t$fcast[, h, ], 2, median)
      actual <- Y_full[t + h, 1:3]
      errs[t - t_start + 1, ] <- (fc_med[1:3] - actual)^2
    }
    
    errors[[paste0("h", h)]] <- errs
    rmse <- sqrt(colMeans(errs, na.rm = TRUE))
    cat(sprintf("  h=%d: RMSE gdp=%.4f cpi=%.4f bond=%.4f\n",
                h, rmse[1], rmse[2], rmse[3]))
  }
  
  bvar_ext_fc_errors[[res$id]] <- errors
}

# =============================================================================
# SECTION 5: RMSE COMPARISON TABLE
# =============================================================================
cat("\n=== SECTION 5: RMSE COMPARISON ===\n\n")

rmse_compare <- bind_rows(lapply(names(bvar_ext_fc_errors), function(mid) {
  errs <- bvar_ext_fc_errors[[mid]]
  bind_rows(lapply(HORIZONS, function(h) {
    e <- errs[[paste0("h", h)]]
    if (is.null(e)) return(NULL)
    rmse <- sqrt(colMeans(e, na.rm = TRUE))
    tibble(
      model     = mid,
      h         = h,
      rmse_gdp  = round(rmse["gdp_g"],    4),
      rmse_cpi  = round(rmse["cpi_g"],    4),
      rmse_bond = round(rmse["bond_dif"], 4)
    )
  }))
}))

if (nrow(rmse_compare) > 0) {
  cat("Extended BVAR RMSE by model and horizon:\n")
  rmse_compare %>%
    arrange(h, rmse_gdp) %>%
    print(n = Inf, width = Inf)
}

# =============================================================================
# SECTION 6: FORECAST PLOTS — BEST MODEL PER SET
# =============================================================================
cat("\n=== SECTION 6: FORECAST PLOTS ===\n\n")

for (i in seq_len(nrow(bvar_ext_best))) {
  row <- bvar_ext_best[i, ]
  res <- bvar_ext_results[[row$set]][[row$sample]][[paste0("p", row$p)]]
  if (is.null(res)) next
  
  s           <- samples[[res$sample]]
  var_labels  <- c("GDP growth (QoQ %)", "CPI inflation (QoQ %)",
                   "Bond yield change (pp)")
  cols        <- c(COL_GDP, COL_CPI, COL_BOND)
  
  hist_dat <- var_adj %>%
    filter(date >= s$start_date) %>%
    tail(20) %>%
    select(date, gdp_g_adj, cpi_g_adj, bond_dif_adj)
  
  last_date <- max(hist_dat$date)
  fc_dates  <- seq(last_date + months(3), by = "quarter", length.out = N_AHEAD)
  
  panels <- lapply(1:3, function(vi) {
    adj_col <- c("gdp_g_adj", "cpi_g_adj", "bond_dif_adj")[vi]
    fc_med  <- list(res$fc_gdp_med, res$fc_cpi_med, res$fc_bond_med)[[vi]]
    fc_lo   <- if (vi == 1) res$fc_gdp_q16 else rep(NA, N_AHEAD)
    fc_hi   <- if (vi == 1) res$fc_gdp_q84 else rep(NA, N_AHEAD)
    col     <- cols[vi]
    
    hist_df <- hist_dat %>% select(date, value = !!sym(adj_col))
    fc_df   <- tibble(date  = fc_dates,
                      fcst  = fc_med,
                      lower = fc_lo,
                      upper = fc_hi)
    
    p <- ggplot() +
      geom_hline(yintercept = 0, linetype = "dashed",
                 color = COL_GREY, linewidth = 0.3)
    
    if (vi == 1) {
      p <- p + geom_ribbon(data = fc_df,
                           aes(x = date, ymin = lower, ymax = upper),
                           fill = col, alpha = 0.15)
    }
    
    p +
      geom_line(data = hist_df, aes(x = date, y = value),
                colour = col, linewidth = 0.8) +
      geom_point(data = hist_df, aes(x = date, y = value),
                 colour = col, size = 1.5) +
      geom_line(data = fc_df, aes(x = date, y = fcst),
                colour = col, linewidth = 0.8, linetype = "dashed") +
      annotate("text",
               x     = fc_df$date[1],
               y     = fc_df$fcst[1] + 0.08,
               label = sprintf("%+.2f", fc_df$fcst[1]),
               size  = 3.2, colour = col, fontface = "bold") +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
      labs(title = var_labels[vi], x = NULL, y = NULL)
  })
  
  fig <- wrap_plots(panels, ncol = 3) +
    plot_annotation(
      title    = sprintf("BVAR forecast — %s", res$id),
      subtitle = sprintf("%s | p=%d | %d-var%s | Q1 2026 GDP: %+.2f%%",
                         s$label, res$p, res$n_endo,
                         if (res$oil) " (Brent exog.)" else "",
                         res$fc_gdp_med[1]),
      theme = theme(
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "gray40")
      )
    )
  print(fig)
  
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("bvar_04_fc_%s.png", res$id)),
  #        fig, width = 12, height = 4, dpi = 150)
}

# =============================================================================
# SECTION 7: SAVE OUTPUTS
# =============================================================================
# write.csv(bvar_ext_summary,
#           file.path(here("output","tables"), "bvar_04_summary.csv"),
#           row.names = FALSE)

cat("\nScript complete.\n")
cat(sprintf("Objects: bvar_ext_results (%d sets) | bvar_ext_summary (%d rows)\n",
            length(bvar_ext_results), nrow(bvar_ext_summary)))
cat(sprintf("         bvar_ext_best (%d rows) | bvar_ext_fc_errors (%d models)\n",
            nrow(bvar_ext_best), length(bvar_ext_fc_errors)))
cat("Next: 05_bvar_extended_evaluation.R\n")