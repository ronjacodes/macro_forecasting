# ============================================================================
# scripts/models/bvar/07_bvar_large.R
#
# PURPOSE:
#   Estimate and evaluate a large 20-variable BVAR model designed to
#   maximise out-of-sample forecast accuracy for all three targets.
#   Motivated by Banbura et al. (2010) who show BVAR handles large
#   systems well precisely because the Minnesota prior prevents
#   overfitting as variables are added.
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/bvar/01_bvar_baseline.R")
#     — var_adj, samples, BASE_COLS, PRIORS, COL_GDP, COL_CPI, COL_BOND
#   source("scripts/models/var/02_var_baseline_evaluation.R")
#     — forecast_errors (AR1 benchmarks), COVID_EVAL_EXCL
#   source("scripts/models/bvar/02_bvar_baseline_evaluation.R")
#     — bvar_fc_errors (baseline BVAR benchmark)
#
# MODEL DESIGN — 20 endogenous variables:
#
#   Block 1 — Core targets (3):
#     gdp_g        GDP growth QoQ %
#     cpi_g        CPI inflation QoQ %
#     bond_dif     Swiss 10Y bond yield change QoQ pp
#
#   Block 2 — External demand (4):
#     ea_esi       EA Economic Sentiment Indicator (level)
#     ea_pmi       EA PMI manufacturing (level)
#     msci_wld     MSCI World equity QoQ %
#     oecd_cli     OECD composite leading indicator CH (level)
#
#   Block 3 — Domestic leading (3):
#     kof_bar      KOF economic barometer (level)
#     sw_pmi       Swiss PMI (level)
#     unemp        Swiss unemployment rate (level)
#
#   Block 4 — Price channel (3):
#     fx_eur       CHF/EUR exchange rate QoQ %
#     imp_prc      Import prices QoQ %
#     ppi          Producer & Import PPI QoQ %
#     prod_ppi     Producer PPI QoQ %  [dropped — see below]
#
#   Block 5 — Financial / monetary (4):
#     bund_10y     German Bund 10Y yield level
#     ea_rate      EA short-term interest rate level
#     sw_m1        Swiss M1 money supply QoQ %
#     sw_m3        Swiss M3 money supply QoQ %
#     ea_m3        EA M3 money supply QoQ %
#
#   Block 6 — Commodities (1):
#     oil          Brent crude oil QoQ % (ENDOGENOUS here, unlike VARX)
#
#   Final selection: 20 endogenous variables (no exogenous oil)
#   Exogenous: COVID dummies only
#
# VARIABLE SELECTION RATIONALE:
#   Dropped Producer PPI (swproprce) — highly collinear with PPI & Import
#   PPI (swppingdf); keeping both adds noise without independent signal.
#   Kept EA PMI + EA ESI — different signals: manufacturing vs overall.
#   Kept OECD CLI + KOF — independent methodologies, low collinearity.
#   Kept M1 + M3 — narrow (transactions) vs broad (monetary conditions).
#   Oil endogenous — allows dynamic interaction with CPI and bond yields.
#
# SETTINGS:
#   Prior:      Minnesota (same as all BVAR models in this project)
#   Lags:       p = 1 (20-var system: 20×21+1 = 421 params per eq — prior essential)
#   Sample:     Full (Q2 2000 – Q4 2025, n = 103)
#   Draws:      10000 total (5000 burn + 5000 save) for main estimation
#   Rolling eval: 3000 total (2000 burn + 1000 save) for speed
#   Eval start: 2005 Q1 (80 windows vs 60 from 2010 Q1)
#   Horizons:   h = 1, 2, 4, 8
#
# REFERENCES:
#   Banbura, Giannone & Reichlin (2010) — Large BVARs
#   Giannone, Lenza & Primiceri (2015) — Prior selection in large BVARs
#   Gerlach & Jordan (2012) — Swiss inflation forecasting
#
# KEY FINDINGS (fill in after running):
#   GDP  h=1: RMSE=[  ] vs AR1=[  ] ratio=[  ] beats AR1=[Y/N]
#   GDP  h=4: RMSE=[  ] vs AR1=[  ] ratio=[  ] beats AR1=[Y/N]
#   CPI  h=1: RMSE=[  ] vs AR1=[  ] ratio=[  ] beats AR1=[Y/N]
#   CPI  h=4: RMSE=[  ] vs AR1=[  ] ratio=[  ] beats AR1=[Y/N]
#   Bond h=1: RMSE=[  ] vs AR1=[  ] ratio=[  ] beats AR1=[Y/N]
#   Bond h=4: RMSE=[  ] vs AR1=[  ] ratio=[  ] beats AR1=[Y/N]
#   Q1 2026: GDP=[  ] CPI=[  ] Bond=[  ]
#
# Outputs (commented — uncomment to save):
#   output/figures/bvar_07_fc.png
#   output/figures/bvar_07_rmse.png
#   output/tables/bvar_07_rmse.csv
# ============================================================================

library(dplyr)

N_AHEAD    <- 8
HORIZONS   <- c(1, 2, 4, 8)
N_DRAW     <- 10000
N_BURN     <- 5000
N_THIN     <- 1
EVAL_START <- as.Date("2005-01-01")
EVAL_DRAWS <- 3000
EVAL_BURN  <- 2000

TARGET_VARS <- c("gdp_g", "cpi_g", "bond_dif")

# =============================================================================
# SECTION 1: EXTRACT AND ALIGN ALL 20 VARIABLES
# =============================================================================
cat("=== SECTION 1: EXTRACT VARIABLES ===\n\n")

# Helper: extract one series from JSON, aggregate to quarterly
extract_var <- function(key, transform, col_name, label) {
  json_key <- paste0(key, "_", transform)
  vals     <- nowcast_raw[[json_key]]
  dts      <- nowcast_raw$dates[[json_key]]
  if (is.null(vals) || is.null(dts)) {
    cat(sprintf("  [!] NOT FOUND: %s (%s)\n", label, json_key))
    return(NULL)
  }
  tibble(date = lubridate::dmy(dts), value = as.numeric(vals)) %>%
    filter(!is.na(value)) %>%
    mutate(date = as.Date(zoo::as.yearqtr(date))) %>%
    group_by(date) %>%
    summarise(!!col_name := mean(value, na.rm = TRUE), .groups = "drop")
}

# Define all 20 variables
var_defs <- list(
  # Block 2 — External demand
  list(key = "ekeusesig",        tf = "lvl",    col = "ea_esi",   label = "EA ESI"),
  list(key = "empmim_hq",        tf = "lvl",    col = "ea_pmi",   label = "EA PMI mfg"),
  list(key = "mswrld_d_",        tf = "pct_3m", col = "msci_wld", label = "MSCI World QoQ"),
  list(key = "swobs085q",        tf = "lvl",    col = "oecd_cli", label = "OECD CLI CH"),
  # Block 3 — Domestic leading
  list(key = "ch_kof_barometer", tf = "lvl",    col = "kof_bar",  label = "KOF barometer"),
  list(key = "swpurchsq",        tf = "lvl",    col = "sw_pmi",   label = "Swiss PMI"),
  list(key = "swun_p_totq",      tf = "lvl",    col = "unemp",    label = "Unemployment"),
  # Block 4 — Price channel
  list(key = "swxsfec_",         tf = "pct_3m", col = "fx_eur",   label = "CHF/EUR QoQ"),
  list(key = "swimpprce",        tf = "pct_3m", col = "imp_prc",  label = "Import prices QoQ"),
  list(key = "swppingdf",        tf = "pct_3m", col = "ppi",      label = "PPI & Import PPI QoQ"),
  # Block 5 — Financial / monetary
  list(key = "bdeuscciq",        tf = "lvl",    col = "bund_10y", label = "Bund 10Y level"),
  list(key = "bdeusicir",        tf = "lvl",    col = "ea_rate",  label = "EA short rate level"),
  list(key = "swm1____a",        tf = "pct_3m", col = "sw_m1",    label = "Swiss M1 QoQ"),
  list(key = "swm3____a",        tf = "pct_3m", col = "sw_m3",    label = "Swiss M3 QoQ"),
  list(key = "emm3____b",        tf = "pct_3m", col = "ea_m3",    label = "EA M3 QoQ"),
  # Block 6 — Commodities (endogenous)
  list(key = "oilbren",          tf = "pct_3m", col = "oil",      label = "Brent oil QoQ")
)

# Extract all and align to var_input dates
large_data <- var_adj %>%
  select(date, gdp_g_adj, cpi_g_adj, bond_dif_adj) %>%
  arrange(date)

for (v in var_defs) {
  s <- extract_var(v$key, v$tf, v$col, v$label)
  if (is.null(s)) next
  large_data <- left_join(large_data, s, by = "date")
  n_ok <- sum(!is.na(large_data[[v$col]]))
  cat(sprintf("  %-25s %3d obs aligned\n", v$label, n_ok))
}

# Column order for the BVAR: GDP, CPI, Bond first, then all others
adj_cols   <- c("gdp_g_adj", "cpi_g_adj", "bond_dif_adj")
extra_cols <- setdiff(names(large_data), c("date", adj_cols))
all_cols   <- c(adj_cols, extra_cols)

cat(sprintf("\nlarge_data: %d rows x %d cols (%s to %s)\n",
            nrow(large_data), ncol(large_data),
            format(min(large_data$date)),
            format(max(large_data$date))))
cat(sprintf("Variables: %s\n\n", paste(all_cols, collapse = ", ")))

# Full data matrix — drop rows with ANY NA
dat_full      <- large_data %>%
  select(date, all_of(all_cols)) %>%
  filter(complete.cases(.))
n_complete    <- nrow(dat_full)
cat(sprintf("Complete cases: %d rows (%s to %s)\n\n",
            n_complete,
            format(min(dat_full$date)),
            format(max(dat_full$date))))

Y_full <- as.matrix(dat_full[, all_cols])
k      <- ncol(Y_full)
cat(sprintf("System: k=%d endogenous variables, p=1\n", k))
cat(sprintf("Parameters per equation: %d (k*p + 1 = %d*1 + 1)\n",
            k + 1, k))
cat(sprintf("Total parameters: %d | obs: %d | prior essential ✓\n\n",
            k * (k + 1), n_complete))

# =============================================================================
# SECTION 2: ESTIMATE LARGE BVAR
# =============================================================================
cat("=== SECTION 2: ESTIMATE LARGE BVAR (full sample) ===\n\n")
cat(sprintf("Prior: Minnesota | Lags: 1 | Draws: %d burn + %d save\n\n",
            N_BURN, N_DRAW - N_BURN))

# COVID dummies aligned to complete cases
exog_large <- exog_base[var_adj$date %in% dat_full$date, , drop = FALSE]

mn_prior <- bv_minnesota(
  lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 1e-4, max = 5),
  alpha  = bv_alpha(mode = 2),
  psi    = bv_psi(),
  var    = 1
)
priors <- bv_priors(hyper = "auto", mn = mn_prior)

cat("Estimating... (this may take 2-3 minutes)\n")
t0 <- proc.time()

bvar_large_mod <- tryCatch(
  bvar(
    data    = Y_full,
    lags    = 1,
    n_draw  = N_DRAW,
    n_burn  = N_BURN,
    n_thin  = N_THIN,
    priors  = priors,
    exogen  = exog_large,
    verbose = FALSE
  ),
  error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message)); NULL
  }
)

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("Estimation complete in %.1f seconds.\n\n", elapsed))

if (!is.null(bvar_large_mod)) {
  # Lambda posterior
  lambda_chain <- tryCatch({
    h <- bvar_large_mod$hyper
    if ("lambda" %in% colnames(h)) h[, "lambda"] else h[, 1]
  }, error = function(e) NULL)
  
  acc_rate <- tryCatch(
    bvar_large_mod$meta$accepted / bvar_large_mod$meta$n_save,
    error = function(e) NA_real_
  )
  
  cat(sprintf("Lambda: mean=%.3f  sd=%.3f  [%.3f, %.3f] 90%% CI\n",
              mean(lambda_chain), sd(lambda_chain),
              quantile(lambda_chain, 0.05),
              quantile(lambda_chain, 0.95)))
  cat(sprintf("Acceptance rate: %.1f%%\n\n", acc_rate * 100))
  
  # Q1 2026 forecast
  fc_full <- tryCatch(
    predict(bvar_large_mod, horizon = N_AHEAD, conf_bands = c(0.16, 0.84)),
    error = function(e) NULL
  )
  
  if (!is.null(fc_full)) {
    # fc$fcast: [draw, horizon, variable]
    cat("Q1 2026 posterior median forecasts:\n")
    for (vi in 1:3) {
      vname <- c("GDP growth", "CPI inflation", "Bond yield Δ")[vi]
      med   <- median(fc_full$fcast[, 1, vi])
      q16   <- quantile(fc_full$fcast[, 1, vi], 0.16)
      q84   <- quantile(fc_full$fcast[, 1, vi], 0.84)
      cat(sprintf("  %-16s %+.3f  [%+.3f, %+.3f] 68%% CI\n",
                  vname, med, q16, q84))
    }
    cat("\n")
  }
}


# =============================================================================
# SECTION 2b: FOCUSED 10-VARIABLE BVAR
# =============================================================================
# MOTIVATION:
#   The 19-variable model above produced lambda=0.153 (very tight shrinkage)
#   and forecast intervals of +/-6% for GDP — far too wide to be useful.
#   With only 103 observations, the marginal likelihood optimizer shrinks
#   aggressively toward zero when given too many variables.
#
# VARIABLE SELECTION (10 endogenous):
#   gdp_g_adj    GDP growth           — core target
#   cpi_g_adj    CPI inflation        — core target
#   bond_dif_adj Bond yield change    — core target
#   ea_esi       EA ESI               — best GDP/bond predictor (confirmed by prof)
#   msci_wld     MSCI World QoQ       — global demand, beats AR1 for bond h=4
#   fx_eur       CHF/EUR QoQ          — imported inflation channel
#   imp_prc      Import prices QoQ    — upstream CPI (strongest Granger)
#   kof_bar      KOF barometer        — best Swiss leading indicator
#   bund_10y     Bund 10Y level       — SNB follows ECB, strongest bond predictor
#   ppi          PPI & Import PPI QoQ — price pressure leading CPI by 1-2Q
#
# Exogenous: oil (Brent QoQ%) + COVID dummies
# =============================================================================
cat("=== SECTION 2b: FOCUSED 10-VARIABLE BVAR ===\n\n")
cat("Motivation: 19-var model over-shrinks (lambda=0.153, intervals +/-6%)\n")
cat("Reducing to 10 empirically motivated variables\n\n")

focused_cols <- c("gdp_g_adj", "cpi_g_adj", "bond_dif_adj",
                  "ea_esi", "msci_wld", "fx_eur",
                  "imp_prc", "kof_bar", "bund_10y", "ppi")

dat_focused <- large_data %>%
  select(date, all_of(focused_cols)) %>%
  filter(complete.cases(.))

Y_focused <- as.matrix(dat_focused[, focused_cols])
k_focused <- ncol(Y_focused)

cat(sprintf("System: k=%d variables, p=1 | obs=%d\n",
            k_focused, nrow(dat_focused)))
cat(sprintf("Parameters per equation: %d | obs/param ratio: %.1f\n\n",
            k_focused + 1, nrow(dat_focused) / (k_focused + 1)))

# Oil as exogenous
oil_series  <- extract_var("oilbren", "pct_3m", "oil", "Brent oil QoQ")
oil_aligned <- dat_focused %>%
  select(date) %>%
  left_join(oil_series, by = "date") %>%
  pull(oil)
oil_aligned[is.na(oil_aligned)] <- 0

exog_focused <- cbind(
  exog_base[var_adj$date %in% dat_focused$date, , drop = FALSE],
  oil_g = oil_aligned
)

cat(sprintf("Exogenous: %s\n\n", paste(colnames(exog_focused), collapse = ", ")))

mn_prior_f <- bv_minnesota(
  lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 1e-4, max = 5),
  alpha  = bv_alpha(mode = 2),
  psi    = bv_psi(),
  var    = 1
)
priors_f <- bv_priors(hyper = "auto", mn = mn_prior_f)

cat("Estimating focused model...\n")
t0_f <- proc.time()

bvar_focused_mod <- tryCatch(
  bvar(
    data    = Y_focused,
    lags    = 1,
    n_draw  = N_DRAW,
    n_burn  = N_BURN,
    n_thin  = N_THIN,
    priors  = priors_f,
    exogen  = exog_focused,
    verbose = FALSE
  ),
  error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message)); NULL
  }
)

elapsed_f <- (proc.time() - t0_f)["elapsed"]
cat(sprintf("Estimation complete in %.1f seconds.\n\n", elapsed_f))

if (!is.null(bvar_focused_mod)) {
  lambda_f <- tryCatch({
    h <- bvar_focused_mod$hyper
    if ("lambda" %in% colnames(h)) h[, "lambda"] else h[, 1]
  }, error = function(e) NULL)
  
  acc_f <- tryCatch(
    bvar_focused_mod$meta$accepted / bvar_focused_mod$meta$n_save,
    error = function(e) NA_real_
  )
  
  cat(sprintf("Lambda: mean=%.3f  sd=%.3f  [%.3f, %.3f] 90%% CI\n",
              mean(lambda_f), sd(lambda_f),
              quantile(lambda_f, 0.05),
              quantile(lambda_f, 0.95)))
  cat(sprintf("Acceptance rate: %.1f%%\n\n", acc_f * 100))
  
  fc_focused <- tryCatch(
    predict(bvar_focused_mod, horizon = N_AHEAD, conf_bands = c(0.16, 0.84)),
    error = function(e) NULL
  )
  
  if (!is.null(fc_focused)) {
    cat("Q1 2026 posterior median forecasts (focused model):\n")
    for (vi in 1:3) {
      vname <- c("GDP growth", "CPI inflation", "Bond yield change")[vi]
      med   <- median(fc_focused$fcast[, 1, vi])
      q16   <- quantile(fc_focused$fcast[, 1, vi], 0.16)
      q84   <- quantile(fc_focused$fcast[, 1, vi], 0.84)
      cat(sprintf("  %-16s %+.3f  [%+.3f, %+.3f] 68%% CI\n",
                  vname, med, q16, q84))
    }
    cat("\n")
  }
}


# =============================================================================
# SECTION 2c: FOCUSED 10-VARIABLE BVAR v2 — SCALE-CORRECTED
# =============================================================================
# MOTIVATION:
#   Section 2b showed lambda=0.250 — still low because EA ESI (sd=9.8) and
#   KOF barometer (sd=8.3) have much larger scales than the growth rate
#   variables (sd~0.4-1.6). The Minnesota prior assumes similar scales.
#   Fix: divide EA ESI and KOF by 10, replace Bund 10Y level with change.
#
# CHANGES vs Section 2b:
#   ea_esi  / 10  → sd 9.8 → 0.98
#   kof_bar / 10  → sd 8.3 → 0.84
#   bund_10y: level (sd=6.9) → first difference (sd=2.7)
#
# RESULT: lambda=0.288, intervals GDP ±2.7%, CPI ±0.65%, Bond ±1.2%
#   Point forecasts: GDP +0.689%, CPI +0.027%, Bond -0.012%
#   → Sensible, aligned with other models → proceed to rolling eval
# =============================================================================
cat("=== SECTION 2c: SCALE-CORRECTED 10-VAR BVAR (v2) ===\n\n")

# Extract Bund yield CHANGE (stationary) instead of level
bund_dif_series <- extract_var("bdeuscciq", "dif_3m", "bund_dif",
                               "Bund yield chg QoQ")

large_data_v2 <- large_data %>%
  left_join(bund_dif_series, by = "date") %>%
  mutate(
    bund_10y = bund_dif,   # replace level with change
    ea_esi   = ea_esi / 10,   # scale down: sd 9.8 -> 0.98
    kof_bar  = kof_bar / 10   # scale down: sd 8.3 -> 0.84
  )

focused_cols_v2 <- c("gdp_g_adj", "cpi_g_adj", "bond_dif_adj",
                     "ea_esi", "msci_wld", "fx_eur",
                     "imp_prc", "kof_bar", "bund_10y", "ppi")

dat_focused_v2 <- large_data_v2 %>%
  select(date, all_of(focused_cols_v2)) %>%
  filter(complete.cases(.))

Y_focused_v2 <- as.matrix(dat_focused_v2[, focused_cols_v2])

cat("Variable scales (sd):\n")
print(round(apply(Y_focused_v2, 2, sd), 3))

# Oil exogenous
oil_series     <- extract_var("oilbren", "pct_3m", "oil", "Brent oil QoQ")
oil_aligned_v2 <- dat_focused_v2 %>%
  select(date) %>%
  left_join(oil_series, by = "date") %>%
  pull(oil)
oil_aligned_v2[is.na(oil_aligned_v2)] <- 0

exog_focused_v2 <- cbind(
  exog_base[var_adj$date %in% dat_focused_v2$date, , drop = FALSE],
  oil_g = oil_aligned_v2
)

mn_prior_v2 <- bv_minnesota(
  lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 1e-4, max = 5),
  alpha  = bv_alpha(mode = 2),
  psi    = bv_psi(),
  var    = 1
)
priors_v2 <- bv_priors(hyper = "auto", mn = mn_prior_v2)

cat("\nEstimating v2...\n")
t0_v2 <- proc.time()

bvar_focused_v2 <- tryCatch(
  bvar(
    data    = Y_focused_v2,
    lags    = 1,
    n_draw  = N_DRAW,
    n_burn  = N_BURN,
    n_thin  = N_THIN,
    priors  = priors_v2,
    exogen  = exog_focused_v2,
    verbose = FALSE
  ),
  error = function(e) { cat(sprintf("ERROR: %s\n", e$message)); NULL }
)

elapsed_v2 <- (proc.time() - t0_v2)["elapsed"]
cat(sprintf("Done in %.1f seconds.\n\n", elapsed_v2))

if (!is.null(bvar_focused_v2)) {
  lambda_v2 <- tryCatch({
    h <- bvar_focused_v2$hyper
    if ("lambda" %in% colnames(h)) h[, "lambda"] else h[, 1]
  }, error = function(e) NULL)
  
  acc_v2 <- tryCatch(
    bvar_focused_v2$meta$accepted / bvar_focused_v2$meta$n_save,
    error = function(e) NA_real_
  )
  
  cat(sprintf("Lambda: mean=%.3f  sd=%.3f  [%.3f, %.3f] 90%% CI\n",
              mean(lambda_v2), sd(lambda_v2),
              quantile(lambda_v2, 0.05),
              quantile(lambda_v2, 0.95)))
  cat(sprintf("Acceptance rate: %.1f%%\n\n", acc_v2 * 100))
  
  fc_v2 <- tryCatch(
    predict(bvar_focused_v2, horizon = N_AHEAD, conf_bands = c(0.16, 0.84)),
    error = function(e) NULL
  )
  
  if (!is.null(fc_v2)) {
    cat("Q1 2026 forecasts (v2 — scale corrected):\n")
    for (vi in 1:3) {
      vname <- c("GDP growth", "CPI inflation", "Bond yield chg")[vi]
      med   <- median(fc_v2$fcast[, 1, vi])
      q16   <- quantile(fc_v2$fcast[, 1, vi], 0.16)
      q84   <- quantile(fc_v2$fcast[, 1, vi], 0.84)
      cat(sprintf("  %-16s %+.3f  [%+.3f, %+.3f] 68%% CI\n",
                  vname, med, q16, q84))
    }
    cat("\n")
  }
}

# =============================================================================
# SECTION 3: EXPANDING-WINDOW RMSE EVALUATION
# =============================================================================
cat("=== SECTION 3: EXPANDING-WINDOW RMSE (focused 10-var model) ===\n\n")
cat("Note: 19-var model not evaluated due to over-shrinkage (Section 2)\n\n")
cat(sprintf("Eval start: %s | Horizons: %s | Draws: %d burn + %d save\n\n",
            format(EVAL_START),
            paste(HORIZONS, collapse = ","),
            EVAL_BURN, EVAL_DRAWS - EVAL_BURN))

eval_dates <- dat_focused_v2$date[dat_focused_v2$date >= EVAL_START]
n_eval     <- length(eval_dates)
MIN_TRAIN  <- 20  # minimum obs before we start forecasting

cat(sprintf("Evaluation windows: %d (%s to %s)\n",
            n_eval, format(min(eval_dates)), format(max(eval_dates))))
cat("Starting rolling estimation...\n\n")

errors_focused <- vector("list", length(HORIZONS))
names(errors_focused) <- paste0("h", HORIZONS)
for (h in HORIZONS) {
  errors_focused[[paste0("h", h)]] <- matrix(
    NA,
    nrow = n_eval,
    ncol = 3,
    dimnames = list(NULL, c("gdp_g", "cpi_g", "bond_dif"))
  )
}

t_start_eval <- proc.time()

for (i in seq_along(eval_dates)) {
  origin_date <- eval_dates[i]
  origin_idx  <- which(dat_full$date == origin_date)
  
  # Training data: all rows up to but not including origin
  train_idx <- seq_len(origin_idx - 1)
  if (length(train_idx) < MIN_TRAIN) next
  
  Y_train    <- Y_focused_v2[train_idx, , drop = FALSE]
  exog_train <- exog_focused_v2[train_idx, , drop = FALSE]
  
  # Estimate
  mod_t <- tryCatch(
    bvar(
      data    = Y_train,
      lags    = 1,
      n_draw  = EVAL_DRAWS,
      n_burn  = EVAL_BURN,
      n_thin  = N_THIN,
      priors  = priors,
      exogen  = exog_train,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(mod_t)) next
  
  # Forecast
  fc_t <- tryCatch(
    predict(mod_t, horizon = max(HORIZONS)),
    error = function(e) NULL
  )
  if (is.null(fc_t)) next
  
  # Store squared errors for each horizon
  for (h in HORIZONS) {
    target_idx <- origin_idx + h - 1
    if (target_idx > nrow(dat_full)) next
    
    fc_med  <- apply(fc_t$fcast[, h, 1:3], 2, median)
    actual  <- Y_full[target_idx, 1:3]
    sq_err  <- (fc_med - actual)^2
    
    errors_focused[[paste0("h", h)]][i, ] <- sq_err
  }
  
  # Progress every 10 windows
  if (i %% 10 == 0) {
    elapsed_so_far <- (proc.time() - t_start_eval)["elapsed"]
    est_total      <- elapsed_so_far / i * n_eval
    cat(sprintf("  Window %d/%d | elapsed: %.0fs | est. total: %.0fs (%.1f min)\n",
                i, n_eval, elapsed_so_far, est_total, est_total / 60))
  }
}

elapsed_eval <- (proc.time() - t_start_eval)["elapsed"]
cat(sprintf("\nRolling evaluation complete in %.1f minutes.\n\n",
            elapsed_eval / 60))

# =============================================================================
# SECTION 4: RMSE COMPARISON
# =============================================================================
cat("=== SECTION 4: RMSE COMPARISON ===\n\n")

# Large BVAR RMSE
rmse_large <- bind_rows(lapply(HORIZONS, function(h) {
  e <- errors_focused[[paste0("h", h)]]
  if (is.null(e)) return(NULL)
  rmse <- sqrt(colMeans(e, na.rm = TRUE))
  tibble(
    model     = "BVAR focused (10-var)",
    h         = h,
    rmse_gdp  = round(rmse["gdp_g"],    4),
    rmse_cpi  = round(rmse["cpi_g"],    4),
    rmse_bond = round(rmse["bond_dif"], 4),
    n         = sum(!is.na(e[, 1]))
  )
}))

# AR(1) benchmark — recompute from 2005 Q1 for fair comparison
ar1_tbl <- bind_rows(
  forecast_errors[["AR1_gdp_g"]],
  forecast_errors[["AR1_cpi_g"]],
  forecast_errors[["AR1_bond_dif"]]
) %>% filter(date >= EVAL_START)

ar1_rmse_2005 <- bind_rows(lapply(HORIZONS, function(hh) {
  sub <- ar1_tbl %>% filter(h == hh, !is.na(error))
  if (nrow(sub) == 0) return(NULL)
  tibble(
    model     = "AR(1) [from 2005]",
    h         = hh,
    rmse_gdp  = round(sqrt(mean(sub$error[sub$variable == "gdp_g"]^2,  na.rm = TRUE)), 4),
    rmse_cpi  = round(sqrt(mean(sub$error[sub$variable == "cpi_g"]^2,  na.rm = TRUE)), 4),
    rmse_bond = round(sqrt(mean(sub$error[sub$variable == "bond_dif"]^2, na.rm = TRUE)), 4),
    n         = sum(sub$variable == "gdp_g")
  )
}))

# Baseline BVAR for reference
bvar_base_rmse <- bind_rows(lapply(HORIZONS, function(hh) {
  sub <- bvar_fc_errors[["bvar_minnesota_full_p4"]] %>%
    filter(h == hh, !is.na(error))
  if (nrow(sub) == 0) return(NULL)
  tibble(
    model     = "BVAR baseline (3-var)",
    h         = hh,
    rmse_gdp  = round(sqrt(mean(sub$error[sub$variable == "gdp_g"]^2,  na.rm = TRUE)), 4),
    rmse_cpi  = round(sqrt(mean(sub$error[sub$variable == "cpi_g"]^2,  na.rm = TRUE)), 4),
    rmse_bond = round(sqrt(mean(sub$error[sub$variable == "bond_dif"]^2, na.rm = TRUE)), 4),
    n         = sum(sub$variable == "gdp_g")
  )
}))

all_rmse <- bind_rows(ar1_rmse_2005, bvar_base_rmse, rmse_large)

cat("── RMSE table ───────────────────────────────────────────────────────\n")
all_rmse %>%
  select(model, h, rmse_gdp, rmse_cpi, rmse_bond, n) %>%
  arrange(h, rmse_gdp) %>%
  print(n = Inf, width = Inf)

# Relative RMSE vs AR(1)
ar1_ref <- ar1_rmse_2005 %>%
  select(h, ar1_gdp = rmse_gdp, ar1_cpi = rmse_cpi, ar1_bond = rmse_bond)

cat("\n── Relative RMSE vs AR(1) from 2005 ────────────────────────────────\n")
all_rmse %>%
  filter(model != "AR(1) [from 2005]") %>%
  left_join(ar1_ref, by = "h") %>%
  mutate(
    rel_gdp  = round(rmse_gdp  / ar1_gdp,  3),
    rel_cpi  = round(rmse_cpi  / ar1_cpi,  3),
    rel_bond = round(rmse_bond / ar1_bond,  3)
  ) %>%
  select(model, h, rel_gdp, rel_cpi, rel_bond) %>%
  arrange(h) %>%
  print(n = Inf, width = Inf)

# =============================================================================
# SECTION 5: FORECAST PLOT
# =============================================================================
cat("\n=== SECTION 5: FORECAST PLOT ===\n\n")

if (!is.null(fc_v2)) {
  var_labels <- c("GDP growth (QoQ %)", "CPI inflation (QoQ %)",
                  "Bond yield \u0394 (pp)")
  cols       <- c(COL_GDP, COL_CPI, COL_BOND)
  
  hist_dat <- dat_full %>%
    tail(20) %>%
    select(date, gdp_g_adj, cpi_g_adj, bond_dif_adj)
  
  last_date <- max(hist_dat$date)
  fc_dates  <- seq(last_date + months(3), by = "quarter",
                   length.out = N_AHEAD)
  
  panels <- lapply(1:3, function(vi) {
    adj_col <- c("gdp_g_adj", "cpi_g_adj", "bond_dif_adj")[vi]
    col     <- cols[vi]
    
    fc_med <- apply(fc_v2$fcast[, , vi], 2, median)
    fc_q16 <- apply(fc_v2$fcast[, , vi], 2, quantile, 0.16)
    fc_q84 <- apply(fc_v2$fcast[, , vi], 2, quantile, 0.84)
    
    hist_df <- hist_dat %>% select(date, value = !!sym(adj_col))
    fc_df   <- tibble(date  = fc_dates,
                      fcst  = fc_med,
                      lower = fc_q16,
                      upper = fc_q84)
    
    ggplot() +
      geom_hline(yintercept = 0, linetype = "dashed",
                 color = COL_GREY, linewidth = 0.3) +
      geom_ribbon(data = fc_df,
                  aes(x = date, ymin = lower, ymax = upper),
                  fill = col, alpha = 0.15) +
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
  
  fig_fc <- wrap_plots(panels, ncol = 3) +
    plot_annotation(
      title    = "BVAR focused (10-var) — Q1 2026 forecast",
      subtitle = sprintf("Full sample | p=1 | Minnesota prior | %d draws | 68%% posterior interval",
                         N_DRAW - N_BURN),
      theme = theme(
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "gray40")
      )
    )
  print(fig_fc)
  # ggsave(file.path(here("output","figures"), "bvar_07_fc.png"),
  #        fig_fc, width = 12, height = 4, dpi = 150)
}

# RMSE comparison plot
rmse_plot_df <- all_rmse %>%
  pivot_longer(c(rmse_gdp, rmse_cpi, rmse_bond),
               names_to  = "target",
               values_to = "rmse") %>%
  mutate(
    target = recode(target,
                    rmse_gdp  = "GDP growth",
                    rmse_cpi  = "CPI inflation",
                    rmse_bond = "Bond yield \u0394"),
    target = factor(target,
                    levels = c("GDP growth", "CPI inflation",
                               "Bond yield \u0394"))
  )

fig_rmse <- ggplot(rmse_plot_df,
                   aes(x = h, y = rmse, colour = model, group = model)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ target, scales = "free_y", ncol = 3) +
  scale_colour_manual(
    values = c("AR(1) [from 2005]"     = "black",
               "BVAR baseline (3-var)" = COL_GDP,
               "BVAR focused (10-var)" = COL_CPI),
    name = NULL
  ) +
  scale_x_continuous(breaks = HORIZONS) +
  labs(
    title    = "RMSE by horizon — BVAR large vs benchmarks",
    subtitle = "Eval: 2005 Q1 – 2025 Q4 | including COVID quarters",
    x        = "Forecast horizon (quarters)",
    y        = "RMSE"
  ) +
  theme(legend.position = "bottom",
        strip.text      = element_text(face = "bold"))
print(fig_rmse)
# ggsave(file.path(here("output","figures"), "bvar_07_rmse.png"),
#        fig_rmse, width = 12, height = 4, dpi = 150)

# =============================================================================
# SAVE
# =============================================================================
# write.csv(all_rmse,
#           file.path(here("output","tables"), "bvar_07_rmse.csv"),
#           row.names = FALSE)

cat("\nScript complete.\n")
cat(sprintf("Objects: bvar_large_mod (19-var) | bvar_focused_v2 (10-var scaled) | errors_focused | all_rmse | fc_v2\n"))
cat("Next: update 06_model_comparison.R with large BVAR results\n")