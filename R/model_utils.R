# ============================================================
# R/model_utils.R
# Shared modelling and evaluation utilities
# ============================================================
#
# USAGE:
#   source(here("R", "model_utils.R"))
#
# ============================================================


# ------------------------------------------------------------
# bridge_nowcast()
# Compute a bridge model nowcast from the GDP equation of an
# estimated VAR, conditioning on current-quarter predictor
# averages instead of using the VAR's internal forecast.
#
# WHY NOT USE predict() FROM THE vars PACKAGE:
#   predict() forecasts all variables jointly using lagged
#   values only — it cannot condition on new within-quarter
#   monthly data. The bridge model instead extracts the GDP
#   equation coefficients and plugs in the current-quarter
#   predictor averages directly. This is what makes it a
#   "bridge": the monthly data bridges the gap to the
#   quarterly GDP observation.
#
# Arguments:
#   var_mod    — VAR model object from vars::VAR()
#   train      — data.frame, the training data used to fit
#                var_mod (needs gdp + all predictor columns)
#   pred_avgs  — 1-row data.frame of current-quarter predictor
#                averages from get_partial_avg()
#   pred_names — character vector of predictor names (in order)
#   p          — integer, lag order of the VAR
#
# Returns:
#   numeric scalar — the nowcast for GDP QoQ growth
# ------------------------------------------------------------

bridge_nowcast <- function(var_mod, train, pred_avgs,
                           pred_names, p) {
  coefs <- coef(var_mod)$gdp[, "Estimate"]

  # Last p values of GDP for the lagged dependent variable terms
  gdp_lags <- tail(train$gdp, p)

  # Build regressor vector in VAR coefficient order:
  # gdp.l1, pred1.l1, ..., predN.l1, gdp.l2, ..., const, covid
  x <- c()
  for (lag in 1:p) {
    x <- c(x, gdp_lags[p - lag + 1])
    for (nm in pred_names) {
      x <- c(x, train[[nm]][nrow(train) - lag + 1])
    }
  }
  x <- c(x, 1)  # constant term
  if ("covid" %in% names(coefs)) x <- c(x, 0)  # not a COVID quarter

  # Bridge step: replace lag-1 predictor values with the
  # current-quarter partial averages. This conditions the
  # forecast on the latest available monthly information.
  if (!is.null(pred_avgs) && ncol(pred_avgs) == length(pred_names)) {
    for (j in seq_along(pred_names)) {
      nm  <- pred_names[j]
      idx <- 1 + j  # position of pred.l1 (after gdp.l1)
      if (!is.na(pred_avgs[[nm]])) x[idx] <- pred_avgs[[nm]]
    }
  }

  names(x) <- names(coefs)[seq_along(x)]
  sum(coefs[seq_along(x)] * x)
}


# ------------------------------------------------------------
# compute_metrics()
# Compute RMSE, MAE, and relative RMSE vs AR(1) benchmark
# for all nowcast vintages in an evaluation data frame.
#
# Arguments:
#   df     — data.frame with columns: actual, var_m1, var_m2,
#            var_m3, ar1 (as produced by the evaluation loop)
#   label  — character, label printed in the output header
#
# Prints results to console and invisibly returns a data.frame
# with the metrics for each vintage.
# ------------------------------------------------------------

compute_metrics <- function(df, label) {
  cat(sprintf("\n--- %s (n = %d) ---\n", label, nrow(df)))

  vintages <- c("var_m1", "var_m2", "var_m3", "ar1")
  out <- lapply(vintages, function(v) {
    e    <- df$actual - df[[v]]
    rmse <- sqrt(mean(e^2, na.rm = TRUE))
    mae  <- mean(abs(e),   na.rm = TRUE)
    n    <- sum(!is.na(e))
    cat(sprintf("  %-8s  RMSE = %.4f  MAE = %.4f  n = %d\n",
                v, rmse, mae, n))
    data.frame(vintage = v, rmse = rmse, mae = mae, n = n)
  })

  rmse_var <- sqrt(mean((df$actual - df$var_m3)^2, na.rm = TRUE))
  rmse_ar  <- sqrt(mean((df$actual - df$ar1)^2,    na.rm = TRUE))
  rel      <- rmse_var / rmse_ar
  cat(sprintf("  Relative RMSE (VAR M3 / AR1): %.3f  %s\n",
              rel,
              ifelse(rel < 1,
                     "VAR beats benchmark \u2713",
                     "VAR does not beat benchmark")))

  invisible(dplyr::bind_rows(out) %>%
              dplyr::mutate(rel_rmse = rmse / rmse_ar))
}


# ------------------------------------------------------------
# fit_var()
# Estimate a VAR with AIC lag selection and COVID dummy.
# Wraps VARselect + VAR from the vars package with error
# handling, for use inside recursive evaluation loops.
#
# Arguments:
#   train     — data.frame, estimation sample
#   pred_names — character vector of predictor names
#   lag_max   — integer, maximum lags to consider (default 2)
#
# Returns:
#   list with elements:
#     model  — VAR object (or NULL if estimation failed)
#     p      — integer, selected lag order
# ------------------------------------------------------------

fit_var <- function(train, pred_names, lag_max = 2) {
  var_mat <- train %>%
    dplyr::select(gdp, dplyr::all_of(pred_names)) %>%
    as.matrix()

  exog <- matrix(train$covid, ncol = 1,
                 dimnames = list(NULL, "covid"))

  p <- tryCatch(
    vars::VARselect(var_mat, lag.max = lag_max,
                    type = "const",
                    exogen = exog)$selection["AIC(n)"],
    error = function(e) 1L
  )

  model <- tryCatch(
    vars::VAR(var_mat, p = p, type = "const", exogen = exog),
    error = function(e) NULL
  )

  list(model = model, p = as.integer(p))
}


# ------------------------------------------------------------
# ar1_forecast()
# Compute a one-step-ahead AR(1) forecast for GDP growth.
# Used as the naive benchmark in forecast evaluation.
#
# Arguments:
#   gdp_series — numeric vector, GDP growth training data
#
# Returns:
#   numeric scalar — AR(1) point forecast
# ------------------------------------------------------------

ar1_forecast <- function(gdp_series) {
  mod <- ar(gdp_series, order.max = 1, method = "ols", aic = FALSE)
  as.numeric(tail(gdp_series, 1) * mod$ar +
             mod$x.mean * (1 - mod$ar))
}
