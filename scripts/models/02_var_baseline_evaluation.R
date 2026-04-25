# ============================================================================
# scripts/models/02_var_baseline_evaluation.R
#
# PURPOSE:
#   Evaluate the baseline VAR models from 01_var_baseline.R against
#   univariate AR(p) benchmarks using pseudo out-of-sample expanding-window
#   RMSE. Includes formal Diebold-Mariano tests for forecast superiority.
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")    — data + globals
#   source("scripts/models/01_var_baseline.R")  — results, var_input, samples,
#                                                  subset_sample(), make_ts(),
#                                                  covid_cols, exog_full
#
# STRUCTURE:
#   Section 0 — Shared evaluation parameters
#   Section 1 — AR(p) benchmark: in-sample estimates (p = 1..5, all samples)
#   Section 2 — Expanding-window pseudo out-of-sample forecast loop
#               runs AR(1) + all VAR baseline models
#   Section 3 — RMSE tables (with COVID / without COVID)
#   Section 4 — Relative RMSE table vs AR(1) benchmark
#   Section 5 — Diebold-Mariano tests: best VAR vs AR(1)
#   Section 6 — RMSE plots by horizon
#
# EVALUATION DESIGN:
#   Scheme        : Expanding window (fixed start, growing sample)
#   First forecast: 2010 Q1  (~10 years of evaluation, matching meeting notes)
#   Last forecast : 2025 Q3  (h=1 forecast; shorter for longer horizons)
#   Horizons      : h = 1, 2, 4, 8 quarters
#   Variables     : gdp_g, cpi_g, bond_dif  (all three target variables)
#   COVID treatment:
#     "incl" — error calculated over ALL evaluation periods
#     "excl" — error calculated excluding 2020Q1–2021Q4 from the denominator
#              (the COVID dummies handle them in estimation, but they still
#               dominate RMSE; excluding gives a cleaner view of normal performance)
#
# FINDINGS (fill in after running):
#   AR(1) RMSE GDP   h=1: [  ]   h=4: [  ]   h=8: [  ]
#   Best VAR GDP     h=1: [  ]   h=4: [  ]   h=8: [  ]   (model: [  ])
#   Relative RMSE best VAR vs AR(1) GDP h=1: [  ]%
#   DM test best VAR vs AR(1) GDP h=1: p=[  ]
#   With vs without COVID: RMSE difference for GDP h=1: [  ] pp
#
# Outputs (commented — uncomment to save):
#   output/tables/02_rmse_table_incl_covid.csv
#   output/tables/02_rmse_table_excl_covid.csv
#   output/tables/02_rmse_relative.csv
#   output/tables/02_dm_tests.csv
#   output/figures/02_rmse_by_horizon.png
# ============================================================================

library(lmtest)   # grangertest, already loaded but being explicit
library(forecast) # Acf, auto.arima — install if needed: install.packages("forecast")

# ── 0. Shared evaluation parameters ──────────────────────────────────────────
EVAL_START   <- as.Date("2010-01-01")   # first quarter we forecast
HORIZONS     <- c(1, 2, 4, 8)          # forecast horizons to evaluate
N_AHEAD_MAX  <- max(HORIZONS)          # max horizon needed in predict()

# COVID quarters to exclude from RMSE calculation (not from estimation)
COVID_EVAL_EXCL <- seq(as.Date("2020-01-01"),
                       as.Date("2021-10-01"), by = "quarter")

# All evaluation dates in var_input
eval_dates <- var_input$date[var_input$date >= EVAL_START]
cat("Evaluation period:", format(min(eval_dates)), "–",
    format(max(eval_dates)), "|", length(eval_dates), "quarters\n")
cat("Horizons:", paste(HORIZONS, collapse = ", "), "\n")
cat("COVID quarters excluded from error calc:",
    length(COVID_EVAL_EXCL), "quarters\n\n")

# ── 1. AR(p) benchmark — in-sample estimates (for reference) ─────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 1: AR(p) IN-SAMPLE ESTIMATES\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# We fit AR(1) through AR(5) for each variable × sample
# AR model: y_t = c + a1*y_{t-1} + ... + ap*y_{t-p} + e_t
# No COVID dummies here — AR is the naive benchmark; adding dummies would
# inflate its apparent fit and make comparisons unfair.
# (We do include them in the expanding-window AR forecasts below.)

ar_vars   <- c("gdp_g", "cpi_g", "bond_dif")
ar_labels <- c("GDP QoQ %", "CPI QoQ %", "Bond yield QoQ pp")

ar_insample <- list()

for (vname in ar_vars) {
  cat(sprintf("── %s ──────────────────────────────\n", vname))
  ar_insample[[vname]] <- list()
  
  for (sname in names(samples)) {
    s_data <- subset_sample(samples[[sname]]$start_date)
    y      <- s_data$dat[[vname]]
    n      <- length(y)
    
    cat(sprintf("  Sample: %-8s [n=%d]\n", sname, n))
    for (p in 1:5) {
      if (n < p + 10) next   # skip if too few obs
      fit <- tryCatch(
        arima(y, order = c(p, 0, 0), include.mean = TRUE),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      sigma <- sqrt(fit$sigma2)
      aic   <- AIC(fit)
      cat(sprintf("    AR(%d): sigma=%.4f  AIC=%.2f\n", p, sigma, aic))
      ar_insample[[vname]][[paste(sname, p, sep="_")]] <- list(
        fit = fit, sigma = sigma, aic = aic, p = p, sample = sname
      )
    }
  }
  cat("\n")
}

# ── 2. Expanding-window forecast loop ─────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 2: EXPANDING-WINDOW PSEUDO OUT-OF-SAMPLE FORECASTS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("This may take a minute — forecasting", length(eval_dates),
    "windows × models × horizons...\n\n")

# Storage: forecast_errors[[model_id]][[variable]] = tibble(date, h, error)
forecast_errors <- list()

# ── Helper: AR(1) expanding window forecast ───────────────────────────────────
# Returns tibble with columns: date (forecast target), h, fcst, actual, error
ar_expanding <- function(y_full, dates_full, eval_start,
                         horizons, p_ar = 1, use_dummies = FALSE,
                         exog_full_arg = NULL) {
  rows <- list()
  eval_idx <- which(dates_full >= eval_start)
  
  for (i in eval_idx) {
    for (h in horizons) {
      target_i <- i + h - 1
      if (target_i > length(y_full)) next
      
      # Training data: all obs before the forecast origin
      y_train <- y_full[1:(i - 1)]
      if (length(y_train) < p_ar + 5) next
      
      # Fit AR(p)
      fit <- tryCatch(
        arima(y_train, order = c(p_ar, 0, 0), include.mean = TRUE),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      
      # Forecast h steps ahead
      fc <- tryCatch(
        predict(fit, n.ahead = h),
        error = function(e) NULL
      )
      if (is.null(fc)) next
      
      actual <- y_full[target_i]
      fcst   <- as.numeric(fc$pred[h])
      rows[[length(rows) + 1]] <- tibble(
        date   = dates_full[target_i],
        origin = dates_full[i],
        h      = h,
        fcst   = fcst,
        actual = actual,
        error  = actual - fcst
      )
    }
  }
  bind_rows(rows)
}

# ── Helper: VAR expanding window forecast ────────────────────────────────────
var_expanding <- function(sname, p_var, eval_start, horizons) {
  s       <- samples[[sname]]
  rows    <- list()
  
  # All dates from sample start
  all_idx <- which(var_input$date >= s$start_date)
  eval_idx <- which(var_input$date >= eval_start)
  # Only evaluate within the sample
  eval_idx <- intersect(eval_idx, all_idx)
  
  for (i in eval_idx) {
    for (h in horizons) {
      target_i <- i + h - 1
      if (target_i > nrow(var_input)) next
      
      # Training: from sample start to i-1
      train_idx <- all_idx[all_idx < i]
      if (length(train_idx) < p_var * 3 + 15) next
      
      dat_train  <- var_input[train_idx, ]
      exog_train <- exog_full[train_idx, , drop = FALSE]
      n_train    <- nrow(dat_train)
      
      # Build ts
      ts_start_q <- as.numeric(format(dat_train$date[1], "%Y"))
      ts_start_m <- as.numeric(format(dat_train$date[1], "%m")) %/% 3 + 1
      en <- ts(dat_train[, c("gdp_g","cpi_g","bond_dif")],
               start = c(ts_start_q, ts_start_m), frequency = 4)
      
      # Fit VAR
      mod <- tryCatch(
        VAR(en, p = p_var, type = "const", exogen = exog_train),
        error = function(e) NULL
      )
      if (is.null(mod)) next
      
      # Forecast exog (all zeros — no COVID in forecast period)
      exog_fc_h <- matrix(0, nrow = max(horizons), ncol = 8,
                          dimnames = list(NULL, covid_cols))
      
      fc <- tryCatch(
        predict(mod, n.ahead = max(horizons), ci = 0.95,
                dumvar = exog_fc_h),
        error = function(e) NULL
      )
      if (is.null(fc)) next
      
      # Extract h-step forecast for each variable
      for (vname in c("gdp_g","cpi_g","bond_dif")) {
        actual <- var_input[[vname]][target_i]
        fcst   <- fc$fcst[[vname]][h, "fcst"]
        rows[[length(rows) + 1]] <- tibble(
          date   = var_input$date[target_i],
          origin = var_input$date[i],
          h      = h,
          variable = vname,
          fcst   = fcst,
          actual = actual,
          error  = actual - fcst
        )
      }
    }
  }
  bind_rows(rows)
}

# ── Run AR(1) benchmarks ──────────────────────────────────────────────────────
cat("Fitting AR(1) benchmarks...\n")
for (vname in ar_vars) {
  cat(sprintf("  AR(1) — %s\n", vname))
  fc_ar <- ar_expanding(
    y_full     = var_input[[vname]],
    dates_full = var_input$date,
    eval_start = EVAL_START,
    horizons   = HORIZONS,
    p_ar       = 1
  ) %>% mutate(variable = vname)
  
  forecast_errors[[paste0("AR1_", vname)]] <- fc_ar
}

# ── Run VAR baseline models ───────────────────────────────────────────────────
# We evaluate all 9 estimated models from 01_var_baseline.R
# plus the 3 best models identified there
var_models_to_eval <- list(
  list(sname = "full",   p = 1),
  list(sname = "full",   p = 2),
  list(sname = "full",   p = 3),
  list(sname = "full",   p = 4),
  list(sname = "post08", p = 1),
  list(sname = "post08", p = 2),
  list(sname = "post08", p = 3),
  list(sname = "post15", p = 1),
  list(sname = "post15", p = 2)
)

for (spec in var_models_to_eval) {
  id <- sprintf("VAR_%s_p%d", spec$sname, spec$p)
  cat(sprintf("  VAR — %s\n", id))
  
  fc_var <- tryCatch(
    var_expanding(sname = spec$sname, p_var = spec$p,
                  eval_start = EVAL_START, horizons = HORIZONS),
    error = function(e) {
      message("  ERROR in ", id, ": ", e$message)
      NULL
    }
  )
  
  if (!is.null(fc_var) && nrow(fc_var) > 0) {
    forecast_errors[[id]] <- fc_var
  }
}

cat("\nForecast loop complete.\n")
cat("Models evaluated:", length(forecast_errors), "\n\n")

# ── 3. RMSE tables ────────────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 3: RMSE TABLES\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Helper: compute RMSE from error tibble, optionally excluding COVID quarters
compute_rmse <- function(err_df, excl_covid = FALSE) {
  if (excl_covid) {
    err_df <- err_df %>% filter(!date %in% COVID_EVAL_EXCL)
  }
  err_df %>%
    group_by(variable, h) %>%
    summarise(
      n      = n(),
      rmse   = sqrt(mean(error^2, na.rm = TRUE)),
      mae    = mean(abs(error),   na.rm = TRUE),
      .groups = "drop"
    )
}

# Combine all VAR forecasts (they already have variable column)
# AR forecasts need to be tagged by variable
all_fc <- bind_rows(
  # AR1 benchmarks
  bind_rows(lapply(ar_vars, function(v) {
    forecast_errors[[paste0("AR1_", v)]] %>%
      mutate(model = "AR1")
  })),
  # VAR models
  bind_rows(lapply(names(forecast_errors)[grepl("^VAR_", names(forecast_errors))],
                   function(mid) {
                     forecast_errors[[mid]] %>% mutate(model = mid)
                   }
  ))
)

# RMSE including COVID
rmse_incl <- all_fc %>%
  group_by(model, variable, h) %>%
  summarise(
    n    = n(),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    mae  = mean(abs(error),   na.rm = TRUE),
    .groups = "drop"
  )

# RMSE excluding COVID quarters from error calculation
rmse_excl <- all_fc %>%
  filter(!date %in% COVID_EVAL_EXCL) %>%
  group_by(model, variable, h) %>%
  summarise(
    n    = n(),
    rmse = sqrt(mean(error^2, na.rm = TRUE)),
    mae  = mean(abs(error),   na.rm = TRUE),
    .groups = "drop"
  )

# Print RMSE for GDP (primary target)
cat("── GDP RMSE (including COVID) ──────────────────────────────────────\n")
rmse_incl %>%
  filter(variable == "gdp_g") %>%
  select(model, h, rmse) %>%
  pivot_wider(names_from = h, values_from = rmse,
              names_prefix = "h=") %>%
  arrange(model) %>%
  mutate(across(starts_with("h="), ~ round(.x, 4))) %>%
  print(n = Inf, width = Inf)

cat("\n── GDP RMSE (excluding COVID 2020Q1–2021Q4) ────────────────────────\n")
rmse_excl %>%
  filter(variable == "gdp_g") %>%
  select(model, h, rmse) %>%
  pivot_wider(names_from = h, values_from = rmse,
              names_prefix = "h=") %>%
  arrange(model) %>%
  mutate(across(starts_with("h="), ~ round(.x, 4))) %>%
  print(n = Inf, width = Inf)

cat("\n── CPI RMSE (including COVID) ──────────────────────────────────────\n")
rmse_incl %>%
  filter(variable == "cpi_g") %>%
  select(model, h, rmse) %>%
  pivot_wider(names_from = h, values_from = rmse,
              names_prefix = "h=") %>%
  arrange(model) %>%
  mutate(across(starts_with("h="), ~ round(.x, 4))) %>%
  print(n = Inf, width = Inf)

cat("\n── Bond yield RMSE (including COVID) ───────────────────────────────\n")
rmse_incl %>%
  filter(variable == "bond_dif") %>%
  select(model, h, rmse) %>%
  pivot_wider(names_from = h, values_from = rmse,
              names_prefix = "h=") %>%
  arrange(model) %>%
  mutate(across(starts_with("h="), ~ round(.x, 4))) %>%
  print(n = Inf, width = Inf)

# ── 4. Relative RMSE vs AR(1) ─────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 4: RELATIVE RMSE vs AR(1) (AR(1) = 100)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("Values < 100 mean the model beats AR(1). Bold improvement > 10%.\n\n")

for (covid_label in c("incl", "excl")) {
  rmse_tbl <- if (covid_label == "incl") rmse_incl else rmse_excl
  cat(sprintf("── %s COVID ─────────────────────────────────────────────────\n",
              ifelse(covid_label == "incl", "Including", "Excluding")))
  
  for (vname in ar_vars) {
    # AR1 benchmark RMSE for this variable
    ar1_rmse <- rmse_tbl %>%
      filter(model == "AR1", variable == vname) %>%
      select(h, ar1_rmse = rmse)
    
    rel <- rmse_tbl %>%
      filter(variable == vname, model != "AR1") %>%
      left_join(ar1_rmse, by = "h") %>%
      mutate(rel_rmse = round(rmse / ar1_rmse * 100, 1)) %>%
      select(model, h, rel_rmse) %>%
      pivot_wider(names_from = h, values_from = rel_rmse,
                  names_prefix = "h=") %>%
      arrange(model)
    
    cat(sprintf("\n  %s:\n", vname))
    print(rel, n = Inf, width = Inf)
  }
  cat("\n")
}

# ── 5. Diebold-Mariano tests ──────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 5: DIEBOLD-MARIANO TESTS\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
cat("H0: Model and AR(1) have equal predictive accuracy\n")
cat("Reject H0 (p<0.05) → model is significantly better (or worse) than AR(1)\n\n")

# Helper: DM test for one model vs AR(1) for one variable and horizon
dm_test <- function(err_model, err_ar1, h_val) {
  # Align on dates
  joined <- inner_join(
    err_model %>% filter(h == h_val) %>% select(date, e_model = error),
    err_ar1   %>% filter(h == h_val) %>% select(date, e_ar1   = error),
    by = "date"
  ) %>% filter(!is.na(e_model), !is.na(e_ar1))
  
  if (nrow(joined) < 10) return(tibble(dm_stat = NA, p_value = NA, n = nrow(joined)))
  
  d <- joined$e_model^2 - joined$e_ar1^2   # loss differential (MSE)
  tryCatch({
    # DM test: t-test on d with HAC standard errors
    dm    <- lm(d ~ 1)
    dm_hc <- coeftest(dm, vcov = sandwich::NeweyWest(dm, lag = h_val - 1))
    tibble(
      dm_stat = round(dm_hc[1, "t value"], 3),
      p_value = round(dm_hc[1, "Pr(>|t|)"], 4),
      n       = nrow(joined)
    )
  }, error = function(e) tibble(dm_stat = NA, p_value = NA, n = nrow(joined)))
}

# Run DM tests for best VAR models vs AR(1) at all horizons
best_var_ids <- c("VAR_full_p1", "VAR_post08_p3", "VAR_post15_p1")

dm_results <- list()
for (mid in best_var_ids) {
  if (is.null(forecast_errors[[mid]])) next
  for (vname in ar_vars) {
    err_var <- forecast_errors[[mid]] %>% filter(variable == vname)
    err_ar1 <- forecast_errors[[paste0("AR1_", vname)]]
    for (h_val in HORIZONS) {
      dm_res <- dm_test(err_var, err_ar1, h_val)
      dm_results[[length(dm_results) + 1]] <- tibble(
        model    = mid,
        variable = vname,
        h        = h_val
      ) %>% bind_cols(dm_res)
    }
  }
}

dm_table <- bind_rows(dm_results) %>%
  mutate(
    significant = case_when(
      is.na(p_value)   ~ "n/a",
      p_value <= 0.01  ~ "*** (1%)",
      p_value <= 0.05  ~ "**  (5%)",
      p_value <= 0.10  ~ "*   (10%)",
      TRUE             ~ "—"
    ),
    direction = case_when(
      is.na(dm_stat)  ~ "n/a",
      dm_stat < 0     ~ "VAR better",
      dm_stat > 0     ~ "AR1 better",
      TRUE            ~ "—"
    )
  )

cat("── GDP ─────────────────────────────────────────────────────────────\n")
dm_table %>% filter(variable == "gdp_g") %>%
  select(model, h, dm_stat, p_value, significant, direction) %>%
  print(n = Inf, width = Inf)

cat("\n── CPI ─────────────────────────────────────────────────────────────\n")
dm_table %>% filter(variable == "cpi_g") %>%
  select(model, h, dm_stat, p_value, significant, direction) %>%
  print(n = Inf, width = Inf)

cat("\n── Bond yield ──────────────────────────────────────────────────────\n")
dm_table %>% filter(variable == "bond_dif") %>%
  select(model, h, dm_stat, p_value, significant, direction) %>%
  print(n = Inf, width = Inf)

# ── 6. RMSE plots by horizon ──────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("SECTION 6: RMSE PLOTS BY HORIZON\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Colour palette for models
model_colours <- c(
  "AR1"            = "black",
  "VAR_full_p1"    = COL_GDP,
  "VAR_full_p2"    = "#4393c3",
  "VAR_full_p3"    = "#92c5de",
  "VAR_full_p4"    = "#d1e5f0",
  "VAR_post08_p1"  = COL_CPI,
  "VAR_post08_p2"  = "#f4a582",
  "VAR_post08_p3"  = "#b2182b",
  "VAR_post15_p1"  = COL_BOND,
  "VAR_post15_p2"  = "#4dac26"
)

make_rmse_plot <- function(rmse_tbl, vname, var_label,
                           covid_label = "incl") {
  plot_df <- rmse_tbl %>%
    filter(variable == vname, h %in% HORIZONS) %>%
    mutate(h = factor(h, levels = HORIZONS))
  
  ggplot(plot_df, aes(x = h, y = rmse,
                      colour = model, group = model)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = model_colours,
                        na.value = "grey70", name = NULL) +
    labs(
      title    = sprintf("%s — RMSE by forecast horizon", var_label),
      subtitle = sprintf("%s COVID quarters | Eval: 2010 Q1 – 2025",
                         ifelse(covid_label == "incl", "Including", "Excluding")),
      x = "Forecast horizon (quarters)",
      y = "RMSE"
    ) +
    theme(legend.position = "right",
          legend.text     = element_text(size = 8))
}

# Plot for each variable, including and excluding COVID
for (covid_label in c("incl", "excl")) {
  rmse_tbl <- if (covid_label == "incl") rmse_incl else rmse_excl
  
  p_gdp  <- make_rmse_plot(rmse_tbl, "gdp_g",    "GDP growth",         covid_label)
  p_cpi  <- make_rmse_plot(rmse_tbl, "cpi_g",    "CPI inflation",      covid_label)
  p_bond <- make_rmse_plot(rmse_tbl, "bond_dif", "Bond yield change",  covid_label)
  
  fig <- (p_gdp / p_cpi / p_bond) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = sprintf("Baseline VAR vs AR(1): RMSE by horizon (%s COVID)",
                         ifelse(covid_label == "incl", "including", "excluding")),
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )
  
  print(fig)
  
  # # Save (uncomment)
  # ggsave(
  #   file.path(here("output","figures"),
  #             sprintf("02_rmse_by_horizon_%s_covid.png", covid_label)),
  #   fig, width = 10, height = 12, dpi = 150
  # )
}

# ── Save tables (uncomment to save) ──────────────────────────────────────────
# write.csv(rmse_incl,
#           file.path(here("output","tables"), "02_rmse_incl_covid.csv"),
#           row.names = FALSE)
# write.csv(rmse_excl,
#           file.path(here("output","tables"), "02_rmse_excl_covid.csv"),
#           row.names = FALSE)
# write.csv(dm_table,
#           file.path(here("output","tables"), "02_dm_tests.csv"),
#           row.names = FALSE)

cat("\nEvaluation complete.\n")