# ============================================================================
# scripts/models/bvar/03_bvar_baseline_diagnostics.R
#
# PURPOSE:
#   Diagnostics for the BVAR baseline models from 01_bvar_baseline.R.
#   Focuses on what BVAR diagnostics add over and above the frequentist VAR:
#     1. IRF with genuine posterior uncertainty bands (not bootstrapped CI)
#     2. FEVD with posterior distributions — not just point estimates
#     3. MCMC convergence: trace plots and acceptance rates
#     4. Residual checks via posterior predictive approach
#     5. Prior vs posterior comparison: how much does the data update the prior?
#
# REQUIRES:
#   source("scripts/exploration/00_setup.R")
#   source("scripts/models/var/01_var_baseline.R")
#   source("scripts/models/bvar/01_bvar_baseline.R")
#
# FOCUS MODELS:
#   One representative model per prior (best per prior from script 01),
#   all of which are full_p4. Then separately bvar_minnesota_post08_p1
#   for residual checks (mirrors the frequentist post08_p3).
#
#   Primary:   bvar_minnesota_full_p4    — best Minnesota
#              bvar_normalwishart_full_p4 — best Normal-Wishart
#              bvar_dummyobs_full_p4     — best Dummy observations
#   Secondary: bvar_minnesota_post08_p1  — comparable to frequentist post08_p3
#
# NOTE ON ARRAY DIMENSIONS (verified via dim()):
#   bvar_irf$irf  : [draw, response, horizon, impulse]
#   bvar_fevd$fevd: [draw, response, horizon, impulse]
#   Reduction over draws always uses apply(..., c(2, 3, 4), ...)
#
# NOTE ON MCMC ACCEPTANCE RATE:
#   mod$acc_rate does not exist in BVAR v1.0.5
#   Correct field: mod$meta$accepted / mod$meta$n_save
#   Minnesota / Normal-Wishart use Gibbs sampler — no MH step, rate N/A
#   Dummy obs uses MH for hyperparameters — rate meaningful
#
# KEY FINDINGS:
#
#   MCMC convergence:
#     Minnesota / Normal-Wishart: Gibbs sampler, no MH rejection step
#     Dummy obs: MH acceptance rate 33.1% ✓ (target: 25-45%)
#     Trace plots: all chains well-mixed, no drift or stuck regions
#
#   IRF (68% posterior bands, full sample, p=4):
#     GDP -> CPI:  excludes zero at h=1,2,3 across ALL priors
#                  (median: +0.12 at h=1, decaying to +0.05 at h=3)
#     CPI -> Bond: excludes zero at h=1 across ALL priors (median: ~+0.236)
#     Bond -> CPI: excludes zero at h=2,4,5,6 across ALL priors
#                  (median: +0.03 to +0.08)
#     Bond -> GDP: excludes zero at h=2,4,7 across ALL priors
#     Prior choice barely affects IRF medians or CI coverage —
#     all three priors yield virtually identical impulse responses
#
#   FEVD (posterior median, full sample, p=4):
#     GDP:  own shocks ~89-91% at h=8 | bond explains ~6-7%
#     CPI:  own shocks ~73% at h=8    | GDP explains ~16% | bond ~9%
#     Bond: own shocks ~80% at h=8    | CPI explains ~14% | GDP ~6%
#     Results virtually identical across all three priors
#     CPI explains ~14% of bond variance at all horizons —
#     inflation expectations persistently priced into Swiss long rates
#
#   Residuals (bvar_minnesota_post08_p1):
#     GDP:  non-normal (JB p<0.001) — COVID outliers visible in QQ tails
#           expected given Approach B OLS-demeaning (vs exog dummies in VAR)
#     CPI:  normal ✓ (JB p=0.990) | no autocorrelation ✓ (LB p=0.784)
#     Bond: normal ✓ (JB p=0.077) | no autocorrelation ✓ (LB p=0.990)
#
#   Prior vs posterior — lambda (prior mode = 0.2):
#     Full sample:  0.54-0.55 across priors
#     Post-GFC:     0.59-0.61
#     Post-2015:    0.63-0.70
#     Data strongly updates lambda upward across all 27 models —
#     Swiss macro data is informative and resists random-walk shrinkage
#     Shorter samples push lambda higher: less data, more uncertainty,
#     yet data still dominates the prior
#
#   Prior vs posterior — dummy obs hyperparameters (prior mode = 1.0):
#     soc posterior: 1.42-1.86 across all models
#       -> Data prefers MORE persistence than prior assumed
#     sur posterior: 0.57-1.09 across all models
#       -> Data roughly consistent with prior on common trend component
#
# Outputs (commented — uncomment to save):
#   output/figures/bvar_03_irf_comparison.png
#   output/figures/bvar_03_fevd_<prior>.png
#   output/figures/bvar_03_mcmc_<prior>.png
#   output/figures/bvar_03_prior_posterior.png
# ============================================================================

library(dplyr)

IRF_HORIZON <- 12
TARGET_VARS <- c("gdp_g", "cpi_g", "bond_dif")
VAR_NAMES   <- c("GDP growth", "CPI inflation", "Bond yield \u0394")

# ── helpers: reduce over draw dimension (dim 1) ───────────────────────────────
# Input:  [draw, response, horizon, impulse]
# Output: [response, horizon, impulse]
irf_median <- function(arr) {
  apply(arr, c(2, 3, 4), median, na.rm = TRUE)
}
irf_quantile <- function(arr, probs) {
  apply(arr, c(2, 3, 4), quantile, probs = probs, na.rm = TRUE)
}

# =============================================================================
# SECTION 1a: IRF PER PRIOR
# =============================================================================
cat("=== SECTION 1a: IRF WITH 68% POSTERIOR BANDS ===\n\n")

for (prior_name in names(PRIORS)) {
  best_row <- best_per_prior %>% filter(prior == prior_name)
  res      <- bvar_results[[prior_name]][[best_row$sample]][[
    paste0("p", best_row$p)]]
  if (is.null(res)) next
  
  bvar_irf <- tryCatch(
    irf(res$model, horizon = IRF_HORIZON, conf_bands = c(0.16, 0.84)),
    error = function(e) {
      cat("IRF error:", prior_name, "-", e$message, "\n")
      NULL
    }
  )
  if (is.null(bvar_irf)) next
  
  # Built-in BVAR plot (3x3 grid)
  plot(bvar_irf, area = TRUE,
       col = c(COL_GDP, COL_CPI, COL_BOND),
       mar = c(2, 2, 2, 0.5))
  title(
    main  = sprintf("IRF — %s | %s | 68%% CI", res$id, PRIOR_LABELS[prior_name]),
    outer = TRUE, line = -1, cex.main = 0.85
  )
  
  # Print which off-diagonal channels exclude zero
  tryCatch({
    irf_q16 <- irf_quantile(bvar_irf$irf, 0.16)
    irf_q84 <- irf_quantile(bvar_irf$irf, 0.84)
    irf_med <- irf_median(bvar_irf$irf)
    # all: [response, horizon, impulse]
    vn        <- c("GDP", "CPI", "Bond")
    found_any <- FALSE
    
    cat(sprintf("[%s] 68%% CI excludes zero:\n", prior_name))
    for (imp in 1:3) {
      for (resp in 1:3) {
        if (imp == resp) next
        excl_h <- which(irf_q16[resp, , imp] > 0 | irf_q84[resp, , imp] < 0)
        if (length(excl_h) > 0) {
          found_any <- TRUE
          cat(sprintf("  %s -> %s: h=%s (med: %s)\n",
                      vn[imp], vn[resp],
                      paste(excl_h, collapse = ","),
                      paste(sprintf("%+.3f", irf_med[resp, excl_h, imp]),
                            collapse = ",")))
        }
      }
    }
    if (!found_any) cat("  (none exclude zero at 68%%)\n")
    cat("\n")
  }, error = function(e) cat("  CI extraction error:", e$message, "\n\n"))
}

# =============================================================================
# SECTION 1b: IRF OVERLAY COMPARISON ACROSS PRIORS
# =============================================================================
cat("=== SECTION 1b: IRF COMPARISON ACROSS PRIORS ===\n\n")

irf_compare <- bind_rows(lapply(names(PRIORS), function(prior_name) {
  best_row <- best_per_prior %>% filter(prior == prior_name)
  res      <- bvar_results[[prior_name]][[best_row$sample]][[
    paste0("p", best_row$p)]]
  if (is.null(res)) return(NULL)
  
  bvar_irf <- tryCatch(
    irf(res$model, horizon = IRF_HORIZON, conf_bands = c(0.16, 0.84)),
    error = function(e) NULL
  )
  if (is.null(bvar_irf)) return(NULL)
  
  irf_med <- irf_median(bvar_irf$irf)
  irf_q16 <- irf_quantile(bvar_irf$irf, 0.16)
  irf_q84 <- irf_quantile(bvar_irf$irf, 0.84)
  # all: [response, horizon, impulse]
  
  vnl <- c("GDP", "CPI", "Bond")
  bind_rows(lapply(1:3, function(imp) {
    bind_rows(lapply(1:3, function(resp) {
      if (imp == resp) return(NULL)
      tibble(
        prior    = prior_name,
        impulse  = vnl[imp],
        response = vnl[resp],
        h        = 1:IRF_HORIZON,
        med      = irf_med[resp, , imp],
        q16      = irf_q16[resp, , imp],
        q84      = irf_q84[resp, , imp]
      )
    }))
  }))
}))

if (!is.null(irf_compare) && nrow(irf_compare) > 0) {
  prior_cols <- c(minnesota     = COL_GDP,
                  normalwishart = COL_CPI,
                  dummyobs      = COL_BOND)
  
  channels <- irf_compare %>%
    distinct(impulse, response) %>%
    mutate(panel_title = sprintf("%s -> %s", impulse, response))
  
  panels_irf <- lapply(seq_len(nrow(channels)), function(i) {
    ch  <- channels[i, ]
    dat <- irf_compare %>%
      filter(impulse == ch$impulse, response == ch$response)
    ggplot(dat, aes(x = h, colour = prior, fill = prior)) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 color = COL_GREY, linewidth = 0.35) +
      geom_ribbon(aes(ymin = q16, ymax = q84), alpha = 0.12, colour = NA) +
      geom_line(aes(y = med), linewidth = 0.8) +
      scale_colour_manual(values = prior_cols,
                          labels = PRIOR_LABELS, name = NULL) +
      scale_fill_manual(values = prior_cols,
                        labels = PRIOR_LABELS, name = NULL) +
      scale_x_continuous(breaks = seq(0, IRF_HORIZON, by = 4)) +
      labs(title = ch$panel_title, x = "Quarters", y = "Response") +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, face = "bold"))
  })
  
  legend_panel <- ggplot(
    tibble(
      prior = names(prior_cols),
      label = unname(PRIOR_LABELS[names(prior_cols)]),
      x     = 1,
      y     = seq_along(prior_cols)
    ),
    aes(x = x, y = y, colour = prior, label = label)
  ) +
    geom_point(size = 3) +
    geom_text(aes(x = 1.1), hjust = 0, size = 3.2) +
    scale_colour_manual(values = prior_cols) +
    xlim(0.9, 3) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(size = 9, face = "bold")) +
    labs(title = "Legend")
  
  fig_irf_comp <- wrap_plots(c(panels_irf, list(legend_panel)),
                             nrow = 2, ncol = 4) +
    plot_annotation(
      title    = "BVAR IRF — all three priors | 68% posterior bands",
      subtitle = "Full sample | p=4",
      theme    = theme(
        plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8,  color = "gray40")
      )
    )
  print(fig_irf_comp)
  cat("IRF comparison plot printed.\n\n")
  # ggsave(file.path(here("output","figures"), "bvar_03_irf_comparison.png"),
  #        fig_irf_comp, width = 14, height = 7, dpi = 150)
}

# =============================================================================
# SECTION 2: FEVD
# =============================================================================
cat("=== SECTION 2: FEVD ===\n\n")

for (prior_name in names(PRIORS)) {
  best_row <- best_per_prior %>% filter(prior == prior_name)
  res      <- bvar_results[[prior_name]][[best_row$sample]][[
    paste0("p", best_row$p)]]
  if (is.null(res)) next
  
  bvar_fevd <- tryCatch(
    fevd(res$model, horizon = IRF_HORIZON),
    error = function(e) {
      cat("FEVD error:", prior_name, "-", e$message, "\n")
      NULL
    }
  )
  if (is.null(bvar_fevd)) next
  
  fevd_arr <- bvar_fevd$fevd
  cat(sprintf("[%s] dim(fevd$fevd): %s\n",
              prior_name, paste(dim(fevd_arr), collapse = " x ")))
  
  # Reduce over draws (dim 1) -> [response, horizon, impulse]
  fevd_med <- apply(fevd_arr, c(2, 3, 4), median, na.rm = TRUE)
  n_vars   <- dim(fevd_med)[1]
  n_h_fevd <- dim(fevd_med)[2]
  
  fevd_tbl <- bind_rows(lapply(seq_len(n_vars), function(resp_i) {
    mat <- fevd_med[resp_i, , ]   # [horizon, impulse]
    as_tibble(mat) %>%
      setNames(TARGET_VARS[seq_len(n_vars)]) %>%
      mutate(h        = seq_len(n_h_fevd),
             response = TARGET_VARS[resp_i]) %>%
      pivot_longer(-c(h, response),
                   names_to  = "impulse",
                   values_to = "share")
  }))
  
  cat(sprintf("[%s] FEVD at h=1,4,8:\n", prior_name))
  fevd_tbl %>%
    filter(h %in% c(1, 4, 8)) %>%
    mutate(share = round(share * 100, 1)) %>%
    pivot_wider(names_from = impulse, values_from = share) %>%
    arrange(response, h) %>%
    print(n = Inf, width = Inf)
  cat("\n")
  
  col_fill <- c(gdp_g    = COL_GDP,
                cpi_g    = COL_CPI,
                bond_dif = COL_BOND)
  lbl      <- c(gdp_g    = "GDP growth",
                cpi_g    = "CPI inflation",
                bond_dif = "Bond yield \u0394")
  
  fig_fevd <- fevd_tbl %>%
    filter(response %in% TARGET_VARS) %>%
    mutate(
      response = factor(response, levels = TARGET_VARS,
                        labels = c("GDP growth", "CPI inflation",
                                   "Bond yield \u0394")),
      impulse  = lbl[impulse]
    ) %>%
    ggplot(aes(x = h, y = share * 100, fill = impulse)) +
    geom_area(alpha = 0.85, colour = "white", linewidth = 0.2) +
    facet_wrap(~ response, ncol = 1) +
    scale_fill_manual(values = setNames(unname(col_fill), unname(lbl)),
                      name = "Shock from:") +
    scale_x_continuous(breaks = seq(1, IRF_HORIZON, by = 2)) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title    = sprintf("FEVD — %s | %s", res$id, PRIOR_LABELS[prior_name]),
      subtitle = "Posterior median | Cholesky: GDP -> CPI -> Bond",
      x        = "Horizon (quarters)",
      y        = "% of forecast error variance"
    ) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  print(fig_fevd)
  # ggsave(file.path(here("output","figures"),
  #                  sprintf("bvar_03_fevd_%s.png", prior_name)),
  #        fig_fevd, width = 9, height = 10, dpi = 150)
}

# =============================================================================
# SECTION 3: MCMC CONVERGENCE
# =============================================================================
cat("\n=== SECTION 3: MCMC CONVERGENCE ===\n\n")
cat(sprintf("%-32s %-14s %-12s %-10s\n",
            "Model", "Acc.rate", "Lambda_mean", "Lambda_sd"))
cat(strrep("-", 70), "\n")

for (prior_name in names(PRIORS)) {
  best_row <- best_per_prior %>% filter(prior == prior_name)
  res      <- bvar_results[[prior_name]][[best_row$sample]][[
    paste0("p", best_row$p)]]
  if (is.null(res)) next
  mod <- res$model
  
  # Acceptance rate: stored as count in mod$meta$accepted
  acc_rate <- tryCatch(
    mod$meta$accepted / mod$meta$n_save,
    error = function(e) NA_real_
  )
  acc_flag <- if (!is.na(acc_rate) && acc_rate >= 0.25 && acc_rate <= 0.45) "v" else "!"
  
  lambda_chain <- tryCatch({
    h <- mod$hyper
    if ("lambda" %in% colnames(h)) h[, "lambda"] else h[, 1]
  }, error = function(e) NULL)
  
  cat(sprintf("%-32s %5.1f%% [%s]  %-12.3f %-10.3f\n",
              res$id,
              acc_rate * 100,
              acc_flag,
              if (!is.null(lambda_chain)) mean(lambda_chain) else NA_real_,
              if (!is.null(lambda_chain)) sd(lambda_chain)   else NA_real_))
  
  if (prior_name == "dummyobs") {
    for (hp in c("soc", "sur")) {
      chain <- tryCatch(mod$hyper[, hp], error = function(e) NULL)
      if (!is.null(chain)) {
        cat(sprintf("  |-- %-28s %-14s %-12.3f %-10.3f\n",
                    hp, "", mean(chain), sd(chain)))
      }
    }
  }
  
  # Trace plots
  tryCatch({
    h       <- mod$hyper
    n_hyper <- min(ncol(h), 3)
    par(mfrow = c(n_hyper, 1), mar = c(3, 4, 2, 1))
    for (j in seq_len(n_hyper)) {
      plot(h[, j], type = "l", col = COL_GDP,
           main = sprintf("Trace: %s (%s)", colnames(h)[j], res$id),
           xlab = "Draw", ylab = colnames(h)[j], cex.main = 0.85)
      abline(h = mean(h[, j]), col = COL_CPI, lty = 2, lwd = 1.5)
    }
    par(mfrow = c(1, 1))
  }, error = function(e) cat("  Trace plot error:", e$message, "\n"))
  cat("\n")
}

# =============================================================================
# SECTION 4: POSTERIOR PREDICTIVE RESIDUALS
# =============================================================================
cat("=== SECTION 4: POSTERIOR PREDICTIVE RESIDUALS ===\n")
cat("Model: bvar_minnesota_post08_p1 (comparable to frequentist post08_p3)\n\n")

res_post08 <- bvar_results[["minnesota"]][["post08"]][["p1"]]

if (!is.null(res_post08)) {
  mod       <- res_post08$model
  beta_mean <- apply(mod$beta, c(2, 3), mean)  # [K*p+1, n_vars]
  
  dat <- var_adj %>%
    filter(date >= samples[["post08"]]$start_date) %>%
    select(date,
           gdp_g    = gdp_g_adj,
           cpi_g    = cpi_g_adj,
           bond_dif = bond_dif_adj) %>%
    arrange(date)
  
  p      <- res_post08$p
  n      <- nrow(dat)
  Y      <- as.matrix(dat[, TARGET_VARS])
  X      <- do.call(rbind, lapply((p + 1):n, function(t) {
    c(1, as.vector(t(Y[(t - 1):(t - p), ])))
  }))
  Y_     <- Y[(p + 1):n, ]
  resids <- Y_ - X %*% beta_mean
  
  cat(sprintf("%-12s %-18s %-18s %-8s\n",
              "Variable", "JB p-val", "LB(10) p-val", "Normal?"))
  cat(strrep("-", 58), "\n")
  
  panels_acf <- list()
  panels_qq  <- list()
  
  for (vi in seq_along(TARGET_VARS)) {
    vname <- TARGET_VARS[vi]
    r     <- resids[, vi]
    col   <- c(COL_GDP, COL_CPI, COL_BOND)[vi]
    
    jb   <- tryCatch(tseries::jarque.bera.test(r), error = function(e) NULL)
    lb   <- Box.test(r, lag = 10, type = "Ljung-Box")
    jb_p <- if (!is.null(jb)) jb$p.value else NA_real_
    
    cat(sprintf("%-12s %-18s %-18s %-8s\n",
                vname,
                if (!is.na(jb_p))
                  sprintf("%.4f [%s]", jb_p, if (jb_p > 0.05) "ok" else "!!")
                else "NA",
                sprintf("%.4f [%s]", lb$p.value,
                        if (lb$p.value > 0.05) "ok" else "!!"),
                if (!is.na(jb_p) && jb_p > 0.05) "Yes" else "No"))
    
    # ACF panel
    acf_v  <- acf(r, lag.max = 12, plot = FALSE)
    ci     <- qnorm(0.975) / sqrt(length(r))
    acf_df <- tibble(lag = acf_v$lag[-1], acf = acf_v$acf[-1], ci = ci)
    
    panels_acf[[vname]] <- ggplot(acf_df, aes(x = lag, y = acf)) +
      geom_hline(yintercept = 0, color = COL_GREY) +
      geom_hline(yintercept =  ci, linetype = "dashed",
                 color = COL_CPI, linewidth = 0.4) +
      geom_hline(yintercept = -ci, linetype = "dashed",
                 color = COL_CPI, linewidth = 0.4) +
      geom_segment(aes(xend = lag, yend = 0),
                   colour = col, linewidth = 0.8) +
      geom_point(colour = col, size = 2) +
      scale_x_continuous(breaks = seq(1, 12, by = 2)) +
      ylim(-0.6, 0.6) +
      labs(title = sprintf("ACF: %s", vname), x = "Lag", y = "ACF")
    
    # QQ panel
    qq_df <- tibble(th = qnorm(ppoints(length(r))), s = sort(r))
    panels_qq[[vname]] <- ggplot(qq_df, aes(x = th, y = s)) +
      geom_abline(slope = sd(r), intercept = mean(r),
                  color = COL_CPI, linewidth = 0.7) +
      geom_point(colour = col, alpha = 0.6, size = 1.5) +
      labs(title = sprintf("QQ: %s", vname),
           x     = "Theoretical quantiles",
           y     = "Sample quantiles")
  }
  
  fig_resid <- (panels_acf[[1]] | panels_acf[[2]] | panels_acf[[3]]) /
    (panels_qq[[1]]  | panels_qq[[2]]  | panels_qq[[3]]) +
    plot_annotation(
      title    = "BVAR residuals — bvar_minnesota_post08_p1",
      subtitle = "Post-GFC | p=1 | Top: ACF | Bottom: QQ",
      theme    = theme(
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "gray40")
      )
    )
  print(fig_resid)
  
} else {
  cat("bvar_minnesota_post08_p1 not found — skipping.\n")
}

# =============================================================================
# SECTION 5: PRIOR VS POSTERIOR
# =============================================================================
cat("\n=== SECTION 5: PRIOR VS POSTERIOR ===\n\n")

cat("Lambda (prior mode = 0.2):\n")
bvar_summary %>%
  group_by(prior, sample) %>%
  summarise(mean_lambda = round(mean(lambda, na.rm = TRUE), 3),
            .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = mean_lambda) %>%
  print(width = Inf)

cat("\nDummy obs — soc / sur (prior mode = 1.0):\n")
bvar_summary %>%
  filter(prior == "dummyobs", !is.na(soc)) %>%
  select(id, sample, p, soc, sur) %>%
  print(n = Inf, width = Inf)

# Lambda dot plot
prior_cols <- c(minnesota     = COL_GDP,
                normalwishart = COL_CPI,
                dummyobs      = COL_BOND)

fig_pp <- ggplot(bvar_summary,
                 aes(x = lambda, y = prior, colour = prior, shape = sample)) +
  geom_vline(xintercept = 0.2, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_point(size = 3, alpha = 0.8,
             position = position_jitter(height = 0.15, seed = 42)) +
  annotate("text", x = 0.22, y = 0.6, label = "Prior mode (0.2)",
           size = 3, color = "black", hjust = 0) +
  scale_colour_manual(values = prior_cols,
                      labels = PRIOR_LABELS, name = "Prior") +
  scale_y_discrete(labels = PRIOR_LABELS) +
  scale_shape_manual(values = c(full = 16, post08 = 17, post15 = 15),
                     name = "Sample") +
  xlim(0, 1) +
  labs(
    title    = "Lambda: prior mode vs posterior means",
    subtitle = "All 27 models | vertical line = prior mode (0.2)",
    x        = "Posterior mean of lambda",
    y        = NULL
  )
print(fig_pp)

# SOC / SUR dot plot
fig_soc_sur <- bvar_summary %>%
  filter(prior == "dummyobs") %>%
  select(sample, p, soc, sur) %>%
  pivot_longer(c(soc, sur), names_to = "hp", values_to = "post_mean") %>%
  ggplot(aes(x = post_mean, y = hp, colour = sample, shape = factor(p))) +
  geom_vline(xintercept = 1.0, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_point(size = 3, alpha = 0.8,
             position = position_jitter(height = 0.1, seed = 42)) +
  annotate("text", x = 1.05, y = 1.5, label = "Prior mode (1.0)",
           size = 3, color = "black", hjust = 0) +
  scale_colour_manual(
    values = c(full = COL_GDP, post08 = COL_CPI, post15 = COL_BOND),
    name   = "Sample"
  ) +
  labs(
    title    = "Dummy obs: soc and sur posterior means",
    subtitle = "Dummy observations prior | vertical line = prior mode (1.0)",
    x        = "Posterior mean",
    y        = "Hyperparameter",
    shape    = "Lag order"
  )
print(fig_soc_sur)

cat("\nDiagnostics complete.\n")
cat("Next: bvar/04_bvar_extended.R\n")