# ============================================================
# analysis/04_var_diagnostics.R
# VAR model diagnostics — IRFs and FEVD
# ============================================================
#
# PURPOSE:
#   Interpret the estimated bridge VAR economically. Having
#   estimated the model (script 02) and evaluated its forecast
#   accuracy (script 03), we now ask: what does the model
#   tell us about the Swiss economy? Two tools are used:
#
#   1. IMPULSE RESPONSE FUNCTIONS (IRFs)
#      How does Swiss GDP respond over time to a one-time
#      shock in each predictor? Traces the dynamic
#      transmission channels — e.g. how a drop in German IP
#      propagates to Swiss GDP over subsequent quarters.
#      This is the core of the "story" behind the forecast.
#
#   2. FORECAST ERROR VARIANCE DECOMPOSITION (FEVD)
#      Of all the uncertainty in the GDP forecast h quarters
#      ahead, what fraction is attributable to shocks in
#      each variable? Shows which predictors drive GDP
#      uncertainty most at short vs long horizons.
#
# NOTE ON PLOTTING:
#   All plots use ggplot2. The vars package plot() method is
#   avoided because it uses base R graphics and triggers
#   interactive "press Enter" prompts in RStudio. Instead we
#   extract the IRF data manually and plot with ggplot2 so
#   all plots render immediately without user input.
#
# REQUIRES:
#   Self-contained — re-estimates model from scratch.
#
# KEY FINDINGS:
#   IRFs (GDP response to a one-SD shock in each predictor):
#     All predictors show the same oscillating pattern —
#     positive at h=1, reversal at h=2, dissipating by h=4.
#     Only the CLI response at h=1 (~+0.45pp) is close to
#     statistically significant (band nearly excludes zero).
#     All other responses have wide confidence bands that
#     include zero throughout, reflecting the fundamental
#     difficulty of predicting quarterly GDP with a short
#     sample. CHF/EUR and KOF show the weakest responses.
#     The CLI reversal at h=2 reflects its leading indicator
#     nature — the signal is front-loaded into h=1.
#
#   FEVD (share of GDP forecast variance, %):
#     h=1:  GDP own = 100% (Cholesky: by construction)
#     h=2:  GDP = 69.8  CLI = 15.2  EA ESI = 4.5
#           DE IP = 4.4  CHF/EUR = 3.2  KOF = 2.8
#     h=4:  GDP = 67.7  CLI = 14.2  EA ESI = 7.1
#           DE IP = 4.3  CHF/EUR = 3.8  KOF = 2.9
#     h=8:  GDP = 67.3  CLI = 14.6  EA ESI = 7.0
#           DE IP = 4.3  CHF/EUR = 3.8  KOF = 2.9
#
#   The CLI is the dominant external driver of GDP uncertainty
#   at all horizons (~15%), consistent with its design as a
#   leading indicator of the Swiss business cycle. EA Sentiment
#   and German IP contribute roughly equally (~4-7%) reflecting
#   Switzerland's export dependence on Europe. The picture
#   stabilises quickly after h=2 — medium-term uncertainty
#   is overwhelmingly driven by GDP's own dynamics (67%).
#
# ============================================================

library(here)
library(jsonlite)
library(dplyr)
library(tidyr)
library(vars)
library(ggplot2)

select <- dplyr::select

source(here("R", "data_utils.R"))
source(here("R", "model_utils.R"))

json_text   <- paste(readLines(here("data", "swiss_nowcast_data.json"),
                               encoding = "UTF-8"), collapse = "")
nowcast_raw <- fromJSON(json_text, flatten = TRUE)


# ============================================================
# 1. REBUILD MODEL
# ============================================================

gdp <- pull_series("ch_kof_modelinput_gdpos_pct_3m") %>%
  rename(gdp = value)

predictors <- list(
  cli    = pull_series("swobs085q_lvl"),
  kof    = pull_series("swcnfbusq_lvl"),
  de_ip  = pull_series("bdiptot_g_pct_1m"),
  ea_esi = pull_series("ekeusesig_lvl"),
  chfeur = pull_series("swxsfec__pct_1m")
)

pred_names <- names(predictors)
quarterly  <- build_quarterly_panel(gdp, predictors)

est_data <- quarterly %>%
  filter(qdate >= as.Date("2004-01-01"),
         qdate <= as.Date("2025-10-01"))

result  <- fit_var(est_data, pred_names, lag_max = 4)
var_mod <- result$model
p_aic   <- result$p

cat("VAR(", p_aic, ") estimated on", nrow(est_data), "quarters\n")


# ============================================================
# 2. IMPULSE RESPONSE FUNCTIONS — ggplot2
# ============================================================

# Orthogonalised IRFs via Cholesky decomposition.
# Variable ordering:  gdp → cli → kof → de_ip → ea_esi → chfeur
# GDP is ordered first — it does not respond to predictors
# within the same quarter (reasonable given publication lags).
# Swiss variables precede foreign ones.
#
# We extract the IRF data from the vars object and reshape
# to long format for ggplot2 — this avoids base R plot()
# and its interactive prompting entirely.

cat("\n--- Impulse Response Functions ---\n")
cat("Horizon: 8 quarters | 95% bootstrap CI | 200 runs\n\n")

# Pretty labels for plot panels
pred_labels <- c(
  cli    = "CH OECD CLI",
  kof    = "KOF Business Situation",
  de_ip  = "DE Ind. Production (MoM%)",
  ea_esi = "EA Economic Sentiment",
  chfeur = "CHF/EUR (MoM%)"
)

# Compute all IRFs at once
irf_obj <- irf(var_mod,
               impulse  = pred_names,
               response = "gdp",
               n.ahead  = 8,
               boot     = TRUE,
               ci       = 0.95,
               runs     = 200)

# Extract into a tidy data frame
irf_df <- bind_rows(lapply(pred_names, function(nm) {
  data.frame(
    predictor = nm,
    label     = pred_labels[nm],
    horizon   = 0:8,
    irf       = as.numeric(irf_obj$irf[[nm]]),
    lower     = as.numeric(irf_obj$Lower[[nm]]),
    upper     = as.numeric(irf_obj$Upper[[nm]])
  )
}))

# Plot all 5 IRFs in one faceted figure — no prompting
ggplot(irf_df, aes(x = horizon)) +
  geom_hline(yintercept = 0, color = "red",
             linewidth = 0.4, linetype = "solid") +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = irf), color = "black", linewidth = 0.7) +
  facet_wrap(~ label, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 0:8) +
  labs(
    title    = "Impulse Response Functions — GDP response",
    subtitle = paste("Response of GDP to a one-SD shock in each predictor",
                     "\nShaded band = 95% bootstrap CI (200 runs)"),
    x        = "Horizon (quarters)",
    y        = "Response (pp)"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.text  = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank())


# ============================================================
# 3. FORECAST ERROR VARIANCE DECOMPOSITION — ggplot2
# ============================================================

cat("\n--- Forecast Error Variance Decomposition ---\n")
cat("Share of GDP forecast variance attributable to each variable\n\n")

fevd_result <- fevd(var_mod, n.ahead = 8)
gdp_fevd    <- fevd_result$gdp

# Print table
cat("Horizon  gdp    cli    kof    de_ip  ea_esi chfeur\n")
cat(strrep("-", 60), "\n")
for (h in c(1, 2, 4, 8)) {
  row <- round(gdp_fevd[h, ] * 100, 1)
  cat(sprintf("h = %-4d  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n",
              h, row["gdp"], row["cli"], row["kof"],
              row["de_ip"], row["ea_esi"], row["chfeur"]))
}

# Plot stacked area chart
fevd_df <- as.data.frame(gdp_fevd) %>%
  mutate(horizon = 1:nrow(gdp_fevd)) %>%
  pivot_longer(-horizon,
               names_to  = "variable",
               values_to = "share") %>%
  mutate(variable = factor(variable,
                           levels = c("gdp", "cli", "kof",
                                      "de_ip", "ea_esi", "chfeur"),
                           labels = c("GDP (own)", "CH OECD CLI",
                                      "KOF Survey", "DE Ind. Prod.",
                                      "EA Sentiment", "CHF/EUR")))

ggplot(fevd_df, aes(x = horizon, y = share * 100,
                    fill = variable)) +
  geom_area(alpha = 0.85) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  scale_x_continuous(breaks = 1:8) +
  labs(title    = "Forecast Error Variance Decomposition — GDP",
       subtitle = "Share of GDP forecast variance explained by each variable (%)",
       x        = "Horizon (quarters)",
       y        = "Share (%)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")