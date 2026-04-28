# Macroeconomic Forecasting — Swiss GDP, CPI & Bond Yield
**Ronja Hegelbach & Naéla Gruber**  
Macroeconomic Forecasting Seminar, Spring 2026 | KOF ETH Zurich / University of Zurich  
Supervisor: Alexander Rathke

---

## Overview

Joint VAR forecasting of three Swiss macroeconomic variables:
- **GDP growth** — QoQ %, KOF sports-event corrected
- **CPI inflation** — QoQ %, headline
- **10Y bond yield change** — QoQ pp (level is I(1); first-differenced)

Baseline models vary lag order (p=1–5) and sample period (full 2000–2025, post-GFC 2009–2025, post-2015 2015–2025). Evaluated pseudo out-of-sample against AR(1) benchmark via expanding window (2010 Q1 – 2025 Q4, horizons h=1,2,4,8). COVID quarters (2020 Q1 – 2021 Q4) handled via exogenous dummies in the VAR and excluded from RMSE comparisons.

---

## Data

`data/` is not tracked by Git. Place the following files there locally:

- `swiss_nowcast_data.json` — KOF nowcast dataset (~4340 series)
- `swiss_nowcast_metadata.xlsx` — variable descriptions (366 base series)

**VAR dataset:** Q2 2000 – Q4 2025, 103 quarterly observations.

---

## Repository Structure

```
macro_forecasting/
├── scripts/
│   ├── 00 exploration/
│   │   ├── 00_setup.R                       # packages, data loading, variable extraction
│   │   └── 01_exploration_target_variables.R # plots, ADF, outliers, Granger tests
│   └── 01 var models/
│       ├── 01_var_baseline.R                # VAR(1-5) × 3 samples, forecasts
│       ├── 02_var_baseline_evaluation.R     # AR(1) benchmark, RMSE, DM tests
│       ├── 03_var_baseline_diagnostics.R    # IRF, FEVD, rolling stability, residuals
│       ├── 04_var_extended_variable_selection.R    #
│       ├── 05_var_extended_estimation.R    #
│       ├── 06_var_extended_evaluation.R    #
│       ├── 07_var_extended_diagnostics.R   # 
│   └── 02 bvar models/
│       ├── 01_bvar_baseline.R               #
│       
├── output/
│   ├── figures/
│   └── tables/
└── README.md
```

---

## How to Run

Always source `00_setup.R` first — it loads all packages and data. Run scripts in numbered order within each folder.

```r
source(here("scripts", "exploration", "00_setup.R"))
```

> `dplyr` is loaded last in `00_setup.R` to prevent masking by `vars`/`MASS`. All `ggsave()` and `write.csv()` calls are commented out — uncomment to save outputs.

---

## Key Results

**Q1 2026 GDP nowcast:** +0.39% QoQ (Nowcasting Lab benchmark: +0.30%)

**Out-of-sample RMSE (ex-COVID):**
- GDP: AR(1) beats all VAR models — classic result, simple benchmarks hard to beat
- CPI: VAR beats AR(1) at h=4,8 (~70% of AR(1) RMSE) via bond yield channel
- Bond: VAR better only at long horizons (h=8, ~45% of AR(1) RMSE)

**Preferred model:** `post08_p3` (Post-GFC sample, p=3) — lowest GDP residual SD, only model with normal residuals in all three equations, best out-of-sample GDP at h=4.

**Granger causality:** Bond yield → GDP (p=0.04) and Bond yield → CPI (p=0.004). No reverse causality.

---

## Planned Extensions

- Extended VAR with CHF/EUR and oil price as exogenous predictors
- Bayesian VAR with Minnesota prior
- Full model comparison table and BoE-style fan chart

---

## Authors

**Ronja Hegelbach & Naéla Gruber** — University of Zurich  
Supervisor: Alexander Rathke, KOF ETH Zurich
