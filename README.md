# Macro Forecasting — Swiss GDP Nowcasting

**Ronja Hegelbach & Naéla Gruber**  
Macroeconomic Forecasting Seminar, Spring 2026  
KOF ETH Zurich / University of Zurich

---

## Project Overview

This project nowcasts Swiss GDP growth using Vector Autoregression (VAR) models estimated on the KOF nowcast dataset. We target **Swiss GDP excluding sporting events** (FIFA/UEFA revenues create large one-off spikes in headline GDP that cannot be predicted by economic indicators — KOF themselves use the sports-excluded series in their own nowcasting model).

Our baseline is a **bridge VAR**: monthly indicators are aggregated to quarterly frequency and a standard VAR is estimated. We evaluate pseudo out-of-sample performance against a naive AR(1) benchmark and compare our Q1 2026 nowcast to the KOF benchmark (+0.30% QoQ).

Planned extensions include a Bayesian VAR and a mixed-frequency VAR.

---

## Data

The `data/` folder is not tracked by Git (too large). Download the following files and place them in a local `data/` folder:

- `swiss_nowcast_data.json` — KOF nowcast predictor series (~4340 series, each with 14 transformation variants)
- `swiss_nowcast_metadata.xlsx` — variable descriptions, sources, and frequencies (366 rows)

### Data structure

The JSON stores each series as a numeric vector under a key like `swobs085q_lvl`. Dates are stored separately under `nowcast_raw$dates[[key]]` in `DD.MM.YYYY` format — there is no shared date index. Each base series has 14 variants: `_lvl`, `_lvl_detrended`, `_pct_1m/3m/1y`, `_dif_1m/3m/1y`, and detrended versions of each.

### Target variable

`ch_kof_modelinput_gdpos_pct_3m` — Swiss real GDP excluding sporting events, quarter-on-quarter growth rate (%). Source: SECO via KOF. Last observed: Q4 2025 (+0.155%). Nowcast target: Q1 2026.

### Selected predictors

| Key | Description | Transformation |
|-----|-------------|----------------|
| `swobs085q` | CH OECD composite leading indicator | Level |
| `swcnfbusq` | KOF business situation survey | Level |
| `bdiptot_g` | German industrial production | MoM growth rate |
| `ekeusesig` | EA economic sentiment indicator | Level |
| `swxsfec_` | CHF/EUR exchange rate | MoM growth rate |

Estimation sample: 2004 Q1 – 2025 Q4 (constrained by `swcnfbusq` which starts January 2004). COVID dummy added for 2020 Q2 and Q3.

---

## Repository Structure

```
macro_forecasting/
│
├── R/                           # shared utility functions
│   ├── data_utils.R             # pull_series, to_quarterly,
│   │                            # get_partial_avg, build_quarterly_panel
│   └── model_utils.R           # bridge_nowcast, compute_metrics,
│                                # fit_var, ar1_forecast
│
├── scripts/
│   ├── exploration/             # initial data inspection
│   │   ├── 01_inspect.R         # raw data structure
│   │   ├── 02_inspect_json.R    # JSON naming convention
│   │   ├── 02_explore_nowcast.R # plots, stationarity, coverage
│   │   └── 03_inspect_gdp.R    # GDP target variable
│   │
│   ├── analysis/                # main modelling pipeline
│   │   ├── 04_variable_selection.R  # predictor selection
│   │   ├── 05_bridge_var.R          # model estimation + Q1 2026 nowcast
│   │   └── 06_evaluation.R          # pseudo out-of-sample evaluation
│   │
│   └── extensions/              # planned extensions
│       ├── 07_bvar.R            # Bayesian VAR (Minnesota prior)
│       ├── 08_mfvar.R           # mixed-frequency VAR
│       └── 09_comparison.R      # model comparison + final table
│
├── output/
│   ├── figures/                 # saved plots
│   └── tables/                  # saved result tables
│
├── .gitignore
├── README.md
└── macro_forecasting.Rproj
```

---

## How to Run

1. Clone the repository and place data files in `data/`
2. Open `macro_forecasting.Rproj` in RStudio
3. Install required packages (or restore with `renv::restore()` if lockfile present):

```r
install.packages(c("readxl", "jsonlite", "here", "dplyr", "tidyr",
                   "lubridate", "zoo", "ggplot2", "tseries", "vars"))
```

4. Run scripts in order within each folder — exploration first, then analysis

**Note:** Load `R/data_utils.R` and `R/model_utils.R` at the top of each analysis script:

```r
source(here("R", "data_utils.R"))
source(here("R", "model_utils.R"))
```

---

## Key Results (baseline bridge VAR)

| Sample | Relative RMSE (VAR / AR1) | VAR better? |
|--------|--------------------------|-------------|
| Full 2015–2025 | 0.889 | Yes |
| Excluding COVID | 2.315 | No |
| COVID quarters only | 0.526 | Yes |

Q1 2026 nowcast: **+1.07% QoQ** (KOF benchmark: +0.30%). Wide 95% interval [−0.74%, +2.88%] — KOF estimate lies within our uncertainty range. Gap likely reflects KOF's access to more timely and granular data.

---

## Authors

Ronja Hegelbach & Naéla Gruber — University of Zurich  
Macroeconomic Forecasting Seminar, Spring 2026  
Supervisor: Alexander Rathke, KOF ETH Zurich
