# Macro Forecasting

## Project Overview
This project analyses Swiss export dynamics using Vector Autoregression (VAR) 
to forecast trade flows across key economic partners including Germany, China, 
Europe, and the Americas. Data spans January 2005 to January 2026.

## Data
The `data/` folder is not tracked by Git (too large).
Download the following files and place them in a local `data/` folder:
- `ons_data_hierarchy_FS26.xlsx` — Swiss export time series (monthly)
- `swiss_nowcast_data.json` — Swiss nowcast metadata

## Repository Structure
```
macro_forecasting/
├── scripts/
│   ├── exploratory.R       # initial data exploration
│   └── var_analysis.R      # main VAR model
├── .gitignore
├── README.md
└── macro_forecasting.Rproj
```

## Requirements
- R version 4.x+
- Packages: vars, readxl, jsonlite, here, tidyverse, tseries, lubridate

## Authors
Ronja Hegelbach & Naéla Gruber
University of Zurich — Macro Forecasting, Spring 2026
