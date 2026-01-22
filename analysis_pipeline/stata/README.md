# Stata Pipeline for AlphaFold2 Impact Analysis

Cleaned-up Stata implementation based on Danilo Messinese's original code, with added diagnostics.

## Quick Start

1. Open Stata
2. Edit `00_master.do` and set `PROJECT_PATH` to your project root
3. Run: `do "analysis_pipeline/stata/00_master.do"`

## Files

| File | Purpose |
|------|---------|
| `00_master.do` | Master runner - sets globals, runs all steps |
| `01_cem_matching.do` | CEM matching and panel construction |
| `02_event_study.do` | Event study regressions (LEVELS, DY, DASINH) |
| `03_diagnostics.do` | **NEW** - Honest assessment and diagnostics |

## What's Different from Danilo's Code

### Improvements:
1. **Relative paths** - No hardcoded Dropbox paths
2. **Diagnostics section** - Shows raw trends, % changes, data truncation warnings
3. **Honest assessment** - Explains what the data actually shows
4. **Better organization** - Modular .do files, clear output structure

### Key Diagnostic Outputs:
- Raw means by treatment/period (not just regression coefficients)
- Percentage changes (shows curves have same shape)
- Data truncation check for 2023
- Simple DiD calculation showing why LEVELS vs DY differ

## Output Structure

```
analysis_pipeline/stata/
├── output/
│   ├── derived/
│   │   ├── matched_panel_semester.dta
│   │   ├── gene_level_features.dta
│   │   └── matched_genes.dta
│   ├── eventstudy_LEVELS_*.csv
│   ├── eventstudy_DY_*.csv
│   ├── eventstudy_DASINH_*.csv
│   └── diagnostics_means_by_period.csv
└── figures/
    └── eventstudy_*.png
```

## Requirements

- Stata 15+
- Packages (auto-installed): `cem`, `reghdfe`, `ftools`, `estout`

## Key Finding from Diagnostics

The "positive effect" in the ΔY specification is driven by:
1. Both groups having **same-shaped curves** (just different scales)
2. Control declining more in **absolute** terms (not percentage)
3. ΔY comparing absolute changes, making control's bigger drop look like "treated doing better"

This is measuring **deceleration difference**, not growth.
