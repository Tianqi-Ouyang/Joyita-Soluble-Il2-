# Soluble IL-2R × ICI

Quarto pipeline counting soluble IL-2 receptor (sIL2R, CD25) draws at Mass General
Brigham relative to immune checkpoint inhibitor (ICI) start.

Site: <https://tianqi-ouyang.github.io/Joyita-Soluble-Il2-/>

## Layout

| Path | Role |
|---|---|
| `R/functions.R` | Definitions: result classification, ICI start date, cohort classification, small-cell helpers |
| `R/build.R` | **Local-only** PHI step: reads RPDR + ICI files, writes `cache/` and `output/` (both gitignored) |
| `R/paths_template.R` | Copy to `R/paths_local.R` (gitignored) and set the local data paths |
| `index.qmd`, `data_management.qmd`, `analysis.qmd` | Site pages; read only the aggregate `cache/summary.rds` |

## Rebuild

```bash
cp R/paths_template.R R/paths_local.R   # edit paths
Rscript R/build.R
quarto render                           # or: quarto publish gh-pages
```

## Data privacy

No patient data is in this repository. The RPDR extract, the ICI data pull, the
EMPI-keyed cache and the chart-review workbook stay local. Rendered pages show
aggregate counts only, with counts from 1 to 10 shown as `<11`.
