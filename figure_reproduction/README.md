# Reproduce Supplementary Figures S1–S5

The archived R script draws the five supplementary figures from saved numerical results. The ten CSV tables in `data/` retain the columns read by the plotting code. The four finite-element tables share a header so that R can combine their residuals. Blank cells in those tables correspond to columns unused for that data series.

The archived environment used R 4.4.1 with ggplot2 4.0.2, patchwork 1.3.2, svglite 2.2.2 and ragg 1.5.1. Package and platform details are recorded in `archived_R_session.txt`. The script requests the Helvetica font.

From this directory, run the following commands.

```sh
mkdir -p figures
Rscript --vanilla figures.R .
```

The script writes one PDF, PNG, SVG and TIFF for each figure in `figures/` and records the current R environment in `R_session_info.txt`.
