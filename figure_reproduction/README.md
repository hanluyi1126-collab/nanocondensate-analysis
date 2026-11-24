# Supplementary Figures S1–S5

Run this command from the repository root with R, ggplot2, patchwork, svglite and ragg installed.

```sh
Rscript --vanilla figure_reproduction/figures.R
```

From this directory, `Rscript --vanilla figures.R .` also works. The script reads the ten CSVs under `data/`, creates `figures/` and writes each figure as PDF, SVG, 300-dpi PNG and 600-dpi TIFF. A second run replaces these outputs. The script writes the local R and package versions to `R_session_info.txt`. The archived rendering used Helvetica, R 4.4.1, ggplot2 4.0.2, patchwork 1.3.2, svglite 2.2.2 and ragg 1.5.1. See `archived_R_session.txt` for the archived platform and library versions.

## Inputs and calculations

The CSVs contain the plotting columns from Data S1 and Data S2. The [supplementary data directory](../supplementary-data/) contains the full archives and a Python cross-check that also repeats the S4 bootstrap from the trajectory means. Some membrane CSVs retain a common header with blank cells in columns unused for that series. The required columns are checked separately for each panel.

| Figure | Contents |
| --- | --- |
| S1 | Encounter rates, capillary energies, electrostatic contact energies and contact-arrival rates from `drop/fusion_verified.csv`. |
| S2 | Model classifications from `drop/state_grid.csv` and the charge and coverage budgets calculated as functions of diameter. Both S1 and S2 use the same plotting code as the standalone droplet script. |
| S3 | Centre and peak tension, the peak position relative to the contact radius, mesh refinement and the empirical distribution of force residuals. Each centre tension is compared with the 3,840-element result at the same diameter. Relative errors are plotted as percentages. The residual panel includes all 201 diameter, 35 refinement, 1,554 cohesive-grid and 2,048 ensemble results. Its dashed vertical line marks the reference value of 10⁻⁶ pN. |
| S4 | Unweighted linear fit of holding force divided by 2π against restrained pore radius using four trajectory means at each of six radii. The intercept estimates the line tension. The histogram reads the 20,000 archived bootstrap intercepts, and dashed lines mark their 2.5th and 97.5th percentiles (R quantile type 7). |
| S5 | Counts across 336 parameter combinations classified by the two diameter endpoints. The middle panel shows the first 50% crossing only for the small-pore / large-wrap class. Its dashed line marks the subset median. The last panel shows three selected probability curves at the indicated defect reference times. |

For the archived data the S4 intercept is 5.174 pN, with a bootstrap percentile interval of 3.047 to 7.181 pN. The S5 endpoint counts are 42 pore/pore, 7 wrap/pore, 128 pore/wrap and 159 wrap/wrap (small/large). The median first crossing for the 128 pore/wrap combinations is 82.862 nm. These summaries describe the saved model outputs. The figures use precomputed finite-element solutions, kinetic calculations and bootstrap estimates.
