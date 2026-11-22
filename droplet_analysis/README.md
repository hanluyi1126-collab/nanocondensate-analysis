# Droplet models and Supplementary Figures S1–S2

This module reproduces the manuscript’s droplet fusion and condensate state-map figures. All quantities are model calculations. Figure S1 covers encounter rates, capillary free energy and electrostatic contact rates. Figure S2 shows a dimensionless state map and interfacial molecule budgets.

## Reproduce the two figures

Run these commands from this directory. The recorded environment is R 4.4.1 with ggplot2 4.0.2, patchwork 1.3.2, svglite 2.2.2 and ragg 1.5.1. The scripts use Helvetica for plot text.

```sh
Rscript --vanilla plot_s1_s2.R
```

The script reads the archived numerical integration results in `data/drop/fusion_verified.csv` and the illustrative model grid in `data/drop/state_grid.csv`. It writes `Fig_S1` and `Fig_S2` as PDF, SVG, PNG and TIFF files to `outputs/figures/`.

## Run the original model scripts

```sh
Rscript --vanilla legacy_fusion_figure.R
Rscript --vanilla legacy_state_model.R
```

The first script performs the original midpoint integration and writes its calculation table and diagnostic figure to `outputs/fusion/`. Supplementary Fig. S1 uses the independently evaluated integration results supplied in `fusion_verified.csv`. The second script generates the state grid, a summary table and a diagnostic plot in `outputs/state/`. The model’s classifications are illustrative. Use `plot_s1_s2.R` for the manuscript figures.
