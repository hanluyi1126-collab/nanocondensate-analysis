# Droplet calculations and supplementary figure reproduction

This repository contains the R droplet models, the figure inputs for Supplementary Figs. S1–S5, and the archived computational tables in Data S1 and S2.

Run these commands from the repository root.

```sh
Rscript --vanilla droplet_analysis/verify_droplet.R
Rscript --vanilla figure_reproduction/figures.R
python3 supplementary-data/check_tables.py
```

The first command recalculates the droplet integrals and checks the state classifications with base R. Figure export requires ggplot2, patchwork, svglite and ragg. The Python check uses NumPy to repeat the 20,000 radius-stratified bootstrap fits from the 24 trajectory means.

| Files | What can be reproduced |
| --- | --- |
| [Droplet models](droplet_analysis/) | Fusion-model contact rates and the illustrative state grid. The module documents the 50-point difference between the original floating-point grid and the inclusive composition boundary. |
| [Supplementary figures](figure_reproduction/) | Figs. S1–S5 in PDF, SVG, PNG and TIFF from the figure tables. |
| [Supplementary data](supplementary-data/) | Full data tables, per-archive catalogues and checksums. The membrane figures use archived simulation outputs. |

The original figure code cited in the manuscript is available at [commit 80cf584](https://github.com/hanluyi1126-collab/nanocondensate-analysis/tree/80cf5840d59df7349260d9f7fa11b6532ef740aa). [Commit 450ea41](https://github.com/hanluyi1126-collab/nanocondensate-analysis/tree/450ea41b394bc3adc7e580e84eff7b32ebb06c8e) added the two data archives. The archive documentation has since been corrected. The numerical inputs to the published figures are preserved.
