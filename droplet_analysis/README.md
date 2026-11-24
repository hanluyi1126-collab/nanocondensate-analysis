# Droplet models and Supplementary Figures S1–S2

Figure S1 reports diffusion, capillary energy and arrival at contact for equal spherical droplets. Figure S2 shows a phenomenological state map and the molecular fractions required for two interfacial budgets. The inputs and equations are given in Supplementary Data 1, pp. S2–S7.

## Reproduce the figures and check the tables

Run from this directory, or pass the path to the script from another directory. Outputs are placed under this module's `outputs/`. Figure export uses ggplot2, patchwork, svglite and ragg. The archived versions are R 4.4.1, ggplot2 4.0.2, patchwork 1.3.2, svglite 2.2.2 and ragg 1.5.1. The numerical check uses base R.

```sh
Rscript --vanilla verify_droplet.R
Rscript --vanilla plot_s1_s2.R
```

`verify_droplet.R` reconstructs the 180 radii at all three surface potentials and compares Eqs. S1–S8 with `data/drop/fusion_verified.csv`. It evaluates the Fuchs integral with both 24,000 midpoint nodes and adaptive quadrature (absolute tolerance 10⁻¹², relative tolerance 10⁻¹⁰). Subtracting the contact energy before taking the exponential keeps the largest barriers within floating-point range. The maximum midpoint–adaptive difference is 0.00023817 in log₁₀ W. The script also checks the four composition scores, ordered classes and three threshold-sensitivity counts against the archived state grid and Supplementary Table S3.

`plot_s1_s2.R` draws the archived tables and writes PDF, SVG, PNG and TIFF files under `outputs/figures/`. It shares the S1–S2 plotting code with the S1–S5 entry point. Figure S1 uses the full adaptive log rates. Figure S2B evaluates Eqs. S11–S12 for diameters of 50–2,000 nm. At a diameter of 100 nm, the charge and coverage fractions are 0.072 and 1.2, respectively.

## Recalculate the model tables

```sh
Rscript --vanilla legacy_fusion_figure.R
Rscript --vanilla legacy_state_model.R
```

The first script writes the original midpoint calculation and a diagnostic figure under `outputs/fusion/`. Its historical CSV field `log10_fusion_s` denotes log₁₀ of the calculated per-droplet contact rate in s⁻¹. Fusion would additionally require a post-contact success probability. Only the diagnostic plot uses a display floor at −50. The electrostatic calculation prescribes a constant surface potential in idealized water. Its κR range is 0.052–2.082, outside the thin-double-layer regime assumed by the proximity approximation. The 25 mV curve also approaches the thermal-voltage scale.

The state model evaluates Eq. S10 at 150 mM using prescribed coefficients and thresholds. It writes the archived grid, candidate-window ranges and a diagnostic plot under `outputs/state/`. Classes are assigned in the order given in Table S3. They represent conditions on dimensionless scores.

### Composition grid boundary

Figure S2A and Table S3 use 12,106 points. The original R filter compares the floating-point sum of two fractions directly with 0.85. At a spacing of 0.005, 50 pairs on the boundary round above 0.85. The inclusive limit specified in the supplementary methods gives 12,156 pairs when evaluated on integer grid indices. All 50 additional pairs fall in the network-arrest class at the nominal thresholds.

| Grid | Dispersed | Coarsening | Candidate window | Network arrest |
| --- | ---: | ---: | ---: | ---: |
| Archived figure and Table S3 | 5,412 | 2,459 | 2,318 | 1,917 |
| Inclusive composition boundary | 5,412 | 2,459 | 2,318 | 1,967 |

The state script exports `additional_boundary_points.csv` and `grid_class_counts.csv` alongside the archived-grid outputs. The candidate-window fraction is 2,318/12,106 (19.15%) on the archived grid and 2,318/12,156 (19.07%) on the inclusive grid.
