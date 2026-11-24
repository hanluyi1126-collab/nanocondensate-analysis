# Supplementary computational data

[Data S1](Data_S1.zip) contains the droplet calculations for Figs. S1 and S2. [Data S2](Data_S2.zip) contains finite-element outputs, trajectory means, bootstrap fits, kinetic tables and post-pore time courses for Figs. S3–S5 and Table S8. Each archive has its own README and CSV catalogue with paths relative to that archive.

To compare the plotting tables with the full archives and regenerate the 20,000 radius-stratified bootstrap fits for Fig. S4, run the following command from the repository root. It uses Python 3 and NumPy.

```sh
python3 supplementary-data/check_tables.py
```

The comparison checks the row count, row order and values in every plotted column. The four finite-element plotting tables use blank padding for fields that are read only in other panels. The full data for those fields are in Data S2.

## Archive SHA-256

```text
aab8697e84a81e8aa31eee20f5361062a895b40d51ba4b783a9d380386686d32  Data_S1.zip
7424c950e86204377e029edcbd2fc811890af8b020c9e8e71f93659f26fa201e  Data_S2.zip
```

The documentation and catalogues inside these archives were revised after commit `450ea41`. All 39 enclosed numerical and metadata files are byte-identical to the corresponding members at that commit.
