"""Cross-check the plotting tables against the data archives and refit Fig. S4."""

import csv
import io
import json
from itertools import zip_longest
from pathlib import Path
from zipfile import ZipFile

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FEM_PADDING = {
    "diameter.csv": {"membrane_elements"},
    "convergence.csv": {"maximum_tension_position_fraction"},
    "cohesive_grid.csv": {
        "diameter_nm", "membrane_elements", "centre_tension_pN_per_nm",
        "maximum_tension_pN_per_nm", "maximum_tension_position_fraction",
    },
}
FEM_PADDING["ensemble.csv"] = FEM_PADDING["cohesive_grid.csv"]


def read_csv(archive, name):
    with archive.open(name) as member, io.TextIOWrapper(member, encoding="utf-8") as text:
        reader = csv.DictReader(text)
        return reader.fieldnames, list(reader)


def check_catalog(archive, prefix):
    _, catalog = read_csv(archive, "data_catalog.csv")
    files = {p for p in archive.namelist() if p.startswith(prefix) and p.endswith(".csv")}
    if len(catalog) != len(files) or {r["file"] for r in catalog} != files:
        raise ValueError(f"{archive.filename}: catalogue does not match archive")
    for entry in catalog:
        columns, rows = read_csv(archive, entry["file"])
        if columns != entry["columns"].split("; ") or len(rows) != int(entry["rows"]):
            raise ValueError(f"Incorrect catalogue entry for {entry['file']}")


def check_projection(path, archive, member):
    columns, full_rows = read_csv(archive, member)
    padding = FEM_PADDING.get(path.name, set())
    with path.open(encoding="utf-8", newline="") as handle:
        projected = csv.DictReader(handle)
        if not set(projected.fieldnames).issubset(columns):
            raise ValueError(f"Unrecognized columns in {path}")
        for line, (row, source) in enumerate(zip_longest(projected, full_rows), start=2):
            if row is None or source is None:
                raise ValueError(f"Row-count mismatch: {path}")
            for field, value in row.items():
                expected = "" if field in padding else source[field]
                if value != expected:
                    raise ValueError(f"{path}:{line}: mismatch in {field}")


def check_lipid_fit(archive):
    _, trajectories = read_csv(archive, "Data_S2/lipid_replicates.csv")
    _, bootstrap = read_csv(archive, "Data_S2/bootstrap_verified.csv")
    if len({r["run"] for r in trajectories}) != 24 or len(trajectories) != 24:
        raise ValueError("Expected 24 distinct archived trajectories")
    radii = np.array([float(r["radius_nm"]) for r in trajectories])
    force = np.array([float(r["force_over_2pi_pN"]) for r in trajectories])
    unique_radii, counts = np.unique(radii, return_counts=True)
    if len(unique_radii) != 6 or not np.all(counts == 4):
        raise ValueError("Expected four trajectories at each of six radii")

    # In each resample all radius groups have equal size, so least squares on
    # the six group means gives the same coefficients as fitting all 24 points.
    by_radius = np.array([force[radii == r] for r in unique_radii])
    rng = np.random.Generator(np.random.PCG64(20260803))
    draws = rng.integers(0, 4, size=(20_000, 6, 4))
    means = np.take_along_axis(by_radius[None, :, :], draws, axis=2).mean(axis=2)
    design = np.column_stack([np.ones(6), -unique_radii])
    fits = np.linalg.lstsq(design, means.T, rcond=None)[0].T
    archived = np.array([[float(r["line_tension_pN"]), float(r["residual_tension_pN_per_nm"])]
                         for r in bootstrap])
    np.testing.assert_allclose(fits, archived, atol=1e-10, rtol=0)
    summary = json.loads(archive.read("Data_S2/lipid_summary.json"))
    fit = np.linalg.lstsq(design, by_radius.mean(axis=1), rcond=None)[0]
    np.testing.assert_allclose(fit, [summary["line_tension_pN"],
                                    summary["residual_membrane_tension_pN_per_nm"]],
                               atol=1e-10, rtol=0)
    intervals = np.quantile(fits, [0.025, 0.975], axis=0).T
    np.testing.assert_allclose(intervals, [summary["line_tension_bootstrap_95ci_pN"],
                                          summary["residual_membrane_tension_bootstrap_95ci_pN_per_nm"]],
                               atol=1e-10, rtol=0)
    print(f"Fig. S4: {fit[0]:.6f} pN; 95% bootstrap interval {intervals[0]} pN")


def main():
    for number, folder in [(1, "drop"), (2, "membrane")]:
        stem = f"Data_S{number}"
        with ZipFile(ROOT / "supplementary-data" / f"{stem}.zip") as archive:
            check_catalog(archive, stem + "/")
            locations = [ROOT / "figure_reproduction" / "data" / folder]
            if number == 1:
                locations.append(ROOT / "droplet_analysis" / "data" / folder)
            for location in locations:
                files = sorted(location.glob("*.csv"))
                expected = 2 if number == 1 else 8
                if len(files) != expected:
                    raise ValueError(f"Expected {expected} plotting tables in {location}")
                for path in files:
                    check_projection(path, archive, f"{stem}/{path.name}")
                print(f"Matched {len(files)} tables in {location.relative_to(ROOT)}")
            if number == 2:
                check_lipid_fit(archive)


if __name__ == "__main__":
    main()
