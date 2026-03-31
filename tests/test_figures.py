"""Tests for figure generation."""
from pathlib import Path
from xml.etree import ElementTree

import pytest

from rse.cohort import build_cohort, classify_events
from rse.survival import fit_km, fit_cox
from rse.competing_risks import aalen_johansen
from rse.figures import (
    plot_km_curve,
    plot_cox_forest,
    plot_cif_stacked,
)

OUTPUT_DIR = Path(__file__).parent / "test_outputs"


@pytest.fixture(autouse=True)
def clean_output_dir():
    OUTPUT_DIR.mkdir(exist_ok=True)
    yield


class TestFigures:
    def test_km_curve_png_nonempty(self, synthetic_cohort):
        """KM curve produces a non-empty PNG."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        km = fit_km(events)
        path = OUTPUT_DIR / "km_overall.png"
        plot_km_curve(km, output_path=path)
        assert path.exists()
        assert path.stat().st_size > 1000

    def test_km_curve_svg_valid(self, synthetic_cohort):
        """KM curve SVG is valid XML."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        km = fit_km(events)
        path = OUTPUT_DIR / "km_overall.svg"
        plot_km_curve(km, output_path=path, fmt="svg")
        assert path.exists()
        ElementTree.parse(path)

    def test_forest_plot_rows_match_covariates(self, synthetic_cohort):
        """Cox forest plot has one row per covariate."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        cox = fit_cox(events, covariates=["lead_sponsor_class", "phase_label"])
        path = OUTPUT_DIR / "cox_forest.png"
        n_rows = plot_cox_forest(cox, output_path=path)
        assert n_rows == len(cox["coefficients"])
