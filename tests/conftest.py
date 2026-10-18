"""Shared test fixtures for RSE test suite."""
import json
from pathlib import Path

import pandas as pd
import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def synthetic_cohort() -> pd.DataFrame:
    """Load the 1000-row synthetic cohort."""
    df = pd.read_csv(FIXTURES_DIR / "synthetic_cohort.csv")
    # The frozen synthetic fixture was generated for the 2026-03-29 snapshot.
    # Test annotations are explicit so production code need not guess dates.
    df["source_as_of_date"] = "2026-03-29"
    df["primary_completion_date_precision"] = df["primary_completion_date"].astype(str).str.fullmatch(
        r"\d{4}-\d{2}-\d{2}"
    ).map({True: "DAY", False: "PARTIAL"})
    df["results_first_post_date_precision"] = df["results_first_post_date"].astype(str).str.fullmatch(
        r"\d{4}-\d{2}-\d{2}"
    ).map({True: "DAY", False: "MISSING"})
    df["results_status_known"] = True
    df["results_status_conflict"] = (
        ~df["has_results"].astype(bool) & df["results_first_post_date"].notna()
    )
    return df


@pytest.fixture
def r_reference() -> dict:
    """Load R survival package reference values."""
    with open(FIXTURES_DIR / "r_reference.json") as f:
        return json.load(f)
