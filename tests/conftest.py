"""Shared test fixtures for RSE test suite."""
import json
from pathlib import Path

import pandas as pd
import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def synthetic_cohort() -> pd.DataFrame:
    """Load the 1000-row synthetic cohort."""
    return pd.read_csv(FIXTURES_DIR / "synthetic_cohort.csv")


@pytest.fixture
def r_reference() -> dict:
    """Load R survival package reference values."""
    with open(FIXTURES_DIR / "r_reference.json") as f:
        return json.load(f)
