"""Tests for cohort builder and event classification."""
from datetime import date

import numpy as np
import pandas as pd
import pytest

from rse.cohort import build_cohort, classify_events, SNAPSHOT_DATE


class TestBuildCohort:
    """Eligibility filtering tests."""

    def test_filters_to_interventional_closed_actual(self, synthetic_cohort):
        """Only interventional, closed, actual-primary-completion studies pass."""
        cohort = build_cohort(synthetic_cohort)
        assert cohort["is_interventional"].all()
        assert cohort["is_closed"].all()
        assert (cohort["primary_completion_date_type"] == "ACTUAL").all()
        assert cohort["primary_completion_date"].notna().all()

    def test_excludes_future_completion(self, synthetic_cohort):
        """Studies with primary completion after snapshot are excluded."""
        df = synthetic_cohort.copy()
        df.loc[0, "primary_completion_date"] = "2027-01-01"
        cohort = build_cohort(df)
        assert "NCT90000000" not in cohort["nct_id"].values

    def test_cohort_size_within_expected_range(self, synthetic_cohort):
        """Synthetic data is all eligible, so cohort should be ~1000."""
        cohort = build_cohort(synthetic_cohort)
        assert len(cohort) == len(synthetic_cohort)


class TestClassifyEvents:
    """Event classification tests."""

    def test_event_codes_mutually_exclusive(self, synthetic_cohort):
        """Each study gets exactly one event code (1, 2, 3) or 0 (censored)."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        assert events["event_code"].isin([0, 1, 2, 3]).all()

    def test_time_to_event_positive(self, synthetic_cohort):
        """All time-to-event values must be positive."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        assert (events["time_days"] > 0).all()
