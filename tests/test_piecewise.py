"""Tests for piecewise Cox regression."""
import numpy as np
import pandas as pd
import pytest

from rse.cohort import build_cohort, classify_events
from rse.piecewise import fit_piecewise_cox, DEFAULT_INTERVALS


class TestPiecewiseCox:
    def test_interval_boundaries_correct(self):
        """Default intervals cover 0 to infinity."""
        assert DEFAULT_INTERVALS[0][0] == 0
        assert DEFAULT_INTERVALS[-1][1] == float("inf")
        for i in range(len(DEFAULT_INTERVALS) - 1):
            assert DEFAULT_INTERVALS[i][1] == DEFAULT_INTERVALS[i + 1][0]

    def test_all_intervals_have_results(self, synthetic_cohort):
        """Each interval should produce a Cox result."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        results = fit_piecewise_cox(events, covariates=["lead_sponsor_class"])
        assert len(results) == len(DEFAULT_INTERVALS)
        for interval, res in results.items():
            assert "coefficients" in res or "error" in res

    def test_covariate_counts_sum_correctly(self, synthetic_cohort):
        """Events across all intervals should sum to total events."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        results = fit_piecewise_cox(events, covariates=["lead_sponsor_class"])
        total_events = sum(r.get("n_events", 0) for r in results.values() if "n_events" in r)
        actual_events = (events["event_code"] == 1).sum()
        assert total_events == actual_events
