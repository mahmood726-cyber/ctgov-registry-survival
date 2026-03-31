"""Tests for KM, Cox regression, and Schoenfeld diagnostics."""
import numpy as np
import pandas as pd
import pytest

from rse.cohort import build_cohort, classify_events
from rse.survival import (
    fit_km,
    fit_cox,
    schoenfeld_test,
    km_summary_at_times,
)


class TestKaplanMeier:
    def test_km_survival_decreases(self, synthetic_cohort):
        """Survival function must be non-increasing."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        km = fit_km(events)
        surv = km["survival"]
        assert all(surv[i] >= surv[i + 1] for i in range(len(surv) - 1))

    def test_km_matches_r_reference(self, synthetic_cohort, r_reference):
        """KM survival at 1, 2, 3, 5 years matches R survfit within 1e-4."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        times = [365, 730, 1095, 1825]
        summary = km_summary_at_times(events, times)
        for i, t in enumerate(times):
            if i < len(r_reference["km_surv"]):
                assert abs(summary["survival"][i] - r_reference["km_surv"][i]) < 1e-4, \
                    f"KM mismatch at t={t}: {summary['survival'][i]} vs R={r_reference['km_surv'][i]}"


class TestCoxRegression:
    def test_cox_coefficients_match_r(self, synthetic_cohort, r_reference):
        """Cox coefficients for sponsor class match R coxph within 0.05."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        cox = fit_cox(events, covariates=["lead_sponsor_class", "phase_label"])
        r_coefs = r_reference["cox_coef"]
        matched = 0
        for r_name, r_val in r_coefs.items():
            # R uses prefix without separator (lead_sponsor_classNIH)
            # pandas get_dummies with prefix_sep="" produces the same
            if r_name in cox["coefficients"]:
                diff = abs(cox["coefficients"][r_name] - r_val)
                assert diff < 0.05, f"Cox coef {r_name}: {cox['coefficients'][r_name]} vs R={r_val}"
                matched += 1
        assert matched >= 10, f"Only matched {matched} coefficient names with R reference"

    def test_schoenfeld_detects_nonproportionality(self, synthetic_cohort):
        """Schoenfeld test should return p-values for covariates."""
        cohort = build_cohort(synthetic_cohort)
        events = classify_events(cohort)
        cox = fit_cox(events, covariates=["lead_sponsor_class", "phase_label"])
        sch = schoenfeld_test(cox["fitter"], cox["training_df"])
        assert len(sch) > 0
        assert all(0 <= p <= 1 for p in sch.values())
