"""Regression tests for identifiable event times and genuine competing-risk math."""
from datetime import date

import pandas as pd
import pytest

from rse.cohort import classify_events
from rse.competing_risks import aalen_johansen, cif_at_times, fine_gray


def row(study_id, *, result=False, posted="", publication=False, stopped=False):
    return {
        "nct_id": study_id,
        "source_as_of_date": "2026-03-29",
        "primary_completion_date_precision": "DAY",
        "results_first_post_date_precision": "DAY" if len(posted) == 10 else "MISSING",
        "results_status_known": True,
        "results_status_conflict": bool(not result and posted),
        "pcd": pd.Timestamp("2023-01-01"),
        "results_first_post_date": posted,
        "has_results": result,
        "pmid_reference_count": 1 if publication else 0,
        "result_reference_count": 0,
        "is_stopped": stopped,
    }


def test_undated_publication_and_stopping_are_censored_not_fake_events():
    df = pd.DataFrame([
        row("NCT00000001", result=True, posted="2023-03-01"),
        row("NCT00000002", publication=True),
        row("NCT00000003", stopped=True),
    ])
    result = classify_events(df, snapshot=date(2026, 3, 29))
    assert dict(zip(result["nct_id"], result["event_code"])) == {
        "NCT00000001": 1, "NCT00000002": 0, "NCT00000003": 0,
    }
    assert result.loc[result["nct_id"] == "NCT00000001", "time_days"].iloc[0] == 59
    assert result.loc[result["nct_id"] == "NCT00000002", "publication_reference_observed"].iloc[0]
    assert result.loc[result["nct_id"] == "NCT00000003", "stopped_without_posted_results"].iloc[0]


def test_undatable_prevalent_and_future_postings_are_excluded_with_audit():
    df = pd.DataFrame([
        row("NCT00000001", result=True, posted="2023-02-01"),
        row("NCT00000002", result=True, posted=""),
        row("NCT00000003", result=True, posted="2022-12-31"),
        row("NCT00000004", result=True, posted="2030-01-01"),
        row("NCT00000005", result=False, posted="2024-01-01"),
    ])
    result = classify_events(df, snapshot=date(2026, 3, 29))
    assert result["nct_id"].tolist() == ["NCT00000001"]
    assert result.attrs["total_excluded"] == 4
    audit = result.attrs["excluded_missing_or_prevalent_results"]
    assert all(n == 1 for n in audit.values())


def test_legacy_proxy_events_warn_and_are_opt_in_only():
    df = pd.DataFrame([
        row("NCT00000001", publication=True),
        row("NCT00000002", stopped=True),
    ])
    with pytest.warns(UserWarning, match="INVALID"):
        result = classify_events(df, legacy_proxy_events=True)
    assert result["event_code"].tolist() == [2, 3]


def test_aj_accounts_for_competing_causes_not_requested_in_cif():
    data = pd.DataFrame({
        "time_days": [1, 2, 3],
        "event_code": [1, 2, 0],
    })
    cif = aalen_johansen(data, event_codes=[1])
    assert cif["times"].tolist() == [1.0, 2.0]
    assert cif["cif"][1].tolist() == pytest.approx([1/3, 1/3])
    assert cif["survival"].tolist() == pytest.approx([2/3, 1/3])


def test_survival_is_one_before_first_event():
    df = pd.DataFrame({"time_days": [30, 60], "event_code": [1, 0]})
    est = cif_at_times(df, [0, 10, 30, 60], event_codes=[1])
    assert est["survival"] == pytest.approx([1.0, 1.0, 0.5, 0.5])
    assert est["cif"][1] == pytest.approx([0.0, 0.0, 0.5, 0.5])


def test_invalid_aj_inputs_rejected():
    with pytest.raises(ValueError):
        aalen_johansen(pd.DataFrame({"time_days": [-1], "event_code": [1]}))
    with pytest.raises(ValueError):
        aalen_johansen(pd.DataFrame({"time_days": [2], "event_code": [1]}), [1, 1])


def test_unvalidated_fine_gray_must_not_emit_hazard_ratios():
    df = pd.DataFrame({"time_days": [1, 2], "event_code": [1, 2], "exposure": [0, 1]})
    with pytest.raises(NotImplementedError, match="arbitrary 0.5"):
        fine_gray(df, ["exposure"])
