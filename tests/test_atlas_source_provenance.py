"""Prevent scientifically invalid censoring dates and inflated trial counts."""
from datetime import date

import pandas as pd
import pytest

from rse.cohort import (
    build_cohort, classify_events, load_atlas, source_snapshot_date,
)


def feature(nct_id="NCT00000001", *, source="2026-03-29",
            primary="2023-01-01", primary_precision="DAY",
            posted="", post_precision="MISSING", has_results=False,
            status_known=True, status_conflict=False):
    return {
        "nct_id": nct_id,
        "source_as_of_date": source,
        "primary_completion_date": primary,
        "primary_completion_date_type": "ACTUAL",
        "primary_completion_date_precision": primary_precision,
        "results_first_post_date": posted,
        "results_first_post_date_precision": post_precision,
        "results_status_known": status_known,
        "results_status_conflict": status_conflict,
        "is_interventional": True,
        "is_closed": True,
        "is_stopped": False,
        "has_results": has_results,
        "pmid_reference_count": 0,
        "result_reference_count": 0,
    }


def test_snapshot_date_comes_from_source_and_rejects_wallclock_override():
    df = pd.DataFrame([feature()])
    assert source_snapshot_date(df) == date(2026, 3, 29)
    assert source_snapshot_date(df, date(2026, 3, 29)) == date(2026, 3, 29)
    with pytest.raises(ValueError, match="does not match"):
        build_cohort(df, snapshot=date(2026, 10, 18))


def test_missing_or_mixed_snapshot_dates_fail_closed():
    df = pd.DataFrame([feature()])
    with pytest.raises(ValueError, match="source_as_of_date"):
        build_cohort(df.drop(columns="source_as_of_date"))
    mixed = pd.DataFrame([feature(), feature("NCT00000002", source="2026-10-18")])
    with pytest.raises(ValueError, match="Mixed"):
        build_cohort(mixed)
    with pytest.raises(ValueError, match="calendar"):
        source_snapshot_date(pd.DataFrame([feature(source="2026-02-30")]))


def test_only_exact_actual_primary_completion_dates_qualify():
    df = pd.DataFrame([
        feature("NCT00000001", primary="2022-01-01", primary_precision="DAY"),
        feature("NCT00000002", primary="2022-01", primary_precision="PARTIAL"),
    ])
    cohort = build_cohort(df)
    assert cohort["nct_id"].tolist() == ["NCT00000001"]


def test_unknown_hasresults_is_not_silently_assumed_false():
    df = pd.DataFrame([
        feature("NCT00000001", status_known=True),
        feature("NCT00000002", status_known=False),
    ])
    cohort = build_cohort(df)
    assert cohort["nct_id"].tolist() == ["NCT00000001"]


def test_partial_results_date_excluded_rather_than_first_of_month_event():
    df = pd.DataFrame([
        feature("NCT00000001", posted="2024-03", post_precision="PARTIAL", has_results=True),
        feature("NCT00000002", posted="2024-03-10", post_precision="DAY", has_results=True),
    ])
    events = classify_events(build_cohort(df))
    assert events["nct_id"].tolist() == ["NCT00000002"]
    assert events["event_code"].tolist() == [1]
    assert events.attrs["total_excluded"] == 1
    assert events.attrs["excluded_missing_or_prevalent_results"]["missing_posted_results_date"] == 1


def test_duplicate_trial_ids_rejected_even_before_merging():
    df = pd.DataFrame([feature(), feature()])
    with pytest.raises(ValueError, match="duplicate NCT IDs"):
        build_cohort(df)


def test_duplicate_condition_rows_rejected(monkeypatch, tmp_path):
    df = pd.DataFrame([feature()])
    cond = pd.DataFrame({
        "nct_id": ["NCT00000001", "NCT00000001"],
        "condition_family": ["cardiovascular", "cardiovascular"]
    })
    source_path = tmp_path / "study_features.parquet"
    condition_path = tmp_path / "condition.parquet"
    condition_path.touch()

    def fake_read(path):
        return cond.copy() if path == condition_path else df.copy()
    monkeypatch.setattr(pd, "read_parquet", fake_read)
    with pytest.raises(ValueError, match="duplicate NCT IDs"):
        load_atlas(source_path, condition_path)


def test_legacy_parquet_without_precision_fields_rejected(monkeypatch, tmp_path):
    df = pd.DataFrame([feature()]).drop(columns=[
        "primary_completion_date_precision", "results_first_post_date_precision"
    ])
    monkeypatch.setattr(pd, "read_parquet", lambda path: df)
    with pytest.raises(ValueError, match="regenerate"):
        load_atlas(tmp_path / "old.parquet", tmp_path / "missing.parquet")


def test_legacy_fixed_snapshot_date_is_not_used_as_default():
    df = pd.DataFrame([feature(source="2026-10-18")])
    cohort = build_cohort(df)
    assert source_snapshot_date(cohort) == date(2026, 10, 18)
    events = classify_events(cohort)
    assert events.loc[0, "time_days"] == (date(2026, 10, 18) - date(2023, 1, 1)).days
