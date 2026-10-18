"""Cohort builder and event classification for registry survival analysis.

Filters the Hiddenness Atlas study_features dataset to eligible, source-dated
interventional studies. The corrected model uses dated posted-results events
(1) and censoring (0); undated publication/termination proxy events are legacy
only and are not valid for scientific inference.
"""
import os
import re
from datetime import date
from pathlib import Path

import pandas as pd


# This date is ONLY the bundled March-2026 synthetic benchmark's reference.
# Production cohort dates must come from the companion atlas source metadata.
SNAPSHOT_DATE = date(2026, 3, 29)


def _assert_unique_nct_ids(df: pd.DataFrame, source: str) -> None:
    if "nct_id" not in df.columns or df["nct_id"].isna().any():
        raise ValueError(f"{source} has missing NCT IDs; cannot construct a trial-level cohort")
    if df["nct_id"].duplicated().any():
        duplicate = str(df.loc[df["nct_id"].duplicated(), "nct_id"].iloc[0])
        raise ValueError(f"{source} contains duplicate NCT IDs ({duplicate}); refusing inflated study counts")


def source_snapshot_date(df: pd.DataFrame, snapshot: date | None = None) -> date:
    """Bind follow-up to the atlas release, not today's date or a hard-coded era."""
    if "source_as_of_date" not in df.columns:
        raise ValueError(
            "Atlas has no source_as_of_date: regenerate features from a dated "
            "registry source before time-to-event analysis"
        )
    if df.empty or df["source_as_of_date"].isna().any():
        raise ValueError("Atlas source_as_of_date is missing for some/all studies")
    values = set(df["source_as_of_date"].astype(str))
    if len(values) != 1:
        raise ValueError("Mixed CT.gov snapshot dates; cannot define one censoring time")
    raw = next(iter(values))
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", raw):
        raise ValueError("Atlas source_as_of_date must be a day-precise YYYY-MM-DD")
    try:
        observed = date.fromisoformat(raw)
    except ValueError:
        raise ValueError("Atlas source_as_of_date is not a valid calendar date") from None
    if snapshot is not None and observed != snapshot:
        raise ValueError(
            f"Requested censoring date {snapshot} does not match registry "
            f"source_as_of_date {observed}; obtain a matching export"
        )
    return observed


REQUIRED_PRECISION_COLUMNS = (
    "primary_completion_date_precision",
    "results_first_post_date_precision",
    "results_status_known",
    "results_status_conflict",
)


def require_precision_fields(df: pd.DataFrame) -> None:
    absent = [col for col in REQUIRED_PRECISION_COLUMNS if col not in df.columns]
    if absent:
        raise ValueError(
            "Atlas feature table has no auditable date/results precision columns "
            f"({', '.join(absent)}); regenerate the source-date-aware atlas"
        )

# Source data lives in the companion ctgov-hiddenness-atlas project. Resolve the
# atlas root from the ATLAS_DATA_DIR environment variable so the path is not
# hardcoded to one machine; fall back to a sibling "data/processed" directory.
_ATLAS_DATA_DIR = Path(os.environ.get("ATLAS_DATA_DIR", "data/processed"))

ATLAS_PARQUET = _ATLAS_DATA_DIR / "study_features.parquet"
CONDITION_PARQUET = _ATLAS_DATA_DIR / "study_condition_family.parquet"

TERMINATED_THRESHOLD_DAYS = 1095  # 3 years


def load_atlas(path: Path | None = None,
               condition_path: Path | None = None) -> pd.DataFrame:
    """Load the Hiddenness Atlas study_features parquet, enriched with condition families."""
    path = path or ATLAS_PARQUET
    condition_path = condition_path or CONDITION_PARQUET
    df = pd.read_parquet(path)
    _assert_unique_nct_ids(df, "Atlas study-features table")
    source_snapshot_date(df)
    require_precision_fields(df)
    if "condition_family" not in df.columns and condition_path.exists():
        cond = pd.read_parquet(condition_path)[["nct_id", "condition_family"]]
        _assert_unique_nct_ids(cond, "Atlas condition-family table")
        df = df.merge(cond, on="nct_id", how="left", validate="one_to_one")
        df["condition_family"] = df["condition_family"].fillna("other")
    return df


def build_cohort(df: pd.DataFrame, snapshot: date | None = None) -> pd.DataFrame:
    """Filter to eligible interventional cohort."""
    snapshot = source_snapshot_date(df, snapshot)
    require_precision_fields(df)
    _assert_unique_nct_ids(df, "Input cohort")
    mask = (
        (df["is_interventional"] == True)
        & (df["is_closed"] == True)
        & (df["primary_completion_date"].notna())
        & (df["primary_completion_date"].astype(str).str.strip() != "")
        & (df["primary_completion_date_type"] == "ACTUAL")
        & (df["primary_completion_date_precision"] == "DAY")
        & (df["results_status_known"] == True)
        & (df["results_status_conflict"] == False)
    )
    cohort = df.loc[mask].copy()
    cohort["pcd"] = pd.to_datetime(cohort["primary_completion_date"], errors="coerce")
    cohort = cohort.loc[cohort["pcd"].notna() & (cohort["pcd"].dt.date <= snapshot)]
    return cohort.reset_index(drop=True)


def classify_events(cohort: pd.DataFrame, snapshot: date | None = None,
                    terminated_threshold: int = TERMINATED_THRESHOLD_DAYS,
                    *, legacy_proxy_events: bool = False) -> pd.DataFrame:
    """Time from *actual* primary completion to first posted CT.gov results.

    Only the registry results-post date can time this event. Publication PMIDs
    in the atlas do not contain a publication date, and current terminated
    status does not supply a termination date. Treating either as an event at
    the snapshot date creates spurious competing risks.

    In the corrected default analysis, undated publication/stopping indicators
    are retained as covariates/flags, but NOT coded as competing events.
    Records with posted results yet unavailable, future or pre-baseline result
    posting dates are excluded from the dated analysis, with exclusion counts.
    The legacy calculation is available ONLY to reproduce historical synthetic
    R comparisons; it is methodologically invalid for scientific estimates.
    """
    import warnings

    snapshot = source_snapshot_date(cohort, snapshot)
    require_precision_fields(cohort)
    _assert_unique_nct_ids(cohort, "Dated event cohort")
    if not cohort["results_status_known"].all() or cohort["results_status_conflict"].any():
        raise ValueError("Event cohort has unknown or conflicting posted-results status")
    out = cohort.copy()
    out["rpd"] = pd.to_datetime(out["results_first_post_date"], errors="coerce")
    posted_date_precise = out["results_first_post_date_precision"] == "DAY"
    days_to_snapshot = (pd.Timestamp(snapshot) - out["pcd"]).dt.days
    results_lag = (out["rpd"] - out["pcd"]).dt.days

    has_results = out["has_results"].astype(bool)
    has_publication = (
        out["pmid_reference_count"].fillna(0).astype(int) > 0
    ) | (
        out["result_reference_count"].fillna(0).astype(int) > 0
    )
    is_stopped = out["is_stopped"].astype(bool)

    out["publication_reference_observed"] = has_publication
    out["stopped_without_posted_results"] = is_stopped & ~has_results

    if legacy_proxy_events:
        warnings.warn(
            "Legacy undated publication/termination proxy events are INVALID "
            "for registry time-to-disclosure inference; use the default model.",
            UserWarning, stacklevel=2,
        )
        event1 = has_results & results_lag.notna() & (results_lag > 0)
        event2 = ~has_results & has_publication
        event3 = (
            is_stopped & ~has_results & ~has_publication
            & (days_to_snapshot >= terminated_threshold)
        )
        out["event_code"] = 0
        out.loc[event3, "event_code"] = 3
        out.loc[event2, "event_code"] = 2
        out.loc[event1, "event_code"] = 1
        out["time_days"] = days_to_snapshot
        out.loc[event1, "time_days"] = results_lag[event1]
        out["time_days"] = out["time_days"].clip(lower=1)
        labels = {0: "censored", 1: "results_posted", 2: "publication_only_undated_proxy",
                  3: "terminated_undated_proxy"}
    else:
        timed_posting = (
            has_results & posted_date_precise & out["rpd"].notna()
            & (results_lag > 0)
            & (out["rpd"] <= pd.Timestamp(snapshot))
        )
        # A result already posted at/before primary completion is prevalent
        # at time zero, and must not be included in a post-completion risk set.
        missing_time = has_results & (~posted_date_precise | out["rpd"].isna())
        prevalent_posting = has_results & posted_date_precise & (results_lag <= 0)
        future_posting = has_results & posted_date_precise & (out["rpd"] > pd.Timestamp(snapshot))
        contradictory_posting = ~has_results & out["rpd"].notna()
        exclude = (
            missing_time | prevalent_posting | future_posting
            | contradictory_posting
        )
        exclusions = {
            "missing_posted_results_date": int(missing_time.sum()),
            "results_at_or_before_primary_completion": int(prevalent_posting.sum()),
            "results_first_post_date_after_snapshot": int(future_posting.sum()),
            "results_date_present_but_has_results_false": int(contradictory_posting.sum()),
        }
        out = out.loc[~exclude].copy()
        # Series are aligned to study indices after exclusion.
        event1 = timed_posting.loc[out.index]
        out["event_code"] = event1.astype(int)
        out["time_days"] = days_to_snapshot.loc[out.index]
        out.loc[event1, "time_days"] = results_lag.loc[out.index][event1]
        if not (out["time_days"] > 0).all():
            raise ValueError("Time-to-disclosure cohort contains zero or negative follow-up")
        labels = {0: "censored_without_dated_results", 1: "results_posted"}
        out.attrs["excluded_missing_or_prevalent_results"] = exclusions
        out.attrs["total_excluded"] = int(exclude.sum())
        out.attrs["model"] = "dated_results_posting_cause_specific"

    out["event_label"] = out["event_code"].map(labels)
    if "era" not in out.columns:
        out["era"] = pd.cut(
            out["pcd"].dt.year,
            bins=[0, 2006, 2017, 9999],
            labels=["pre-2007", "2007-2017", "post-2017"],
        )
    return out


def cohort_summary(events: pd.DataFrame) -> dict:
    """Return a summary dict of cohort characteristics."""
    n = len(events)
    counts = events["event_code"].value_counts().to_dict()
    return {
        "n_total": n,
        "n_results_posted": counts.get(1, 0),
        "n_publication_only": counts.get(2, 0),
        "n_terminated_no_disclosure": counts.get(3, 0),
        "n_censored": counts.get(0, 0),
        "pct_results_posted": round(100 * counts.get(1, 0) / n, 1) if n > 0 else 0,
        "median_days_to_event": int(events.loc[events["event_code"] == 1, "time_days"].median())
            if counts.get(1, 0) > 0 else None,
    }
