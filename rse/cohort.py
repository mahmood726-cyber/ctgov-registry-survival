"""Cohort builder and event classification for registry survival analysis.

Filters the Hiddenness Atlas study_features dataset to eligible interventional
studies and classifies competing events: results posted (1), publication only (2),
terminated without disclosure (3), or censored (0).
"""
from datetime import date
from pathlib import Path

import pandas as pd

SNAPSHOT_DATE = date(2026, 3, 29)

ATLAS_PARQUET = Path(r"C:\Projects\ctgov-hiddenness-atlas\data\processed\study_features.parquet")
CONDITION_PARQUET = Path(r"C:\Projects\ctgov-hiddenness-atlas\data\processed\study_condition_family.parquet")

TERMINATED_THRESHOLD_DAYS = 1095  # 3 years


def load_atlas(path: Path | None = None,
               condition_path: Path | None = None) -> pd.DataFrame:
    """Load the Hiddenness Atlas study_features parquet, enriched with condition families."""
    path = path or ATLAS_PARQUET
    condition_path = condition_path or CONDITION_PARQUET
    df = pd.read_parquet(path)
    if "condition_family" not in df.columns and condition_path.exists():
        cond = pd.read_parquet(condition_path)[["nct_id", "condition_family"]]
        df = df.merge(cond, on="nct_id", how="left")
        df["condition_family"] = df["condition_family"].fillna("other")
    return df


def build_cohort(df: pd.DataFrame, snapshot: date | None = None) -> pd.DataFrame:
    """Filter to eligible interventional cohort."""
    snapshot = snapshot or SNAPSHOT_DATE
    mask = (
        (df["is_interventional"] == True)
        & (df["is_closed"] == True)
        & (df["primary_completion_date"].notna())
        & (df["primary_completion_date"].astype(str).str.strip() != "")
        & (df["primary_completion_date_type"] == "ACTUAL")
    )
    cohort = df.loc[mask].copy()
    cohort["pcd"] = pd.to_datetime(cohort["primary_completion_date"], errors="coerce")
    cohort = cohort.loc[cohort["pcd"].notna() & (cohort["pcd"].dt.date <= snapshot)]
    return cohort.reset_index(drop=True)


def classify_events(cohort: pd.DataFrame, snapshot: date | None = None,
                    terminated_threshold: int = TERMINATED_THRESHOLD_DAYS) -> pd.DataFrame:
    """Classify competing events for each study."""
    snapshot = snapshot or SNAPSHOT_DATE
    out = cohort.copy()
    out["rpd"] = pd.to_datetime(out["results_first_post_date"], errors="coerce")
    results_lag = (out["rpd"] - out["pcd"]).dt.days
    snapshot_dt = pd.Timestamp(snapshot)
    days_to_snapshot = (snapshot_dt - out["pcd"]).dt.days

    has_results = out["has_results"].astype(bool)
    has_publication = (
        out["pmid_reference_count"].fillna(0).astype(int) > 0
    ) | (
        out["result_reference_count"].fillna(0).astype(int) > 0
    )
    is_stopped = out["is_stopped"].astype(bool)

    event1 = has_results & results_lag.notna() & (results_lag > 0)
    event2 = ~has_results & has_publication
    event3 = (is_stopped & ~has_results & ~has_publication & (days_to_snapshot >= terminated_threshold))

    out["event_code"] = 0
    out.loc[event3, "event_code"] = 3
    out.loc[event2, "event_code"] = 2
    out.loc[event1, "event_code"] = 1

    out["time_days"] = days_to_snapshot
    out.loc[event1, "time_days"] = results_lag[event1]
    out["time_days"] = out["time_days"].clip(lower=1)

    event_labels = {0: "censored", 1: "results_posted", 2: "publication_only",
                    3: "terminated_no_disclosure"}
    out["event_label"] = out["event_code"].map(event_labels)

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
