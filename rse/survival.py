"""Kaplan-Meier estimation, Cox PH regression, and Schoenfeld diagnostics.

Wraps lifelines for KM curves and Cox proportional hazards models,
with helpers for stratified KM, survival-at-timepoint queries, and
PH assumption testing via Schoenfeld residuals.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter


# Reference levels to match R's factor ordering
_REFERENCE_LEVELS: dict[str, str] = {
    "lead_sponsor_class": "INDUSTRY",
    "phase_label": "Early Phase 1",
}


def fit_km(
    events: pd.DataFrame,
    label: str = "Overall",
    event_of_interest: int = 1,
) -> dict:
    """Fit a Kaplan-Meier survival curve.

    Parameters
    ----------
    events : DataFrame with ``time_days`` and ``event_code`` columns.
    label : curve label for display.
    event_of_interest : event code treated as the event (others censored).

    Returns
    -------
    dict with keys: times, survival, ci_lower, ci_upper, median, label, fitter.
    """
    event_flag = (events["event_code"] == event_of_interest).astype(int)
    kmf = KaplanMeierFitter()
    kmf.fit(events["time_days"], event_observed=event_flag, label=label)

    sf = kmf.survival_function_at_times(kmf.timeline)
    ci = kmf.confidence_interval_survival_function_

    return {
        "times": kmf.timeline.tolist(),
        "survival": sf.values.flatten().tolist(),
        "ci_lower": ci.iloc[:, 0].values.tolist(),
        "ci_upper": ci.iloc[:, 1].values.tolist(),
        "median": kmf.median_survival_time_,
        "label": label,
        "fitter": kmf,
    }


def km_summary_at_times(
    events: pd.DataFrame,
    times: list[int | float],
    event_of_interest: int = 1,
) -> dict:
    """Return KM survival estimates at specific timepoints.

    For each query time *t*, finds the largest observed time <= t in the
    fitted survival function and returns the corresponding survival estimate.

    Parameters
    ----------
    events : DataFrame with ``time_days`` and ``event_code`` columns.
    times : list of timepoints (e.g. [365, 730, 1095, 1825]).
    event_of_interest : event code treated as the event.

    Returns
    -------
    dict with keys: times, survival, ci_lower, ci_upper.
    """
    km = fit_km(events, event_of_interest=event_of_interest)
    kmf = km["fitter"]

    surv_vals = []
    lower_vals = []
    upper_vals = []
    ci_df = kmf.confidence_interval_survival_function_

    for t in times:
        sf_at_t = kmf.predict(t)
        surv_vals.append(float(sf_at_t.iloc[0]) if hasattr(sf_at_t, "iloc") else float(sf_at_t))

        # CI at timepoint: find the closest timeline point <= t
        valid = kmf.timeline[kmf.timeline <= t]
        if len(valid) > 0:
            closest = valid[-1]  # timeline is sorted ascending
            idx = np.searchsorted(ci_df.index.values, closest)
            idx = min(idx, len(ci_df) - 1)
            lower_vals.append(float(ci_df.iloc[idx, 0]))
            upper_vals.append(float(ci_df.iloc[idx, 1]))
        else:
            lower_vals.append(1.0)
            upper_vals.append(1.0)

    return {
        "times": list(times),
        "survival": surv_vals,
        "ci_lower": lower_vals,
        "ci_upper": upper_vals,
    }


def fit_km_stratified(
    events: pd.DataFrame,
    stratum_col: str,
    event_of_interest: int = 1,
) -> dict[str, dict]:
    """Fit KM curves stratified by a categorical column.

    Parameters
    ----------
    events : DataFrame with ``time_days``, ``event_code``, and *stratum_col*.
    stratum_col : column name to stratify by.
    event_of_interest : event code treated as the event.

    Returns
    -------
    dict mapping each stratum level to its KM result dict.
    """
    results = {}
    for level in sorted(events[stratum_col].dropna().unique()):
        subset = events[events[stratum_col] == level]
        results[str(level)] = fit_km(subset, label=str(level),
                                     event_of_interest=event_of_interest)
    return results


def fit_cox(
    events: pd.DataFrame,
    covariates: list[str],
    event_of_interest: int = 1,
) -> dict:
    """Fit a Cox proportional hazards model.

    Categorical covariates are dummy-encoded with ``drop_first=True``,
    using reference levels that match R's ``coxph`` defaults where
    configured in ``_REFERENCE_LEVELS``.

    Parameters
    ----------
    events : DataFrame with ``time_days``, ``event_code``, and covariate columns.
    covariates : list of column names to include as predictors.
    event_of_interest : event code treated as the event.

    Returns
    -------
    dict with keys: coefficients, hazard_ratios, ci_lower, ci_upper,
    p_values, concordance, fitter, summary_df.
    """
    analysis = events[["time_days"] + covariates].copy()
    analysis["event"] = (events["event_code"] == event_of_interest).astype(int)

    for col in covariates:
        if analysis[col].dtype == object or str(analysis[col].dtype) == "category":
            ref = _REFERENCE_LEVELS.get(col)
            levels = sorted(analysis[col].dropna().unique())
            if ref and ref in levels:
                levels.remove(ref)
                levels.insert(0, ref)
            analysis[col] = pd.Categorical(analysis[col], categories=levels)
            dummies = pd.get_dummies(
                analysis[col], prefix=col, prefix_sep="",
                drop_first=True, dtype=float,
            )
            analysis = pd.concat([analysis, dummies], axis=1)
            analysis = analysis.drop(columns=[col])

    analysis = analysis.dropna()

    cph = CoxPHFitter()
    cph.fit(analysis, duration_col="time_days", event_col="event")

    summary = cph.summary
    return {
        "coefficients": summary["coef"].to_dict(),
        "hazard_ratios": summary["exp(coef)"].to_dict(),
        "ci_lower": summary["exp(coef) lower 95%"].to_dict(),
        "ci_upper": summary["exp(coef) upper 95%"].to_dict(),
        "p_values": summary["p"].to_dict(),
        "concordance": cph.concordance_index_,
        "fitter": cph,
        "training_df": analysis,
        "summary_df": summary,
    }


def schoenfeld_test(
    cph_fitter: CoxPHFitter,
    training_df: pd.DataFrame | None = None,
) -> dict[str, float]:
    """Test the proportional hazards assumption via Schoenfeld residuals.

    Uses the rank time-transform and returns a dict mapping each covariate
    name to its p-value.

    Parameters
    ----------
    cph_fitter : a fitted ``CoxPHFitter`` instance.
    training_df : the DataFrame used to fit the model (required by lifelines
        ``proportional_hazard_test``).

    Returns
    -------
    dict mapping covariate names to p-values.
    """
    from lifelines.statistics import proportional_hazard_test

    try:
        if training_df is not None:
            result = proportional_hazard_test(
                cph_fitter, training_df, time_transform="rank",
            )
        else:
            result = proportional_hazard_test(
                cph_fitter, time_transform="rank",
            )
        p_values: dict[str, float] = {}
        for idx in result.summary.index:
            covariate = idx[0] if isinstance(idx, tuple) else idx
            p_values[covariate] = float(result.summary.loc[idx, "p"])
        return p_values
    except Exception:
        return {}
