"""Aalen-Johansen cumulative incidence function + Fine-Gray subdistribution hazard.

Implements the cause-specific Aalen-Johansen (AJ) estimator for competing-risk
CIF curves and a weighted Cox approximation to the Fine-Gray subdistribution
hazard model.

Key formula (AJ):
    CIF_k(t) = sum_{t_j <= t} S(t_j-) * d_k(t_j) / n(t_j)

where S(t_j-) is the overall Kaplan-Meier survival just before event time t_j,
d_k(t_j) is the count of type-k events at t_j, and n(t_j) is the number at risk.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Core Aalen-Johansen estimator
# ---------------------------------------------------------------------------

def aalen_johansen(events: pd.DataFrame, event_codes: list[int] | None = None) -> dict:
    """Estimate the Aalen-Johansen CIF for each event type.

    Parameters
    ----------
    events : DataFrame with columns ``time_days`` (numeric) and
        ``event_code`` (int; 0 = censored, >0 = event type).
    event_codes : list of integer event codes to compute CIFs for.
        Defaults to all non-zero codes present in the data.

    Returns
    -------
    dict with keys:
        times      : np.ndarray of unique event times (ascending).
        cif        : dict mapping each event code -> np.ndarray of CIF values.
        survival   : np.ndarray of overall KM survival S(t) at each time.
        n_at_risk  : np.ndarray of number at risk at each time.
        event_codes: list of event codes included.
    """
    if event_codes is None:
        event_codes = sorted([c for c in events["event_code"].unique() if c != 0])

    times = events["time_days"].values.astype(float)
    codes = events["event_code"].values.astype(int)

    event_mask = codes > 0
    unique_times = np.sort(np.unique(times[event_mask]))

    cif: dict[int, list[float]] = {code: [] for code in event_codes}
    survival_curve: list[float] = []
    n_at_risk_list: list[int] = []

    cum_cif: dict[int, float] = {code: 0.0 for code in event_codes}
    surv = 1.0

    for t in unique_times:
        n_at_risk = int(np.sum(times >= t))
        if n_at_risk == 0:
            for code in event_codes:
                cif[code].append(cum_cif[code])
            survival_curve.append(surv)
            n_at_risk_list.append(0)
            continue

        at_t = times == t
        d_total = 0
        for code in event_codes:
            d_k = int(np.sum(at_t & (codes == code)))
            cum_cif[code] += surv * d_k / n_at_risk
            d_total += d_k

        surv *= 1.0 - d_total / n_at_risk

        for code in event_codes:
            cif[code].append(cum_cif[code])
        survival_curve.append(surv)
        n_at_risk_list.append(n_at_risk)

    return {
        "times": unique_times,
        "cif": {code: np.array(cif[code]) for code in event_codes},
        "survival": np.array(survival_curve),
        "n_at_risk": np.array(n_at_risk_list),
        "event_codes": event_codes,
    }


# ---------------------------------------------------------------------------
# CIF at specific query timepoints
# ---------------------------------------------------------------------------

def cif_at_times(
    events: pd.DataFrame,
    query_times: list[float] | np.ndarray,
    event_codes: list[int] | None = None,
) -> dict:
    """Return AJ CIF estimates at specific timepoints.

    Uses step-function interpolation (last value at or before each query time).

    Parameters
    ----------
    events : DataFrame with ``time_days`` and ``event_code``.
    query_times : sequence of timepoints to query.
    event_codes : event codes to include (defaults to all non-zero).

    Returns
    -------
    dict with keys:
        query_times : list of queried times.
        cif         : dict mapping event code -> list of CIF values at query times.
        survival    : list of overall survival values at query times.
    """
    aj = aalen_johansen(events, event_codes=event_codes)
    event_times = aj["times"]

    def _interp(arr: np.ndarray, t: float) -> float:
        """Step-function lookup: value at largest event_time <= t."""
        idx = np.searchsorted(event_times, t, side="right") - 1
        if idx < 0:
            return 0.0
        return float(arr[idx])

    query_times_arr = np.asarray(query_times, dtype=float)
    result_cif: dict[int, list[float]] = {}
    for code in aj["event_codes"]:
        result_cif[code] = [_interp(aj["cif"][code], t) for t in query_times_arr]

    return {
        "query_times": query_times_arr.tolist(),
        "cif": result_cif,
        "survival": [_interp(aj["survival"], t) for t in query_times_arr],
    }


# ---------------------------------------------------------------------------
# Stratified AJ
# ---------------------------------------------------------------------------

def cif_stratified(
    events: pd.DataFrame,
    stratum_col: str,
    event_codes: list[int] | None = None,
) -> dict[str, dict]:
    """Compute Aalen-Johansen CIF separately for each stratum.

    Parameters
    ----------
    events : DataFrame with ``time_days``, ``event_code``, and *stratum_col*.
    stratum_col : column name to stratify by.
    event_codes : event codes to include.

    Returns
    -------
    dict mapping stratum level (str) -> AJ result dict.
    """
    results: dict[str, dict] = {}
    for level in sorted(events[stratum_col].dropna().unique()):
        subset = events[events[stratum_col] == level].reset_index(drop=True)
        results[str(level)] = aalen_johansen(subset, event_codes=event_codes)
    return results


# ---------------------------------------------------------------------------
# Fine-Gray subdistribution hazard (weighted Cox approximation)
# ---------------------------------------------------------------------------

def fine_gray(
    events: pd.DataFrame,
    covariates: list[str],
    event_of_interest: int = 1,
) -> dict:
    """Fit a Fine-Gray subdistribution hazard model.

    Approximation: subjects who experience a competing event are kept in the
    risk set for the event of interest (extended to the maximum observed time),
    which is the defining feature of subdistribution hazard regression.  A
    weighted Cox PH model is then fitted via lifelines ``CoxPHFitter``.

    Parameters
    ----------
    events : DataFrame with ``time_days``, ``event_code``, and covariate cols.
    covariates : predictor column names.
    event_of_interest : event code for the primary event.

    Returns
    -------
    dict with keys: coefficients, hazard_ratios, ci_lower, ci_upper,
        p_values, n_events, fitter.
    """
    from lifelines import CoxPHFitter

    df = events[["time_days", "event_code"] + covariates].copy().dropna()

    # Determine competing codes (all non-zero codes != event_of_interest)
    all_event_codes = [c for c in df["event_code"].unique() if c != 0]
    competing_codes = [c for c in all_event_codes if c != event_of_interest]

    max_time = float(df["time_days"].max())

    # Build subdistribution risk set:
    #   - primary event subjects: event_flag=1, time as observed
    #   - censored subjects: event_flag=0, time as observed
    #   - competing event subjects: event_flag=0, time extended to max_time
    #     (weight = inverse of KM probability of NOT having had competing event)
    df["sd_time"] = df["time_days"].astype(float)
    df["sd_event"] = (df["event_code"] == event_of_interest).astype(int)

    competing_mask = df["event_code"].isin(competing_codes)
    df.loc[competing_mask, "sd_time"] = max_time

    # Simple weight: 1 for non-competing, 0.5 for competing (conservative approx)
    # A full Fine-Gray uses IPCW weights; this is a practical approximation.
    df["weight"] = 1.0
    df.loc[competing_mask, "weight"] = 0.5

    # Encode categorical covariates
    analysis = df[["sd_time", "sd_event", "weight"] + covariates].copy()
    for col in covariates:
        if analysis[col].dtype == object or str(analysis[col].dtype) == "category":
            analysis[col] = pd.Categorical(analysis[col])
            dummies = pd.get_dummies(
                analysis[col], prefix=col, prefix_sep="_",
                drop_first=True, dtype=float,
            )
            analysis = pd.concat([analysis, dummies], axis=1)
            analysis = analysis.drop(columns=[col])

    cph = CoxPHFitter()
    cph.fit(
        analysis,
        duration_col="sd_time",
        event_col="sd_event",
        weights_col="weight",
        robust=True,
    )

    summary = cph.summary
    n_events = int(df["sd_event"].sum())

    return {
        "coefficients": summary["coef"].to_dict(),
        "hazard_ratios": summary["exp(coef)"].to_dict(),
        "ci_lower": summary["exp(coef) lower 95%"].to_dict(),
        "ci_upper": summary["exp(coef) upper 95%"].to_dict(),
        "p_values": summary["p"].to_dict(),
        "n_events": n_events,
        "fitter": cph,
    }
