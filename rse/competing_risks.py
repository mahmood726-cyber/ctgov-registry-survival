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
        event_codes = sorted(int(c) for c in events["event_code"].unique() if c != 0)
    elif len(set(event_codes)) != len(event_codes) or any(
        isinstance(c, bool) or not isinstance(c, (int, np.integer)) or c <= 0
        for c in event_codes
    ):
        raise ValueError("event_codes must contain unique positive integers")

    times = events["time_days"].to_numpy(dtype=float)
    raw_codes = events["event_code"].to_numpy(dtype=float)
    if not (np.isfinite(times).all() and np.isfinite(raw_codes).all()):
        raise ValueError("Non-finite follow-up times or event codes")
    if (times < 0).any() or (raw_codes < 0).any() or not np.equal(raw_codes, raw_codes.astype(int)).all():
        raise ValueError("Follow-up times and event codes must be non-negative")
    codes = raw_codes.astype(int)

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
        # All real event types decrease event-free survival, including
        # causes that the caller did not request in this subset of CIFs.
        d_total = int(np.sum(at_t & (codes > 0)))
        for code in event_codes:
            d_k = int(np.sum(at_t & (codes == code)))
            cum_cif[code] += surv * d_k / n_at_risk

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

    def _interp(arr: np.ndarray, t: float, baseline: float = 0.0) -> float:
        """Step-function lookup: value at largest event_time <= t."""
        idx = np.searchsorted(event_times, t, side="right") - 1
        if idx < 0:
            return baseline
        return float(arr[idx])

    query_times_arr = np.asarray(query_times, dtype=float)
    result_cif: dict[int, list[float]] = {}
    for code in aj["event_codes"]:
        result_cif[code] = [_interp(aj["cif"][code], t) for t in query_times_arr]

    return {
        "query_times": query_times_arr.tolist(),
        "cif": result_cif,
        "survival": [_interp(aj["survival"], t, baseline=1.0) for t in query_times_arr],
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

def fine_gray(events: pd.DataFrame, covariates: list[str],
              event_of_interest: int = 1) -> dict:
    """Refuse a misleading Fine-Gray estimate.

    The historical implementation extended competing-event follow-up to the
    last observation and assigned an arbitrary 0.5 weight. This is NOT the
    inverse-probability-of-censoring weighted subdistribution risk set used
    in Fine-Gray regression. There is no scientifically defensible constant
    replacement weight: an independently validated IPCW implementation and
    accurately dated causes would be required.
    """
    raise NotImplementedError(
        "Fine-Gray estimates disabled: historical arbitrary 0.5 "
        "competing-event weights are invalid. Registry publication/termination "
        "dates are also missing. Use cause-specific results-posting Cox/KM "
        "or an Aalen-Johansen CIF only with genuinely dated competing causes."
    )
