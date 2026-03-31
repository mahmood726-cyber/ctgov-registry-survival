"""Piecewise Cox regression for interval-specific hazard ratios."""
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

DEFAULT_INTERVALS = [
    (0, 365),
    (365, 730),
    (730, 1460),
    (1460, float("inf")),
]


def _split_to_interval(events, lower, upper, event_of_interest=1):
    """Create episode dataset for a single time interval.

    Subjects at risk during [lower, upper):
    - Exclude subjects who had event or were censored before lower
    - Exit at min(event_time, upper)
    - Event = 1 only if actual event falls within interval
    """
    df = events.copy()
    df = df.loc[df["time_days"] > lower].copy()
    df["t_entry"] = lower
    df["t_exit"] = df["time_days"].clip(upper=upper)
    df["interval_event"] = (
        (df["event_code"] == event_of_interest)
        & (df["time_days"] <= upper)
    ).astype(int)
    df["interval_duration"] = df["t_exit"] - df["t_entry"]
    df = df.loc[df["interval_duration"] > 0]
    return df


def fit_piecewise_cox(events, covariates, intervals=None, event_of_interest=1):
    """Fit separate Cox models for each time interval.

    Returns dict mapping interval label -> result dict with:
    coefficients, hazard_ratios, ci_lower, ci_upper, p_values, n_subjects, n_events.
    If too few events, returns error dict.
    """
    intervals = intervals or DEFAULT_INTERVALS
    results = {}

    for lower, upper in intervals:
        upper_label = f"{int(upper)}" if upper != float("inf") else "inf"
        label = f"{int(lower)}-{upper_label}"

        episode = _split_to_interval(events, lower, upper, event_of_interest)
        n_subjects = len(episode)
        n_events = int(episode["interval_event"].sum())

        if n_events < 5:
            results[label] = {
                "error": f"Too few events ({n_events}) in interval {label}",
                "n_subjects": n_subjects,
                "n_events": n_events,
            }
            continue

        analysis = episode[covariates + ["interval_duration", "interval_event"]].copy()

        for col in covariates:
            if analysis[col].dtype == object or str(analysis[col].dtype) == "category":
                dummies = pd.get_dummies(analysis[col], prefix=col, drop_first=True, dtype=float)
                analysis = pd.concat([analysis, dummies], axis=1)
                analysis = analysis.drop(columns=[col])

        analysis = analysis.dropna()

        try:
            cph = CoxPHFitter()
            cph.fit(analysis, duration_col="interval_duration", event_col="interval_event")
            summary = cph.summary
            results[label] = {
                "coefficients": summary["coef"].to_dict(),
                "hazard_ratios": summary["exp(coef)"].to_dict(),
                "ci_lower": summary["exp(coef) lower 95%"].to_dict(),
                "ci_upper": summary["exp(coef) upper 95%"].to_dict(),
                "p_values": summary["p"].to_dict(),
                "n_subjects": n_subjects,
                "n_events": n_events,
                "concordance": cph.concordance_index_,
            }
        except Exception as exc:
            results[label] = {
                "error": str(exc),
                "n_subjects": n_subjects,
                "n_events": n_events,
            }

    return results
