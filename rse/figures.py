"""Publication-quality matplotlib figures for registry survival analysis.

Generates KM curves, Cox forest plots, CIF stacked areas, Schoenfeld
residual plots, and piecewise forest plots.  All output uses the Agg
backend so no display is required.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# ── Publication style ────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Helvetica", "Arial"],
    "font.size": 10,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
})

COLORS = [
    "#326891", "#8b6914", "#c03b26", "#2a7f62",
    "#6b4c9a", "#d4831a", "#1a6b8a", "#8c564b",
]

DAYS_PER_MONTH = 30.44

EVENT_LABELS: dict[int, str] = {
    1: "Results posted",
    2: "Publication only",
    3: "Terminated, no disclosure",
}


# ── Helpers ──────────────────────────────────────────────────────────

def _ensure_dir(path: Path) -> None:
    """Create parent directories if needed."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def _save_and_close(fig: plt.Figure, path: Path, fmt: str = "png") -> None:
    """Save figure and close it."""
    _ensure_dir(path)
    fig.savefig(path, format=fmt)
    plt.close(fig)


def _clean_covariate_label(name: str) -> str:
    """Strip R-style covariate prefix to produce a readable label.

    E.g. ``lead_sponsor_classNIH`` -> ``NIH``,
         ``phase_labelPhase 2``    -> ``Phase 2``.
    Handles prefixes defined in survival._REFERENCE_LEVELS keys.
    """
    known_prefixes = ["lead_sponsor_class", "phase_label", "era",
                      "condition_family"]
    for prefix in known_prefixes:
        if name.startswith(prefix) and len(name) > len(prefix):
            return name[len(prefix):]
    return name


# ── 1. Kaplan-Meier curve ────────────────────────────────────────────

def plot_km_curve(
    km_result: dict,
    output_path: str | Path,
    fmt: str = "png",
    title: str = "Kaplan-Meier Survival Curve",
    xlabel: str = "Time (months)",
    ylabel: str = "Survival probability",
) -> None:
    """Plot a single KM curve with CI band, median line, and risk table.

    Parameters
    ----------
    km_result : dict returned by ``fit_km``.
    output_path : file path for the saved figure.
    fmt : image format (``"png"`` or ``"svg"``).
    title, xlabel, ylabel : axis labels.
    """
    output_path = Path(output_path)
    times_months = np.array(km_result["times"]) / DAYS_PER_MONTH
    surv = np.array(km_result["survival"])
    ci_lo = np.array(km_result["ci_lower"])
    ci_hi = np.array(km_result["ci_upper"])
    median = km_result["median"]
    label = km_result.get("label", "Overall")

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Step function
    ax.step(times_months, surv, where="post", color=COLORS[0],
            linewidth=1.5, label=label)
    # CI band
    ax.fill_between(times_months, ci_lo, ci_hi, step="post",
                    alpha=0.15, color=COLORS[0])

    # Median survival line
    if median is not None and np.isfinite(median):
        median_months = median / DAYS_PER_MONTH
        ax.axhline(0.5, color="grey", linestyle="--", linewidth=0.6)
        ax.axvline(median_months, color="grey", linestyle="--", linewidth=0.6)
        ax.annotate(
            f"Median = {median_months:.1f} mo",
            xy=(median_months, 0.5),
            xytext=(median_months + 2, 0.55),
            fontsize=8, color="grey",
            arrowprops=dict(arrowstyle="->", color="grey", lw=0.6),
        )

    ax.set_xlim(left=0)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.legend(loc="lower left", frameon=False)

    # Risk table as text below the plot
    fitter = km_result.get("fitter")
    if fitter is not None:
        tick_months = np.arange(0, times_months.max() + 12, 12)
        tick_days = tick_months * DAYS_PER_MONTH
        n_at_risk = []
        for td in tick_days:
            n_at_risk.append(int((fitter.timeline <= td).sum()
                                 if td == 0
                                 else fitter.event_table["at_risk"].iloc[
                                     max(0, np.searchsorted(
                                         fitter.timeline, td, side="right") - 1)
                                 ]))
        # Place text below axes
        for i, (tm, nr) in enumerate(zip(tick_months, n_at_risk)):
            ax.text(tm, -0.08, str(nr), transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=7, color="grey")
        ax.text(-0.02, -0.08, "At risk", transform=ax.get_yaxis_transform(),
                ha="right", va="top", fontsize=7, color="grey", style="italic")

    fig.tight_layout()
    _save_and_close(fig, output_path, fmt)


# ── 2. Stratified KM curves ─────────────────────────────────────────

def plot_km_stratified(
    km_results_dict: dict[str, dict],
    output_path: str | Path,
    fmt: str = "png",
    title: str = "Stratified Kaplan-Meier Curves",
) -> None:
    """Plot multiple KM curves on the same axes.

    Parameters
    ----------
    km_results_dict : dict mapping stratum label -> KM result dict.
    output_path : file path for the saved figure.
    fmt : image format.
    title : plot title.
    """
    output_path = Path(output_path)
    fig, ax = plt.subplots(figsize=(7, 4.5))

    for i, (stratum, km) in enumerate(km_results_dict.items()):
        color = COLORS[i % len(COLORS)]
        times_months = np.array(km["times"]) / DAYS_PER_MONTH
        surv = np.array(km["survival"])
        ci_lo = np.array(km["ci_lower"])
        ci_hi = np.array(km["ci_upper"])

        ax.step(times_months, surv, where="post", color=color,
                linewidth=1.3, label=stratum)
        ax.fill_between(times_months, ci_lo, ci_hi, step="post",
                        alpha=0.10, color=color)

    ax.set_xlim(left=0)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.legend(loc="lower left", frameon=False, fontsize=8)

    fig.tight_layout()
    _save_and_close(fig, output_path, fmt)


# ── 3. Cox forest plot ──────────────────────────────────────────────

def plot_cox_forest(
    cox_result: dict,
    output_path: str | Path,
    fmt: str = "png",
    title: str = "Cox Proportional Hazards — Hazard Ratios",
) -> int:
    """Plot a horizontal forest plot of hazard ratios from a Cox model.

    Parameters
    ----------
    cox_result : dict returned by ``fit_cox``.
    output_path : file path for the saved figure.
    fmt : image format.
    title : plot title.

    Returns
    -------
    int : number of rows (covariates) plotted.
    """
    output_path = Path(output_path)
    names = list(cox_result["coefficients"].keys())
    hrs = [cox_result["hazard_ratios"][n] for n in names]
    ci_lo = [cox_result["ci_lower"][n] for n in names]
    ci_hi = [cox_result["ci_upper"][n] for n in names]
    pvals = [cox_result["p_values"][n] for n in names]

    n_rows = len(names)
    fig_height = max(3.0, 0.45 * n_rows + 1.5)
    fig, ax = plt.subplots(figsize=(7, fig_height))

    y_positions = list(range(n_rows))

    for i, (name, hr, lo, hi, p) in enumerate(
        zip(names, hrs, ci_lo, ci_hi, pvals)
    ):
        color = COLORS[0] if p < 0.05 else "#999999"
        ax.plot(hr, i, "o", color=color, markersize=6, zorder=3)
        ax.plot([lo, hi], [i, i], "-", color=color, linewidth=1.5, zorder=2)

        # Annotation: HR (CI)
        annot = f"{hr:.2f} ({lo:.2f}\u2013{hi:.2f})"
        ax.text(
            hi * 1.05, i, annot,
            va="center", fontsize=7.5, color=color,
        )

    # Reference line at HR=1
    ax.axvline(1.0, color="grey", linestyle="--", linewidth=0.7, zorder=1)

    # Labels
    clean_names = [_clean_covariate_label(n) for n in names]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(clean_names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Hazard Ratio (95% CI, log scale)")
    ax.set_title(title, fontsize=11, fontweight="bold")

    fig.tight_layout()
    _save_and_close(fig, output_path, fmt)
    return n_rows


# ── 4. CIF stacked area ─────────────────────────────────────────────

def plot_cif_stacked(
    aj_result: dict,
    output_path: str | Path,
    fmt: str = "png",
    title: str = "Cumulative Incidence Functions (Aalen-Johansen)",
) -> None:
    """Plot stacked area CIF chart from Aalen-Johansen results.

    Parameters
    ----------
    aj_result : dict returned by ``aalen_johansen``.
    output_path : file path for the saved figure.
    fmt : image format.
    title : plot title.
    """
    output_path = Path(output_path)
    times_months = aj_result["times"] / DAYS_PER_MONTH
    codes = aj_result["event_codes"]

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Build stacked arrays
    bottoms = np.zeros(len(times_months))
    for i, code in enumerate(codes):
        cif_vals = aj_result["cif"][code]
        label = EVENT_LABELS.get(code, f"Event {code}")
        color = COLORS[i % len(COLORS)]
        ax.fill_between(times_months, bottoms, bottoms + cif_vals,
                        step="post", alpha=0.7, color=color, label=label)
        bottoms = bottoms + cif_vals

    ax.set_xlim(left=0)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Cumulative incidence")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.legend(loc="upper left", frameon=False, fontsize=8)

    fig.tight_layout()
    _save_and_close(fig, output_path, fmt)


# ── 5. Schoenfeld residual plot ──────────────────────────────────────

def plot_schoenfeld(
    cph_fitter,
    covariate: str,
    output_path: str | Path,
    fmt: str = "png",
) -> None:
    """Plot Schoenfeld residuals for a single covariate with LOWESS smooth.

    Parameters
    ----------
    cph_fitter : fitted ``CoxPHFitter`` from lifelines.
    covariate : name of covariate to plot residuals for.
    output_path : file path for the saved figure.
    fmt : image format.
    """
    output_path = Path(output_path)

    try:
        resids = cph_fitter.schoenfeld_residuals_
    except AttributeError:
        # Cannot compute residuals; create placeholder
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "Schoenfeld residuals not available",
                ha="center", va="center", transform=ax.transAxes)
        _save_and_close(fig, output_path, fmt)
        return

    if covariate not in resids.columns:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, f"Covariate '{covariate}' not found",
                ha="center", va="center", transform=ax.transAxes)
        _save_and_close(fig, output_path, fmt)
        return

    times_months = resids.index.values / DAYS_PER_MONTH
    vals = resids[covariate].values

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(times_months, vals, s=8, alpha=0.4, color=COLORS[0],
               edgecolors="none")

    # LOWESS smooth
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        smooth = lowess(vals, times_months, frac=0.3, return_sorted=True)
        ax.plot(smooth[:, 0], smooth[:, 1], color=COLORS[2], linewidth=1.5,
                label="LOWESS")
        ax.legend(frameon=False, fontsize=8)
    except ImportError:
        pass  # statsmodels not available; skip smooth

    ax.axhline(0, color="grey", linestyle="--", linewidth=0.6)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel(f"Schoenfeld residual ({_clean_covariate_label(covariate)})")
    ax.set_title(f"PH Diagnostic: {_clean_covariate_label(covariate)}",
                 fontsize=11, fontweight="bold")

    fig.tight_layout()
    _save_and_close(fig, output_path, fmt)


# ── 6. Piecewise forest plot ────────────────────────────────────────

def plot_piecewise_forest(
    piecewise_results: dict[str, dict],
    covariate_prefix: str,
    output_path: str | Path,
    fmt: str = "png",
) -> None:
    """Plot side-by-side forest plots, one per time interval.

    Parameters
    ----------
    piecewise_results : dict mapping interval label -> result dict
        from ``fit_piecewise_cox``.
    covariate_prefix : prefix to filter covariates (e.g. ``"lead_sponsor_class"``).
    output_path : file path for the saved figure.
    fmt : image format.
    """
    output_path = Path(output_path)

    # Filter to intervals with actual results (no error)
    valid = {k: v for k, v in piecewise_results.items() if "error" not in v}
    n_panels = len(valid)
    if n_panels == 0:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No valid intervals", ha="center", va="center",
                transform=ax.transAxes)
        _save_and_close(fig, output_path, fmt)
        return

    fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, 4),
                             sharey=True, squeeze=False)
    axes = axes.flatten()

    for panel_idx, (interval, res) in enumerate(valid.items()):
        ax = axes[panel_idx]
        # Filter covariates matching the prefix
        cov_names = [k for k in res["hazard_ratios"] if k.startswith(covariate_prefix)]
        if not cov_names:
            cov_names = list(res["hazard_ratios"].keys())

        for i, name in enumerate(cov_names):
            hr = res["hazard_ratios"][name]
            lo = res["ci_lower"][name]
            hi = res["ci_upper"][name]
            p = res["p_values"][name]
            color = COLORS[0] if p < 0.05 else "#999999"

            ax.plot(hr, i, "o", color=color, markersize=5, zorder=3)
            ax.plot([lo, hi], [i, i], "-", color=color, linewidth=1.2, zorder=2)

        ax.axvline(1.0, color="grey", linestyle="--", linewidth=0.6, zorder=1)
        ax.set_xscale("log")
        ax.set_title(f"{interval} days", fontsize=9, fontweight="bold")
        ax.set_xlabel("HR (log)")

        if panel_idx == 0:
            clean = [_clean_covariate_label(n) for n in cov_names]
            ax.set_yticks(range(len(cov_names)))
            ax.set_yticklabels(clean, fontsize=8)
            ax.invert_yaxis()
        else:
            ax.invert_yaxis()

    fig.suptitle(f"Piecewise HRs: {covariate_prefix}", fontsize=11,
                 fontweight="bold", y=1.02)
    fig.tight_layout()
    _save_and_close(fig, output_path, fmt)
