# CT.gov Registry Survival — IMPORTANT SCIENTIFIC CORRECTION (15 October 2026)

**Historical time-to-disclosure estimates and dashboard figures are NOT VALIDATED and must not be used for scientific inference until recomputed.** Three structural errors were discovered in the original analysis:

1. The original competing-risk model treated **any observed publication PMID** as a dated competing event even though its dataset did not include first-publication dates; it allocated the *snapshot date* as the event time.
2. It likewise converted a current terminated/withdrawn registry status plus a 3-year threshold into a competing event at the *snapshot date*, without an actual termination time. That does not establish a valid event history.
3. The function named `fine_gray()` assigned a fixed weight of **0.5** to subjects experiencing competing events. This is **not** Fine–Gray inverse-probability-of-censoring weighting. It is now disabled rather than emitting misleading hazard ratios.

**Corrected method:** `rse.cohort.classify_events` now uses only an observed, dated ClinicalTrials.gov results-first-post event following the ACTUAL primary-completion date; undated publication links and termination flags remain descriptive indicators and do not become invented event times. Studies reporting posted results without an identifiable first-post date, with posting before or on primary completion, or with contradictory/future dates are excluded and counted. Publication timing cannot be recovered from a PMID count. A causal or regulatory claim requires additional data and adjudication. An opt-in `legacy_proxy_events=True` flag is retained *solely* to reproduce historical synthetic R regression checks; its estimates are not valid for publication.

The Aalen–Johansen implementation now subtracts **all** actual event types from the event-free risk set even when a subset of cause-specific CIFs is requested, and returns survival 1 before the first event. Use it for true competing risks **only when cause types AND dates are reliable**. At present the registry features alone support the dated results-posting analysis, with censoring assumptions explicitly stated, not valid publication/termination competing-risk estimates.

To produce revised scientific results: acquire and validate the atlas full-registry release, audit the date exclusions and censoring, recompute all derived plots/tables and independently validate chosen estimands. The historical static dashboard is deliberately left visible with a prominent correction warning; its pre-correction numerical estimates are **superseded pending reanalysis**.

## Project description

Survival analysis of clinical-trial result disclosure on ClinicalTrials.gov.

This repository builds an eligible interventional cohort from the companion
ctgov-hiddenness-atlas study-features dataset, identifies dated registry results-posting events; undated publication and
termination signals are retained as descriptive flags, not competing events,
and runs a corrected cause-specific time-to-posting analysis.

## Modules (`rse/`)

- `cohort.py` — cohort eligibility filtering and competing-event classification.
- `survival.py` — Kaplan-Meier estimation, Cox proportional-hazards regression,
  and Schoenfeld PH-assumption diagnostics (via `lifelines`).
- `competing_risks.py` — Aalen-Johansen cumulative incidence while the previous arbitrary-weight Fine-Gray approximation is disabled.
- `piecewise.py` — piecewise Cox regression for interval-specific hazard ratios.
- `figures.py` — figure helpers.

## Dashboard

`index.html` is a self-contained, offline E156 dashboard. The source data
directory is configured via the `ATLAS_DATA_DIR` environment variable (it
defaults to `data/processed`); no absolute machine paths are hardcoded.

## Tests

```
pip install -r requirements.txt
pytest -q
```

19 tests cover cohort filtering, event classification, KM/Cox fits, competing
risks, and piecewise models against synthetic and R-reference fixtures.
