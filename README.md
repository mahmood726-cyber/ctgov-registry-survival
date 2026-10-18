# CT.gov Registry Survival — IMPORTANT SCIENTIFIC CORRECTION (15 October 2026)

**Historical time-to-disclosure estimates and dashboard figures are NOT VALIDATED and must not be used for scientific inference until recomputed.** Three structural errors were discovered in the original analysis:

1. The original competing-risk model treated **any observed publication PMID** as a dated competing event even though its dataset did not include first-publication dates; it allocated the *snapshot date* as the event time.
2. It likewise converted a current terminated/withdrawn registry status plus a 3-year threshold into a competing event at the *snapshot date*, without an actual termination time. That does not establish a valid event history.
3. The function named `fine_gray()` assigned a fixed weight of **0.5** to subjects experiencing competing events. This is **not** Fine–Gray inverse-probability-of-censoring weighting. It is now disabled rather than emitting misleading hazard ratios.

**Corrected method:** `rse.cohort.classify_events` now uses only an observed, dated ClinicalTrials.gov results-first-post event following the ACTUAL primary-completion date; undated publication links and termination flags remain descriptive indicators and do not become invented event times. Studies reporting posted results without an identifiable first-post date, with posting before or on primary completion, or with contradictory/future dates are excluded and counted. Publication timing cannot be recovered from a PMID count. A causal or regulatory claim requires additional data and adjudication. An opt-in `legacy_proxy_events=True` flag is retained *solely* to reproduce historical synthetic R regression checks; its estimates are not valid for publication.

The Aalen–Johansen implementation now subtracts **all** actual event types from the event-free risk set even when a subset of cause-specific CIFs is requested, and returns survival 1 before the first event. Use it for true competing risks **only when cause types AND dates are reliable**. At present the registry features alone support the dated results-posting analysis, with censoring assumptions explicitly stated, not valid publication/termination competing-risk estimates.

To produce revised scientific results: acquire and validate the atlas full-registry release, audit the date exclusions and censoring, recompute all derived plots/tables and independently validate chosen estimands. The historical static dashboard is deliberately left visible with a prominent correction warning; its pre-correction numerical estimates are **superseded pending reanalysis**.

## October 18, 2026 — source-bound survival eligibility

The upstream [CT.gov Hiddenness Atlas](https://github.com/mahmood726-cyber/ctgov-hiddenness-atlas) now preserves the original date precision and propagates a `source_as_of_date` for every record. **This analysis now requires a newly regenerated source table containing those fields**: `source_as_of_date`, `primary_completion_date_precision`, `results_first_post_date_precision`, `results_status_known` and `results_status_conflict`. Older parquet data without this provenance is rejected, not quietly interpreted as day-precise.

- A source release contains one verified date for all rows, which is the **only default censoring date**. Passing another date raises an error. The March 29, 2026 `SNAPSHOT_DATE` constant is retained solely as the known date of the bundled *synthetic* reference dataset.
- Only `ACTUAL`, day-precise `YYYY-MM-DD` primary completion can enter the risk set. Unknown posted-results status is excluded. Posted results with month/year-only dates or invalid/future/pre-baseline dates are excluded and counted; they are **not** assigned arbitrary first-of-month posting times.
- Duplicate NCT IDs in the study table or condition-family lookup are refused. The optional condition join is validated as one-to-one so it cannot silently multiply trial counts.
- The model uses an explicit exact-day parser, avoiding ambiguous mixed-format parsing by pandas. Existing synthetic R comparisons still run with explicitly annotated source dates; they do not validate the real full-registry analysis.

Use the corrected atlas snapshot-generation and analysis scripts first; then point `ATLAS_DATA_DIR` at the newly derived `study_features.parquet`. Run `pytest -q` before regenerating results. **No historical public dashboards or paper figures were recalculated in this change**, and no legal results-reporting conclusions are implied.

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
