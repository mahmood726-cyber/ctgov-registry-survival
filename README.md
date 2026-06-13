# ctgov-registry-survival

Survival analysis of clinical-trial result disclosure on ClinicalTrials.gov.

This repository builds an eligible interventional cohort from the companion
ctgov-hiddenness-atlas study-features dataset, classifies competing disclosure
events (results posted, publication only, terminated without disclosure, or
censored), and runs time-to-disclosure models.

## Modules (`rse/`)

- `cohort.py` — cohort eligibility filtering and competing-event classification.
- `survival.py` — Kaplan-Meier estimation, Cox proportional-hazards regression,
  and Schoenfeld PH-assumption diagnostics (via `lifelines`).
- `competing_risks.py` — Aalen-Johansen cumulative incidence and a weighted-Cox
  approximation to the Fine-Gray subdistribution hazard.
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
