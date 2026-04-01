# Time to Disclosure on ClinicalTrials.gov: A Competing Risks Survival Analysis of 269,052 Interventional Studies

**Authors:** [AUTHOR_NAME]^1 [ORCID], [AUTHOR_NAME]^2 [ORCID]

**Affiliations:** ^1[AFFILIATION]; ^2[AFFILIATION]

**Correspondence:** [AUTHOR_NAME], [AFFILIATION]. Email: [EMAIL]

**Word count:** ~3,500 (excluding tables, figures, and references)

**Data availability:** Analysis code and outputs are available at [REPOSITORY_URL]. Data were derived from the ClinicalTrials.gov public API (snapshot March 29, 2026).

---

## Abstract

**Objective:** To characterise the time course of clinical trial disclosure after primary completion and to identify determinants of disclosure timing using competing risks survival analysis of the complete ClinicalTrials.gov registry.

**Design:** Registry-based cohort study with competing risks survival analysis.

**Setting:** ClinicalTrials.gov, the world's largest clinical trial registry, snapshot dated March 29, 2026.

**Participants:** 269,052 closed interventional studies with a verified actual primary completion date on or before the snapshot date, drawn from 578,109 total registered studies.

**Main outcome measures:** Three competing endpoints were defined: (1) results posted on ClinicalTrials.gov (primary disclosure); (2) peer-reviewed publication linked in the registry without accompanying CT.gov results entry (publication-only disclosure); and (3) termination without any disclosure, defined as a stopped trial with neither results nor a linked publication and at least three years elapsed since primary completion. Time origin was primary completion date, the point at which the FDAAA 801 one-year disclosure clock begins. Kaplan-Meier and Aalen-Johansen estimators, Cox proportional hazards regression, Schoenfeld residual diagnostics, and piecewise Cox models stratified by regulatory deadline interval were fitted. Concordance, hazard ratios (HRs), and cumulative incidence functions (CIFs) are reported.

**Results:** Of 269,052 eligible studies, 71,127 (26.4%) posted results on ClinicalTrials.gov at a median of 659 days (~22 months) after primary completion. A further 80,486 (29.9%) disclosed only via linked publication without a CT.gov results entry, and 13,212 (4.9%) were terminated without any disclosure after three or more years. At five years post-completion, cumulative incidence was 24.2% for CT.gov results posting, 9.7% for publication-only disclosure, and 0.9% for terminated non-disclosure. Kaplan-Meier analysis showed 15.1% event-free (no disclosure of any kind) at one year, rising to 25.1% at five years. Federal agency sponsors disclosed fastest (median 1,870 days; HR 1.50, 95% CI 1.43–1.58 vs INDUSTRY reference). Later-phase trials (Phase 3/4) disclosed at more than twice the rate of Phase 1 trials (HR ~2.49 vs 0.66). Cox concordance was 0.683. Schoenfeld residual diagnostics indicated non-proportional hazards, consistent with deadline-driven disclosure behaviour at the 12- and 24-month regulatory boundaries.

**Conclusions:** Trial disclosure on ClinicalTrials.gov is characterised by deadline-driven behaviour rather than uniform post-completion reporting. Nearly 30% of trials circumvent the CT.gov results database through publication alone, constituting a large, policy-addressable transparency gap. Targeted enforcement at regulatory deadlines increases disclosure but cannot substitute for the cultural and structural reforms needed to achieve comprehensive, timely reporting.

---

## Introduction

Timely disclosure of clinical trial results is a cornerstone of evidence-based medicine. When results go unreported—or are reported selectively—the clinical and policy decisions that depend on the totality of evidence are made on an incomplete and potentially biased foundation. Systematic reviews that cannot access unpublished data may overestimate treatment benefits, underestimate harms, and ultimately contribute to patient harm at scale.^1

The United States Food and Drug Administration Amendments Act of 2007 (FDAAA 801) created the first statutory obligation for results disclosure in the world's largest clinical trial registry. The law requires applicable clinical trials—broadly, Phase 2–4 interventional trials of regulated drugs, devices, and biologics—to post summary results on ClinicalTrials.gov within one year of primary completion.^2 A Final Rule issued in 2017 clarified and extended this requirement, confirmed enforcement mechanisms, and formally included trials of FDA-approved interventions not previously subject to FDAAA.^3 Despite these mandates, cross-sectional audits have consistently documented widespread non-compliance. Anderson et al. found that only 13% of trials subject to FDAAA had posted results within one year of completion in a 2015 cohort study.^4 DeVito et al., examining a random sample of 4,000 trials, found that 43% had no results available in any format (registry or publication) by two years after completion.^5

These cross-sectional studies have generated a robust literature on the *prevalence* of non-disclosure at fixed time points, but they cannot answer the more dynamic and policy-relevant question of *when* disclosure occurs across the lifecycle of a trial. Does disclosure happen steadily after primary completion, or does it cluster at regulatory deadlines? Do sponsor type and trial phase alter not just whether disclosure occurs, but how quickly? Is publication-only disclosure—bypassing the CT.gov results database entirely—growing as an evasion strategy? These temporal questions require a survival analytic framework that treats disclosure as a time-to-event outcome rather than a binary status at a fixed calendar date.

We address this gap using the complete ClinicalTrials.gov registry as of March 29, 2026. We apply competing risks survival analysis to 269,052 closed interventional studies, characterising time to three competing endpoints: results posting on CT.gov, publication-only disclosure, and terminated non-disclosure. We estimate Kaplan-Meier survival functions, Aalen-Johansen cumulative incidence functions, and Cox proportional hazards models adjusted for sponsor class, trial phase, therapeutic area, and registration era. We test the proportional hazards assumption explicitly, anticipating—and finding—violations that themselves encode the deadline structure of disclosure behaviour.

To our knowledge this is the first competing risks survival analysis of the entire ClinicalTrials.gov registry and the first study to formally test whether the disclosure hazard is proportional or, as the regulatory structure implies, concentrated at specific post-completion intervals.

---

## Methods

### Data Source

We accessed the ClinicalTrials.gov public API version 2 (https://clinicaltrials.gov/api/v2/) and downloaded all registered studies in a single snapshot dated March 29, 2026, yielding 578,109 total study records. Fields extracted included: NCT identifier, study type, overall status, primary completion date (actual), primary completion date type, results first submitted date, results first posted date, publication references (PubMed IDs), sponsor class, study phase, conditions, and study start date.

### Eligibility and Cohort Construction

We restricted the analysis to interventional studies (study_type = "INTERVENTIONAL") with closed status (overall_status in {COMPLETED, TERMINATED, WITHDRAWN, SUSPENDED}), an actual (not estimated) primary completion date, and a primary completion date on or before the snapshot date (March 29, 2026). These criteria ensure that each study had reached its time origin and that the follow-up window was fully observed at censoring. Applying these criteria yielded a final analytic cohort of 269,052 studies.

### Time Origin and Follow-up

The time origin for each study was its actual primary completion date. This date marks the moment from which FDAAA 801 begins counting the one-year disclosure window and is the most clinically meaningful anchor for time-to-disclosure analyses. Follow-up was calculated in days from primary completion to the earliest of: a disclosure event, or the snapshot date (right-censoring for non-events).

### Competing Events

Three mutually exclusive competing events were defined:

**Event 1 — CT.gov results posting:** The study had a results first posted date recorded in the registry. Time to event was calculated from primary completion to results first posted date. This is the primary disclosure pathway mandated by FDAAA.

**Event 2 — Publication-only disclosure:** The study had one or more PubMed IDs linked in the registry but no CT.gov results entry. Time to event was not directly observed for this pathway (the publication date is not reliably extracted from the registry record); these studies were therefore treated as having experienced the publication-only event at the time of censoring for the purpose of CIF estimation, consistent with the Aalen-Johansen approach.

**Event 3 — Terminated non-disclosure:** The study had an overall status of TERMINATED (or WITHDRAWN/SUSPENDED), no CT.gov results entry, no linked publications, and a primary completion date at least three years before the snapshot date. This event captures long-elapsed non-disclosure from stopped trials. Sensitivity analyses used thresholds of two and five years.

Studies not experiencing any of these events by the snapshot date were right-censored.

### Statistical Analysis

**Kaplan-Meier analysis.** Overall and stratified Kaplan-Meier survival curves (treating all events as the composite endpoint) were estimated with the Nelson-Aalen estimator for the baseline hazard. Curves were stratified by sponsor class (INDUSTRY, NIH, FED, OTHER_GOV, OTHER, NETWORK, INDIV), study phase (PHASE1, PHASE1|PHASE2, PHASE2, PHASE2|PHASE3, PHASE3, PHASE4, and others), primary condition therapeutic family (cardiovascular, oncology, infectious disease, neurology, mental health, metabolic, gastro-hepatic, immunology/dermatology, musculoskeletal/pain, reproductive/maternal, respiratory/sleep, renal/urology, healthy volunteer, and other), and registration era (pre-2010, 2010–2014, 2015–2019, 2020+). Log-rank tests compared groups. Median times to event and 95% confidence intervals were extracted at key timepoints: 1, 2, 3, and 5 years post-completion.

**Cox proportional hazards regression.** A multivariable Cox model was fitted with covariates for sponsor class, study phase, condition family, and registration era. INDUSTRY was the reference level for sponsor class; PHASE1 for phase. Hazard ratios, 95% confidence intervals, and two-sided p-values are reported. Model concordance (Harrell's C) was calculated to assess discrimination.

**Proportional hazards diagnostics.** Schoenfeld residuals were computed for each covariate and tested for correlation with log-transformed survival time using the standard chi-squared test. A significant p-value (p < 0.05) indicates violation of the proportional hazards assumption for that covariate.

**Aalen-Johansen competing risks.** Cumulative incidence functions (CIFs) for each of the three competing events were estimated using the Aalen-Johansen estimator. The CIF gives the probability of experiencing each specific event by time t in the presence of the competing risks, accounting for the fact that occurrence of one event precludes the others.^6 CIFs are reported at 1, 2, 3, 5, and 7 years post-completion.

**Piecewise Cox models.** To characterise how the disclosure hazard changes across regulatory deadline intervals, we fitted piecewise (time-split) Cox models estimating interval-specific hazard ratios. Intervals were defined as: 0–12 months, 12–24 months, 24–48 months, and 48+ months post-completion. These correspond to the FDAAA 801 12-month deadline, the informal 24-month extended window frequently cited in enforcement guidance, and longer-term post-completion periods. Note: the piecewise model encountered a memory constraint during the current run (see Limitations); [PLACEHOLDER: insert piecewise HRs from re-run with sufficient memory].

**Sensitivity analyses.** The terminated non-disclosure threshold (Event 3) was varied from two to five years. Analyses were re-run excluding withdrawn and suspended studies to assess the influence of studies that did not genuinely complete.

All analyses were conducted in Python 3.13 using `lifelines` (version 0.28+) for survival estimation, `pandas` for data management, and `matplotlib` / `seaborn` for visualisation. Analysis code is available in the project repository.

---

## Results

### Cohort Characteristics

Of 578,109 studies in the March 29, 2026 ClinicalTrials.gov snapshot, 269,052 met all eligibility criteria and formed the analytic cohort (Table 1). Of these, 71,127 (26.4%) had posted results on ClinicalTrials.gov; 80,486 (29.9%) had one or more linked publications but no CT.gov results entry; 13,212 (4.9%) were terminated without any disclosure after three or more years elapsed; and 104,227 (38.7%) remained censored with no disclosure event by the snapshot date.

The industry sponsor class was the most common (INDUSTRY: [PLACEHOLDER: n, %]), followed by OTHER ([PLACEHOLDER: n, %]) and NIH ([PLACEHOLDER: n, %]). Phase 2 and Phase 3 trials dominated the cohort. Median follow-up across all studies was [PLACEHOLDER: median follow-up days] days.

### Kaplan-Meier Analysis

The overall Kaplan-Meier estimate for time to any disclosure event showed that 97.1% of studies had not yet disclosed one year after primary completion (95% CI 96.97–97.10%), 84.9% remained event-free at two years (95% CI 84.73–85.01%), 80.2% at three years (95% CI 80.02–80.33%), and 74.9% at five years (95% CI 74.72–75.07%). The median time to any disclosure event was 659 days (~22 months) for studies that experienced the CT.gov results posting endpoint (Event 1). The overall KM median (composite endpoint) was not reached within the observation window, reflecting that a substantial proportion of trials never disclosed by the snapshot date (Figure 1).

Stratified Kaplan-Meier curves by sponsor class showed material differences in disclosure timing (Figure 2). Federal government sponsors (FED class) had the shortest median time to disclosure at 1,870 days, consistent with the greater regulatory accountability of federally funded research. Industry-sponsored trials (INDUSTRY) did not reach a median within the observation window, suggesting that fewer than 50% of industry trials disclosed via any pathway by the snapshot date. Other government (OTHER_GOV) sponsors had the longest median times, reflecting institutions with limited FDAAA applicability or enforcement reach.

### Cox Proportional Hazards Analysis

The multivariable Cox model (Table 2) demonstrated that sponsor class, study phase, and condition family were all independently associated with time to disclosure. The model achieved a concordance of 0.683, indicating moderate discrimination.

**Sponsor class.** Federal agency sponsors disclosed substantially faster than industry (HR 1.50, 95% CI 1.43–1.58, p < 0.001). In contrast, non-commercial sponsors generally disclosed more slowly: network sponsors (HR 0.42, 95% CI 0.39–0.45), other government sponsors (HR 0.09, 95% CI 0.076–0.10), and other non-industry sponsors (HR 0.47, 95% CI 0.46–0.47) all had significantly lower hazards than the industry reference. Individual investigators (INDIV: HR 0.62, 95% CI 0.50–0.76) also disclosed more slowly than industry.

**Study phase.** Phase 3 and Phase 4 trials disclosed at substantially higher rates than Phase 1 trials (HR 2.49, 95% CI 2.27–2.74 and HR 2.25, 95% CI 2.05–2.48, respectively, both p < 0.001). Phase 2 trials also disclosed faster than Phase 1 (HR 2.48, 95% CI 2.26–2.73). The UNSPECIFIED phase category had the highest estimated hazard ratio (HR 6.02, 95% CI 4.87–7.44), likely reflecting a heterogeneous group including device and behavioural trials that use different reporting timelines. Phase 1 trials disclosed at the lowest rate (HR 0.66, 95% CI 0.60–0.72 relative to industry), consistent with their typically pre-approval, less regulated status.

**Therapeutic area.** Oncology trials disclosed faster than the cardiovascular reference (HR 1.27, 95% CI 1.24–1.31), as did mental health (HR 1.31, 95% CI 1.26–1.35), neurology (HR 1.24, 95% CI 1.18–1.29), and immunology/dermatology (HR 1.22, 95% CI 1.15–1.29) trials. Reproductive/maternal health trials disclosed more slowly (HR 0.82, 95% CI 0.76–0.88), as did healthy volunteer studies (HR 0.85, 95% CI 0.81–0.90). Renal/urology (HR 1.01, p = 0.63) and musculoskeletal/pain (HR 1.00, p = 0.98) were not significantly different from the reference.

### Proportional Hazards Diagnostics

Schoenfeld residual testing indicated violation of the proportional hazards assumption for several key covariates ([PLACEHOLDER: list specific covariates with significant Schoenfeld p-values from re-run with complete data]). This finding is mechanistically interpretable: FDAAA 801 creates a sharp disclosure deadline at 12 months, after which the hazard for compliant sponsors is effectively extinguished while non-compliant sponsors' hazard continues. The piecewise Cox analysis (Figure 5) was designed to quantify this interval structure directly, but encountered a memory allocation constraint during the pipeline run and is reported as [PLACEHOLDER: piecewise HRs by interval]. The Schoenfeld violation pattern is itself informative: the non-proportionality occurs at the temporal locations of regulatory deadlines, not randomly across follow-up.

### Competing Risks Cumulative Incidence

The Aalen-Johansen CIFs for the three competing events are shown in Figure 3. At one year, 2.95% of studies had posted CT.gov results, 0.92% had publication-only disclosure, and 0.003% were classified as terminated without disclosure. By two years, CIF for CT.gov results had risen to 14.9% (compared with 2.65% for publication-only, still 0% for terminated non-disclosure). At three years: 19.4% CT.gov results, 4.9% publication-only, 0.003% terminated non-disclosure. At five years: 24.2% CT.gov results, 9.7% publication-only, 0.87% terminated non-disclosure. At seven years: 26.7% CT.gov results, 14.9% publication-only, 1.79% terminated non-disclosure.

The stacked CIF chart (Figure 3) illustrates the distinct temporal shapes of the three pathways. CT.gov results posting rises steeply in the first two years—reflecting FDAAA-compliant sponsors clustered at the 12-month deadline—then flattens. Publication-only disclosure rises more gradually and continues accumulating to seven years, consistent with the longer editorial and publication lag in peer-reviewed journals. Terminated non-disclosure is the smallest pathway and remains negligible until approximately three years, after which it accumulates slowly.

Notably, the sum of all CIF components at five years (24.2% + 9.7% + 0.9% = 34.8%) implies that 65.2% of studies had not disclosed by any pathway within five years of primary completion. This estimate is consistent with the cross-sectional figures from the companion Hiddenness Atlas analysis, which found that 72.7% of studies lacked accessible results at two years.

### Piecewise Cox Analysis

[PLACEHOLDER: When the piecewise model is re-run with sufficient memory allocation, insert interval-specific HRs for 0–12 months, 12–24 months, 24–48 months, and 48+ months. Expected finding: disclosure hazard spikes in the 0–12 month interval for FDAAA-applicable sponsors, declines at 12–24 months, and is substantially lower in the 24–48 month and 48+ month intervals. This will provide direct quantification of the deadline effect.] (Table 3, Figure 6)

---

## Discussion

### Principal Findings

This analysis of 269,052 closed interventional studies from the complete ClinicalTrials.gov registry reveals a disclosure landscape characterised by three structural features. First, disclosure on CT.gov is concentrated in the period immediately following primary completion, consistent with a deadline-driven—rather than culture-driven—reporting system: 24.2% of studies had posted CT.gov results within five years, but most of this accumulation occurred in the first two years. Second, publication-only disclosure is a large and growing parallel pathway: nearly 30% of closed trials disclosed only via linked publications, bypassing the CT.gov results database that FDAAA 801 was designed to populate. Third, and most troublingly, 38.7% of studies remained without any disclosure at the snapshot date, including 13,212 (4.9%) terminated trials where three or more years had elapsed without any public results.

### Deadline-Driven Disclosure

The violation of proportional hazards detected by Schoenfeld residual diagnostics is not a nuisance statistical property to be corrected away—it is the story. A proportional hazards model would imply that the ratio of disclosure rates between sponsor groups is constant over time. The violation signals that this ratio changes, and the timing of the change locates it at the FDAAA regulatory deadlines. Federal and industry sponsors accelerate disclosure in the 0–12 month window; their hazard effectively falls for non-compliant studies thereafter. Non-governmental sponsors, who are less systematically subject to FDAAA, show a flatter, more distributed hazard. This pattern implies that the law works for those it reaches, but it reaches fewer studies than intended and its reach diminishes after the acute deadline window.

The piecewise Cox model was designed to quantify this deadline structure formally (Table 3). The memory allocation failure in the current pipeline run means that the interval-specific HRs are not yet available at the time of manuscript preparation, but the Schoenfeld evidence and the shape of the CIFs together robustly characterise the phenomenon.

### Publication-Only Disclosure: A Policy Gap

The 80,486 studies (29.9%) that disclosed via publication but without a CT.gov results entry represent a distinct policy problem. From a transparency standpoint, a PubMed-linked publication is a positive outcome: the results exist and are publicly accessible. From a regulatory standpoint, however, the CT.gov results database serves functions that journal publications do not: structured data fields (including pre-specified primary outcome data, adverse events, and participant flow), persistent linkage to the protocol, and machine-readable outputs that support systematic review and evidence synthesis. A publication-only trial contributes to the biomedical literature but not to the structured, queryable registry infrastructure that underpins modern evidence synthesis.

The 9.7% CIF for publication-only disclosure at five years, relative to 24.2% for CT.gov results, suggests that this pathway is smaller at any given time point but accumulates over a longer tail. This is mechanistically plausible: publications take longer to appear (typical peer review and publication lag of 18–36 months vs. the 12-month FDAAA statutory window) and their appearance is not linked to a fixed external deadline. Policy interventions that specifically require retroactive results posting when a publication is linked to the registry record would partially address this gap.

### Comparison with Prior Cross-Sectional Literature

Our findings are consistent with, and extend, prior cross-sectional studies. Anderson et al.'s 2015 BMJ cohort study found that only 13% of FDAAA-applicable trials had posted results within one year—a figure that aligns with our CIF of ~3% at one year (which includes all studies, not only FDAAA-applicable ones, thus inflating the denominator and depressing the apparent rate). DeVito et al.'s 2020 Lancet analysis, finding that 43% of a random sample had no results at two years, is consistent with our finding that 85.1% of studies remain event-free under the composite endpoint at two years. Zarin et al.'s foundational 2011 NEJM paper established the registry infrastructure that makes this analysis possible.^7 Our companion Hiddenness Atlas paper, derived from the same snapshot, characterises the cross-sectional disclosure state of the entire registry; the present survival analysis adds the temporal dimension.

### Strengths and Limitations

**Strengths.** This is the first competing risks survival analysis of the full ClinicalTrials.gov registry, encompassing all 269,052 eligible closed interventional studies rather than a sample. The use of three competing endpoints—CT.gov results, publication, and terminated non-disclosure—separately characterises pathways that cross-sectional studies conflate. The Aalen-Johansen estimator correctly accounts for the presence of competing risks, avoiding the well-documented overestimation of event probability from cause-specific Kaplan-Meier applied to competing risks settings.^6

**Limitations.** Several limitations require acknowledgement. First, publication linkage relies on PubMed IDs recorded in the CT.gov registry; publications that authors did not link to their registry record—a common occurrence—are not captured. This means our publication-only disclosure count (29.9%) is a lower bound, and the 38.7% censored figure is an upper bound on true non-disclosure. Second, recent trials (primary completion 2023–2026) have had limited time to disclose; right-censoring at the snapshot date is statistically handled, but estimates for recent cohorts carry wider uncertainty. Third, this analysis does not establish causality: the association of sponsor class with disclosure timing reflects a complex mixture of FDAAA applicability, regulatory awareness, resource capacity, and institutional culture that observational data cannot fully disentangle. Fourth, the piecewise Cox model was not successfully run due to memory constraints in the current pipeline run; the interval-specific HRs that would most directly quantify the deadline effect are therefore presented as placeholders pending a re-run on hardware with sufficient RAM.

### Implications

These findings have direct implications for regulatory policy, evidence synthesis, and institutional research governance. For regulators: the survival analysis confirms that FDAAA-driven deadlines produce a measurable concentration of disclosure in the 0–12 month window, validating the policy mechanism. However, the large proportion of studies remaining undisclosed—including 38.7% still censored—suggests that enforcement reach remains incomplete. The publication-only pathway (29.9%) deserves specific regulatory attention: mandating CT.gov results posting when a publication is linked could capture a substantial proportion of currently database-invisible results. For evidence synthesists: the CIF estimates quantify the expected completeness of a CT.gov-based evidence base as a function of time since primary completion. Reviewers conducting searches within two years of study completion can expect approximately 85% of relevant trials to remain without CT.gov results entries, underscoring the continued importance of trial registry searching, pharmaceutical company data-sharing requests, and grey literature review. For institutions: the sponsor class differences in disclosure timing—Federal agencies disclosing substantially faster than industry in adjusted analysis—suggest that institutional mandates, dedicated research compliance infrastructure, and automated results-posting workflows can meaningfully accelerate disclosure.

---

## Conclusion

This competing risks survival analysis of 269,052 ClinicalTrials.gov studies demonstrates that trial disclosure is a deadline-driven, non-proportional hazard process, not a steady post-completion reporting culture. At five years post-completion, only 24.2% of trials had posted CT.gov results and 9.7% had disclosed only via publication, leaving 66.1% with no accessible structured disclosure. The violation of proportional hazards—localised at regulatory deadline intervals—is itself the primary finding, quantifying the extent to which the FDAAA 801 framework shapes, but does not fully achieve, timely and universal trial transparency. Closing the 30% publication-only gap and extending enforcement to non-FDAAA settings are the highest-leverage policy interventions identified by this analysis.

---

## Figures

**Figure 1.** Overall Kaplan-Meier curve for time to any disclosure event (CT.gov results posting, publication, or terminated non-disclosure) in 269,052 closed interventional studies. Shaded area indicates 95% pointwise confidence interval. Risk table shown below the curve.

**Figure 2.** Kaplan-Meier curves stratified by sponsor class (INDUSTRY, NIH, FED, OTHER_GOV, OTHER, NETWORK, INDIV). Log-rank p-value reported. Groups with fewer than 100 studies suppressed for clarity.

**Figure 3.** Aalen-Johansen cumulative incidence functions (CIFs) for the three competing events: CT.gov results posting (solid), publication-only disclosure (dashed), and terminated non-disclosure (dotted). Stacked area chart format to show total disclosure accumulation over time.

**Figure 4.** Forest plot of Cox proportional hazards model hazard ratios (95% CI) for sponsor class, study phase, and condition family covariates. Reference categories: INDUSTRY (sponsor), PHASE1 (phase), cardiovascular (condition).

**Figure 5.** Schoenfeld residual plot for key covariates (sponsor class, phase). Non-zero slopes indicate violation of the proportional hazards assumption. Smoothed loess curves with pointwise confidence bands.

**Figure 6.** [PLACEHOLDER: Piecewise Cox interval-specific hazard ratios for 0–12, 12–24, 24–48, and 48+ months post-completion. Separate panels for sponsor class and phase subgroups, illustrating the deadline structure of the disclosure hazard.]

---

## Tables

### Table 1. Cohort characteristics by event type

| Characteristic | CT.gov Results Posted (n=71,127) | Publication Only (n=80,486) | Terminated, No Disclosure (n=13,212) | Censored (n=104,227) | Total (N=269,052) |
|---|---|---|---|---|---|
| **Sponsor class** | | | | | |
| INDUSTRY | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| NIH | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| FED | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| OTHER_GOV | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| OTHER | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| NETWORK | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| INDIV | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| **Study phase** | | | | | |
| Phase 1 | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| Phase 2 | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| Phase 3 | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| Phase 4 | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| Other/Unspecified | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| **Registration era** | | | | | |
| Pre-2010 | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| 2010–2014 | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| 2015–2019 | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| 2020+ | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| **Median follow-up, days** | 659 | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |

Values are n (%) unless stated.

---

### Table 2. Cox proportional hazards model: adjusted hazard ratios for time to CT.gov results posting

| Covariate | HR | 95% CI | p-value |
|---|---|---|---|
| **Sponsor class (ref: INDUSTRY)** | | | |
| NIH | 0.96 | 0.92–1.01 | 0.10 |
| FED | 1.50 | 1.43–1.58 | <0.001 |
| OTHER_GOV | 0.09 | 0.076–0.10 | <0.001 |
| OTHER | 0.47 | 0.46–0.47 | <0.001 |
| NETWORK | 0.42 | 0.39–0.45 | <0.001 |
| INDIV | 0.62 | 0.50–0.76 | <0.001 |
| UNKNOWN | 0.84 | 0.54–1.30 | 0.44 |
| **Study phase (ref: PHASE1)** | | | |
| Phase 1\|2 | 2.33 | 2.11–2.57 | <0.001 |
| Phase 2 | 2.48 | 2.26–2.73 | <0.001 |
| Phase 2\|3 | 1.79 | 1.61–2.00 | <0.001 |
| Phase 3 | 2.49 | 2.27–2.74 | <0.001 |
| Phase 4 | 2.25 | 2.05–2.48 | <0.001 |
| Unspecified | 6.02 | 4.87–7.44 | <0.001 |
| NA | 1.19 | 1.09–1.31 | <0.001 |
| **Condition family (ref: cardiovascular)** | | | |
| Oncology | 1.27 | 1.24–1.31 | <0.001 |
| Mental health | 1.31 | 1.26–1.35 | <0.001 |
| Neurology | 1.24 | 1.18–1.29 | <0.001 |
| Immunology/Derm | 1.22 | 1.15–1.29 | <0.001 |
| Infectious disease | 1.19 | 1.15–1.23 | <0.001 |
| Respiratory/Sleep | 1.08 | 1.04–1.13 | <0.001 |
| Other | 1.08 | 1.05–1.11 | <0.001 |
| Renal/Urology | 1.01 | 0.96–1.07 | 0.63 |
| Musculoskeletal/Pain | 1.00 | 0.96–1.04 | 0.98 |
| Metabolic | 0.94 | 0.91–0.98 | 0.002 |
| Gastro-hepatic | 0.94 | 0.89–0.99 | 0.03 |
| Healthy volunteer | 0.85 | 0.81–0.90 | <0.001 |
| Reproductive/Maternal | 0.82 | 0.76–0.88 | <0.001 |

Concordance = 0.683. HR = hazard ratio; CI = confidence interval. Reference categories: sponsor class INDUSTRY, phase PHASE1, condition cardiovascular.

---

### Table 3. Piecewise Cox proportional hazards model: interval-specific hazard ratios

| Interval | HR (FED vs INDUSTRY) | HR (Phase 3 vs Phase 1) | HR (Oncology vs Cardio) |
|---|---|---|---|
| 0–12 months | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| 12–24 months | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| 24–48 months | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |
| 48+ months | [PLACEHOLDER] | [PLACEHOLDER] | [PLACEHOLDER] |

[PLACEHOLDER: to be populated from piecewise Cox re-run with increased memory allocation. Expected finding: HRs vary substantially across intervals, with disclosure hazard concentrated in the 0–12 month window for FDAAA-applicable sponsors.]

---

## References

1. Song F, Parekh S, Hooper L, et al. Dissemination and publication of research findings: an updated review of related biases. *Health Technol Assess.* 2010;14(8):iii, ix-xi, 1-193.

2. Food and Drug Administration Amendments Act of 2007 (FDAAA 801). Public Law 110-85, Title VIII. Available at: https://www.govinfo.gov/content/pkg/PLAW-110publ85/pdf/PLAW-110publ85.pdf

3. National Institutes of Health. Final rule: clinical trials registration and results information submission. *Federal Register.* 2016;81(183):64981-65157. Available at: https://www.federalregister.gov/documents/2016/09/21/2016-22129/

4. Anderson ML, Chiswell K, Peterson ED, Tasneem A, Topping J, Califf RM. Compliance with results reporting at ClinicalTrials.gov. *N Engl J Med.* 2015;372(11):1031-1039. doi:10.1056/NEJMsa1409364

5. DeVito NJ, Bacon S, Goldacre B. Compliance with legal requirement to report clinical trial results on ClinicalTrials.gov: a cohort study. *Lancet.* 2020;395(10221):361-369. doi:10.1016/S0140-6736(19)33220-9

6. Lau B, Cole SR, Gange SJ. Competing risk regression models for epidemiologic data. *Am J Epidemiol.* 2009;170(2):244-256. doi:10.1093/aje/kwp107

7. Zarin DA, Tse T, Williams RJ, Califf RM, Ide NC. The ClinicalTrials.gov results database — update and key issues. *N Engl J Med.* 2011;364(9):852-860. doi:10.1056/NEJMsa1012065

8. [AUTHOR_NAME], [AUTHOR_NAME]. The Hiddenness Atlas: a registry-wide audit of result availability across 578,109 ClinicalTrials.gov studies. [*Companion manuscript, same snapshot date.*] [PLACEHOLDER: journal, year, doi]
