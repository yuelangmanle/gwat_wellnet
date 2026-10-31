# gwat_wellnet — Monitoring-network density and the reliability of groundwater-level trends

Research pipeline behind the *Groundwater* (Wiley/NGWA) paper
**"How sparse is too sparse? A closed-form law for the reliability of
groundwater-level trends."**

**Question.** How reliably does a reduced groundwater monitoring network reproduce
the trend conclusion of the full network, and what is the minimum network size
that keeps the conclusion stable?

**Design in one paragraph.** All discrete USGS depth-to-water measurements
(parameter 72019) are pulled from the USGS Water Data API, quality-screened, and
reduced to per-well annual medians for wells with >= 30 measurements spanning
>= 20 calendar years. Well trends are Theil-Sen slopes with Yue-Pilon
pre-whitened Mann-Kendall significance. For each principal aquifer (and the
CONUS-wide network) the reference trend is the slope of the averaged anomaly
hydrograph; networks are then virtually thinned (uniform random, clustered by
0.5-degree blocks, longest-record-first, spatial-spread-first) and the
probability of sign flips, significance loss, and slope bias is mapped against
retained network size, yielding per-aquifer minimum reliable network sizes
n_min = Nz²/(z² + (c·SNR)²N) with c = 1.27 (95% CI 1.04–1.55).

## Pipeline

**Data.**
1. `scripts/download_ogc.py` — bulk cursor-paged download from
   `api.waterdata.usgs.gov/ogcapi/v1` (field-measurements, parameter_code=72019,
   sharded by decade; monitoring-locations for metadata). Checkpointed, resumable.
2. `scripts/build_dataset.py` — QA filters (range, approval, pumping flags,
   CONUS bbox), well-day medians, site join, eligibility screen.
3. `scripts/trends.py` — deseasonalize (month-of-year medians), annual medians,
   Theil-Sen + Yue-Pilon MK.

**Experiment.**

4. `scripts/experiment.py` — deletion experiment (uniform random + clustered),
   1,000 replicates per (stratum, retention).
5. `scripts/strategies.py` — deterministic retention strategies
   (longest-record-first, farthest-point spatial spread).
6. `scripts/sensitivity.py` — aggregation sensitivity (mean vs spatial median).

**The law.**

7. `scripts/law.py` — the finite-population sampling law, its calibration,
   closed-form n_min, and the clustered-loss effective size.
8. `scripts/law_scope.py` — scope test on significance loss (companion form).

**Tests.**

9. `scripts/robustness.py` — density weighting, bootstrap CIs, leave-one-out.
10. `scripts/censored_test.py` — direct thinning test of the censored
    (never-reached) predictions.
11. `scripts/decomposition.py` — benchmark (subset mean) versus pipeline
    estimator, per stratum.
12. `scripts/snr_recipes.py` — matched-SNR decomposition, OLS/mean-of-slopes
    recipe tests, frequency axis under the paper's own recipe.
13. `scripts/coverage_tests.py` — binomial coverage by flip regime, t-link
    competitor, balanced-panel composition test, Moran's I.
14. `scripts/sim_validation.py` — synthetic-population testbed (SI S1) and the
    forward Kish check of the clustered arm (SI S2).
15. `scripts/scope_tests.py` — significance-loss reporting, tolerance grid,
    window sensitivity, confined/unconfined split, per-stratum clustered
    deletion, block-size sensitivity.
16. `scripts/external_validation.py` — hold-out validation of c; OLS
    estimator robustness.
17. `scripts/frequency_axis.py` — one-measurement-per-well-year axis.

**Outputs.**

18. `scripts/tables.py` — Table 1 (`manuscript/tables/`).
19. `scripts/figures.py` — Figs. 1–3 and the SI figure
    (`figures/figS_strategies`), 17.15 cm double-column exports
    (SVG/PDF/PNG/TIFF).
20. `scripts/graphical_abstract.py` — graphical abstract.
21. `scripts/audit_numbers.py` — recomputes every number quoted in the
    manuscript from the archived data files (90 checks).

One command: `bash scripts/run_all.sh` (stage 0 download optional).

## Repository contents

| Path | Contents |
|---|---|
| `scripts/` | the pipeline above (`run_all.sh` runs it end to end) |
| `data/` | analysis-layer products: per-well trends (29,717 wells), annual anomalies, the thinning experiment, strategy runs, law fit and coverage tests, robustness and validation files |
| `figures/` | Figs. 1–3 + `figS_strategies` (SI), graphical abstract, `captions.md` |
| `manuscript/manuscript_draft.md` | manuscript source text (the audit checks quoted numbers against it) |
| `manuscript/tables/` | Table 1 (markdown + CSV) |
| `manuscript/SI_sim_validation.md` | Supporting Information source (Sections S1–S5: synthetic testbed, forward cluster check, significance-loss companion model, retention strategies, estimator decomposition and sensitivity checks) |

The raw measurement snapshot (~9.5 million rows, USGS Water Data API,
October 2026) is not included; `scripts/download_ogc.py` reproduces it from the
public API, and the download is deterministic given the API state.

## Requirements

Python 3 with numpy, pandas, scipy, and matplotlib. Figures additionally
use Arial or Helvetica.

## Data and code sources

Measurements: U.S. Geological Survey National Water Information System
(https://doi.org/10.5066/F7P55KJN), public domain. Code: MIT license (see
`LICENSE`). Analysis code by Zijie Deng; questions: 3335196397vip@gmail.com.
