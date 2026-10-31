#!/bin/bash
# Full pipeline: raw download -> dataset -> trends -> experiment -> law -> tests -> figures
# Order matters; every stage writes its products under data/ and figures/.
set -e
cd "$(dirname "$0")/.."
echo "== stage 0: download (optional; ~9.5M rows, run once) =="
# python3 scripts/download_ogc.py --collection field-measurements --params "parameter_code=72019&datetime=1950-01-01/1989-12-31" --name gw_a
# (repeat for decades 1990s/2000s/2010s-2020s, and monitoring-locations as 'sites')
echo "== stage 1: build dataset =="
python3 scripts/build_dataset.py
echo "== stage 2: per-well trends =="
python3 scripts/trends.py
echo "== stage 3: thinning experiment (1,000 replicates) =="
python3 scripts/experiment.py --replicates 1000
echo "== stage 3b: retention strategies =="
python3 scripts/strategies.py
echo "== stage 3c: aggregation sensitivity =="
python3 scripts/sensitivity.py
echo "== stage 3d: robustness (density/bootstrap/LOO) =="
python3 scripts/robustness.py
echo "== stage 3e: sampling law + n_min + clustered n_eff =="
python3 scripts/law.py
echo "== stage 3f: sigloss companion-form scope =="
python3 scripts/law_scope.py
echo "== stage 3g: direct tests of the predictions =="
python3 scripts/censored_test.py
python3 scripts/decomposition.py
echo "== stage 3h: recipe and coverage tests =="
python3 scripts/snr_recipes.py
python3 scripts/coverage_tests.py
python3 scripts/sim_validation.py
echo "== stage 3i: scope of the archived experiment =="
python3 scripts/scope_tests.py
python3 scripts/external_validation.py
python3 scripts/frequency_axis.py
echo "== stage 4: tables + figures =="
python3 scripts/tables.py
python3 scripts/figures.py
python3 scripts/graphical_abstract.py
echo "== stage 5: audit =="
python3 scripts/audit_numbers.py
echo "pipeline complete"
