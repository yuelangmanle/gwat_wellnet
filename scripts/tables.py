#!/usr/bin/env python3
"""Generate Table 1 (manuscript/tables/table1_strata.csv and .md).

Columns: Stratum | N | Falling (%) | Median well slope | Reference slope |
95% CI | SNR | n_min (sim) | n_min (law) | n_min/N. All slopes in m/yr.
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))
OUT = os.path.abspath(os.path.join(HERE, "..", "manuscript", "tables"))
os.makedirs(OUT, exist_ok=True)
FT = 0.3048

AQ_NAMES = json.load(open(os.path.join(DATA, "aquifer_names.json")))
AGG_DEPENDENT = {"S400FLORDN", "S100RIOGRD", "S100CNRLVL", "N600CMBPLV"}
AQ_NAMES["UNKNOWN"] = "Unassigned wells"
AQ_NAMES["CONUS"] = "CONUS (national network)"


def clean_zero(val, tol=5e-4):
    if pd.isna(val):
        return val
    return 0.0 if abs(float(val)) < tol else float(val)


def build_row(aq, sub, res, th, boot, law):
    rr = res[(res.stratum == aq) & (~res.clustered)]
    ref = clean_zero(float(rr.ref_slope.iloc[0]) * FT) if len(rr) else np.nan
    b = boot.loc[aq] if aq in boot.index else None
    ci = "[%.3f, %.3f]" % (clean_zero(b.lo * FT), clean_zero(b.hi * FT)) if b is not None else "-"
    t = th[(th.stratum == aq) & (~th.clustered)]
    nmin = t.n_min.iloc[0] if len(t) and pd.notna(t.n_min.iloc[0]) else np.nan
    n_all = len(sub)
    dag = aq in AGG_DEPENDENT
    lr = law.loc[aq] if aq in law.index else None
    med_well = clean_zero(sub.slope_ft_yr.median() * FT)
    return dict(
        Stratum=AQ_NAMES.get(aq, aq) + ("\u2020" if dag else ""),
        N=str(n_all),
        Falling_pct="%.1f" % (100 * (sub.slope_ft_yr > 0).mean()),
        Median_well_slope="%.3f" % med_well,
        Reference_slope="%.3f" % ref,
        CI_95=ci,
        SNR=("%.3f" % clean_zero(float(lr.snr)) if lr is not None else "-"),
        n_min=(str(int(nmin)) if pd.notna(nmin) else "-"),
        n_min_law=(str(int(round(float(lr.nmin_pred)))) if lr is not None else "-"),
        n_min_over_N=("-" if pd.isna(nmin) else "%.3f" % (nmin / n_all)),
    )


def main():
    trends = pd.read_csv(os.path.join(DATA, "well_trends.csv"),
                         dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"),
                        dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    wells = wells.join(trends, how="inner")
    res = pd.read_csv(os.path.join(DATA, "experiment_results.csv"))
    th = pd.read_csv(os.path.join(DATA, "thresholds.csv"))
    boot = pd.read_csv(os.path.join(DATA, "robustness_bootstrap.csv")).set_index("stratum")
    law = pd.read_csv(os.path.join(DATA, "law_nmin.csv")).set_index("stratum")

    wells["aq"] = wells.national_aquifer_code.fillna("UNKNOWN")
    counts = wells.aq.value_counts()
    strata = [s for s in counts[counts >= 100].index]

    rows = [build_row(aq, wells[wells.aq == aq], res, th, boot, law) for aq in strata]
    rows.sort(key=lambda r: int(r["N"]), reverse=True)
    rows.append(build_row("CONUS", wells, res, th, boot, law))

    t1 = pd.DataFrame(rows)
    t1.to_csv(os.path.join(OUT, "table1_strata.csv"), index=False)

    md_headers = [
        "Stratum", "N", "Falling (%)", "Median well slope",
        "Reference slope", "95% CI", "SNR", "n_min (sim)", "n_min (law)", "n_min/N"
    ]
    md_lines = [
        "| " + " | ".join(md_headers) + " |",
        "| :--- | ---: | ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: |",
    ]
    for r in rows:
        vals = [
            r["Stratum"], r["N"], r["Falling_pct"], r["Median_well_slope"],
            r["Reference_slope"], r["CI_95"], r["SNR"], r["n_min"],
            r["n_min_law"], r["n_min_over_N"]
        ]
        md_lines.append("| " + " | ".join(vals) + " |")

    caption = (
        "# Table 1. Eligible monitoring wells, reference groundwater-level trends "
        "(1975\u20132025; m yr\u207b\u00b9, where positive values indicate falling water levels), "
        "well-bootstrap 95% confidence intervals (CIs) of the reference slope, "
        "signal-to-noise ratios (SNR), and simulated versus predicted minimum reliable "
        "network sizes (\u2264 10% sign-reversal risk under uniform random removal): "
        "simulated n_min and closed-form prediction n_min(law) = Nz\u00b2/(z\u00b2+(c\u00b7SNR)\u00b2N). "
        "A dash (-) indicates that the stratum does not reach \u2264 10% sign-reversal risk on "
        "the standard retention grid (up to 90% retention); \u2020 indicates strata where the "
        "reference trend sign changes under spatial-median or inverse-density aggregation.\n\n"
    )
    with open(os.path.join(OUT, "table1_strata.md"), "w") as fh:
        fh.write(caption + "\n".join(md_lines) + "\n")
    print("merged table ->", os.path.join(OUT, "table1_strata.csv"))


if __name__ == "__main__":
    main()
