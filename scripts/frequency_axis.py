#!/usr/bin/env python3
"""B: the frequency axis. Restrict every well to ONE measurement per year
(its most common month), recompute well slopes on that reduced schedule,
then stratum SNR and the law's n_min. Compares annual-visit monitoring with
the full archive schedule.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

DATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
Z = float(norm.ppf(0.9))
C = json.load(open(os.path.join(DATA, "law_fit.json")))["c"]


def main():
    meas = pd.read_csv(os.path.join(DATA, "well_measurements.csv"),
                       usecols=["monitoring_location_id", "value", "year", "month", "date"])
    # one visit per well-year: the modal month's median value
    meas["date"] = pd.to_datetime(meas["date"])
    per_wy = meas.groupby(["monitoring_location_id", "year"]).agg(
        month=("month", lambda s: s.mode().iat[0]),
        n_days=("date", "nunique")).reset_index()
    # keep only rows in the modal month
    m = meas.merge(per_wy[["monitoring_location_id", "year", "month"]],
                   on=["monitoring_location_id", "year", "month"], how="inner")
    # one value per well-year (median of that month's days)
    annual = m.groupby(["monitoring_location_id", "year"], as_index=False).agg(
        value=("value", "median"))
    counts = annual.groupby("monitoring_location_id").year.nunique()
    keep = counts[counts >= 20].index
    annual = annual[annual.monitoring_location_id.isin(keep)]
    print(f"wells with >=20 annual visits after thinning: {len(keep):,}")

    def slope(d):
        if len(d) < 5:
            return np.nan
        y = d.value.values
        if y.std() == 0:
            return 0.0
        return float(np.polyfit(d.year, y, 1)[0])
    slopes = annual.groupby("monitoring_location_id").apply(slope)
    sl = pd.DataFrame({"monitoring_location_id": slopes.index, "slope_1v": slopes.values})

    wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"),
                        dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    trends = pd.read_csv(os.path.join(DATA, "well_trends.csv"),
                         dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    base = wells.join(trends, how="inner")
    red = base.join(sl.set_index("monitoring_location_id"), how="inner")
    print(f"matched to eligible wells: {len(red):,} (baseline {len(base):,})")

    out = []
    for aq, g in red.groupby(red.national_aquifer_code.fillna("UNKNOWN")):
        if len(g) < 100:
            continue
        sd1 = g.slope_1v.std(ddof=0)
        sd0 = g.slope_ft_yr.std(ddof=0)
        mu1 = g.slope_1v.mean()
        N = len(g)
        # reference slope: mean anomaly hydrograph unchanged in sign; use
        # SNR = |mean of well slopes| / SD (operator formula) for both schedules
        snr1 = abs(mu1) / sd1 if sd1 > 0 else np.nan
        snr0 = abs(g.slope_ft_yr.mean()) / sd0 if sd0 > 0 else np.nan
        nmin = lambda s: (Z / (C * s)) ** 2 if np.isfinite(s) and s > 0 else np.nan
        out.append(dict(stratum=aq, N=N, snr_full=snr0, snr_1visit=snr1,
                         nmin_full=nmin(snr0), nmin_1visit=nmin(snr1)))
    df = pd.DataFrame(out)
    df["snr_ratio"] = df.snr_1visit / df.snr_full
    df["nmin_ratio"] = df.nmin_1visit / df.nmin_full
    df.to_csv(os.path.join(DATA, "frequency_axis.csv"), index=False)
    strong = df[df.nmin_full < 0.5 * df.N]  # size-independent regime
    print("\nSNR (1 visit/yr vs full): median ratio %.3f" % df.snr_ratio.median())
    print("n_min ratio (1 visit vs full): median %.2f across strata; %.2f in strong-signal strata"
          % (df.nmin_ratio.median(), strong.nmin_ratio.median()))
    print(df[["stratum", "N", "snr_full", "snr_1visit", "nmin_ratio"]].to_string(index=False))
    json.dump(dict(snr_ratio_median=float(df.snr_ratio.median()),
                   nmin_ratio_median=float(df.nmin_ratio.median()),
                   nmin_ratio_strong=float(strong.nmin_ratio.median())),
              open(os.path.join(DATA, "frequency_axis_summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
