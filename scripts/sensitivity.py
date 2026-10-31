#!/usr/bin/env python3
"""Aggregation sensitivity.

The reference network trend uses the mean annual anomaly across wells. This
script recomputes every stratum's reference slope with the spatial median and
reports sign agreement between the two aggregation statistics — the
robustness statement for the Methods' aggregation choice.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from experiment import build_matrix, DATA  # noqa: E402
from trends import yue_pilon_mk  # noqa: E402


def network_stat(mat, stat):
    series = stat(mat, axis=0)
    ok = ~np.isnan(series)
    if ok.sum() < 20:
        return np.nan
    slope, _, _, _, _ = yue_pilon_mk(np.arange(len(series))[ok], series[ok])
    return slope


def main():
    trends = pd.read_csv(os.path.join(DATA, "well_trends.csv"),
                         dtype={"monitoring_location_id": str})
    wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"),
                        dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    anoms = pd.read_csv(os.path.join(DATA, "annual_anoms.csv"),
                        dtype={"monitoring_location_id": str})
    wells = wells.join(trends.set_index("monitoring_location_id"),
                       how="inner").dropna(subset=["slope_ft_yr"])
    anoms = anoms[anoms.monitoring_location_id.isin(wells.index)]
    years = np.sort(anoms.year.unique())
    years = years[(years >= 1975) & (years <= 2025)]

    aquifer = wells.national_aquifer_code.fillna("UNKNOWN")
    counts = aquifer.value_counts()
    strata = [("CONUS", wells.index.values)] + [
        (aq, wells.index[aquifer == aq].values) for aq in counts[counts >= 100].index]

    rows = []
    for name, sel in strata:
        mat = build_matrix(anoms, sel, years)
        s_mean = network_stat(mat, np.nanmean)
        s_med = network_stat(mat, np.nanmedian)
        rows.append(dict(stratum=name, n=len(sel),
                         slope_mean=s_mean, slope_median=s_med,
                         sign_agree=bool(np.sign(s_mean) == np.sign(s_med))))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DATA, "sensitivity_aggregation.csv"), index=False)
    print(df.to_string(max_rows=25))
    print("sign agreement: %d/%d strata" % (df.sign_agree.sum(), len(df)))


if __name__ == "__main__":
    main()
