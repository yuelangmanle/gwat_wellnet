#!/usr/bin/env python3
"""Per-well water-level trend estimation.

Method (documented for the manuscript):
- Deseasonalize each well's discrete measurements by subtracting
  month-of-year medians (months with >= 2 observations); months with fewer
  observations are centred on the well's overall median.
- Aggregate residuals to annual medians; require >= 20 annual values.
- Theil-Sen slope on annual medians (robust to outliers/non-normality).
- Mann-Kendall test with Yue-Pilon trend-free pre-whitening to correct the
  variance for lag-1 serial correlation (AR removed only if significant at
  the 95% level).

Outputs: data/well_trends.csv (one row per eligible well).
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))
MIN_ANNUAL = 20


def deseasonalize_annual(df):
    """Return annual median residuals (Series indexed by year)."""
    med_all = df["value"].median()
    month_med = df.groupby("month")["value"].median()
    counts = df.groupby("month")["value"].size()

    def corr(row):
        m = row["month"]
        if counts.get(m, 0) >= 2:
            return row["value"] - month_med[m]
        return row["value"] - med_all

    df = df.copy()
    df["resid"] = df.apply(corr, axis=1)
    ann = df.groupby("year")["resid"].median()
    return ann


def mk_s_var(s):
    """Mann-Kendall S statistic and variance with tie correction."""
    n = len(s)
    s_stat = 0
    _, tie_counts = np.unique(s, return_counts=True)
    tie_term = np.sum(tie_counts * (tie_counts - 1) * (2 * tie_counts + 5))
    for k in range(n - 1):
        s_stat += np.sum(s[k + 1:] > s[k]) - np.sum(s[k + 1:] < s[k])
    var_s = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0
    return s_stat, var_s


def mk_p(s_series):
    """Two-sided MK p-value from S and var(S)."""
    s_stat, var_s = mk_s_var(s_series)
    if var_s <= 0:
        return np.nan
    if s_stat > 0:
        z = (s_stat - 1) / np.sqrt(var_s)
    elif s_stat < 0:
        z = (s_stat + 1) / np.sqrt(var_s)
    else:
        z = 0.0
    return 2.0 * (1.0 - stats.norm.cdf(abs(z)))


def yue_pilon_mk(years, values):
    """Theil-Sen slope + MK p-value with Yue-Pilon pre-whitening.

    Returns (slope, p, n, r1, ar_corrected).
    """
    years = np.asarray(years, dtype=float)
    values = np.asarray(values, dtype=float)
    n = len(values)

    res = stats.theilslopes(values, years)
    beta = res[0]

    detrended = values - beta * years
    # lag-1 autocorrelation (biased estimator), mean-centred
    dy = detrended - detrended.mean()
    r1 = np.sum(dy[:-1] * dy[1:]) / max(np.sum(dy ** 2), 1e-12)

    ar_sig = abs(r1) > (1.96 / np.sqrt(n))
    if ar_sig:
        pw = np.empty(n)
        pw[0] = detrended[0]
        pw[1:] = detrended[1:] - r1 * detrended[:-1]
        values_test = pw + beta * years
    else:
        values_test = values

    p = mk_p(values_test)
    return float(beta), float(p), n, float(r1), bool(ar_sig)


def main():
    meas = pd.read_csv(os.path.join(DATA, "well_measurements.csv"),
                       usecols=["monitoring_location_id", "value", "year", "month"])
    if meas is None or len(meas) == 0:
        raise SystemExit("run build_dataset.py first")

    rows = []
    anom_frames = []
    for wid, g in meas.groupby("monitoring_location_id", sort=False):
        ann = deseasonalize_annual(g)
        if len(ann) < MIN_ANNUAL:
            continue
        slope, p, n, r1, arc = yue_pilon_mk(ann.index.values, ann.values)
        rows.append(dict(monitoring_location_id=wid, slope_ft_yr=slope, p_value=p,
                         n_annual=n, lag1=r1, ar_corrected=arc,
                         mean_level_ft=g["value"].mean()))
        anom_frames.append(pd.DataFrame(
            {"monitoring_location_id": wid, "year": ann.index.values,
             "resid": ann.values}))

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA, "well_trends.csv"), index=False)
    pd.concat(anom_frames, ignore_index=True).to_csv(
        os.path.join(DATA, "annual_anoms.csv"), index=False)
    print("wells with trends:", len(out))
    print(out["slope_ft_yr"].describe())
    print("declining (slope>0, depth increasing):",
          float((out["slope_ft_yr"] > 0).mean()))


if __name__ == "__main__":
    main()
