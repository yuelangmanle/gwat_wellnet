#!/usr/bin/env python3
"""Direct test of the censored predictions (strata whose grids stop above
their predicted minimum).

Runs random-mode thinning at retention 0.92/0.95/0.98 (and 0.99 for the
largest) for the three never-stable strata plus CONUS clustered at
0.40/0.50/0.70 (law-n_min side of the clustered threshold), then compares
observed flip vs the law's prediction with c = 1.27.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from experiment import build_matrix, network_slope, block_labels, DATA  # noqa: E402
from trends import yue_pilon_mk  # noqa: E402

RNG = np.random.default_rng(20261024)
C = json.load(open(os.path.join(DATA, "law_fit.json")))["c"]
Z = float(norm.ppf(0.9))
YEARS = None


def law_flip(n, N, snr):
    x = C * snr * np.sqrt(n * N / (N - n))
    return float(norm.cdf(-x))


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

    # SD of well slopes (the law's S) per stratum
    sd = wells.groupby(wells.national_aquifer_code.fillna("UNKNOWN")).slope_ft_yr.std()

    print("== censored strata: random-mode flip at high retention (1,000 reps) ==")
    print(f"{'stratum':12s} {'N':>6s} {'n':>6s} {'obs_flip':>8s} {'law':>6s} {'binom95CI':>16s}")
    rows = []
    for aq, rets in [("N100BSNRGB", [0.92, 0.95, 0.98]),
                     ("S100SURFCL", [0.92, 0.95]),
                     ("S400FLORDN", [0.92, 0.95, 0.98])]:
        sel = wells.index[wells.national_aquifer_code.fillna("UNKNOWN") == aq]
        mat = build_matrix(anoms, sel, years)
        N = len(sel)
        ref = network_slope(mat, np.arange(N))[0]
        snr = abs(ref) / sd[aq]
        sign = 1 if ref > 0 else -1
        for f in rets:
            n = int(round(f * N))
            flips = 0
            for _ in range(1000):
                keep = RNG.choice(N, size=n, replace=False)
                s = network_slope(mat, keep)[0]
                if np.isfinite(s) and (1 if s > 0 else -1) != sign:
                    flips += 1
            p = flips / 1000
            from scipy.stats import beta as _beta
            lo = _beta.ppf(0.025, flips, 1000 - flips + 1) if flips > 0 else 0.0
            hi = _beta.ppf(0.975, flips + 1, 1000 - flips) if flips < 1000 else 1.0
            pred = law_flip(n, N, snr)
            print(f"{aq:12s} {N:6d} {n:6d} {p:8.3f} {pred:6.3f} [{lo:.3f},{hi:.3f}]")
            rows.append(dict(stratum=aq, N=N, n=n, obs=p, lo=lo, hi=hi, law=pred,
                             law_min_ret=1 - min(n / N, 1)))
    pd.DataFrame(rows).to_csv(os.path.join(DATA, "censored_test.csv"), index=False)

    print("\n== CONUS clustered above 30%: does flip fall as retention grows? ==")
    sel = wells.index
    meta = wells.loc[sel]
    mat = build_matrix(anoms, sel, years)
    N = len(sel)
    ref = network_slope(mat, np.arange(N))[0]
    snr_conus = abs(ref) / wells.slope_ft_yr.std()
    sign = 1 if ref > 0 else -1
    labels = block_labels(meta)
    uniq, cnt = np.unique(labels, return_counts=True)
    keepb = set(uniq[cnt >= 10])
    mask = np.array([l in keepb for l in labels])
    mat_c, labels_c = mat[mask], labels[mask]
    base_n = int(mask.sum())
    blocks = {b: np.where(labels_c == b)[0] for b in np.unique(labels_c)}
    keys = list(blocks)
    for f in [0.40, 0.50, 0.70, 0.90]:
        n_keep = int(f * base_n)
        flips = done = 0
        for _ in range(500):
            order = RNG.permutation(len(keys))
            selr, taken = [], 0
            for oi in order:
                idx = blocks[keys[oi]]
                if taken + len(idx) <= max(n_keep * 1.15, n_keep + 10):
                    selr.append(idx)
                    taken += len(idx)
                if taken >= n_keep * 0.9:
                    break
            keep_rows = np.concatenate(selr) if selr else np.arange(base_n)
            s = network_slope(mat_c, keep_rows)[0]
            if np.isnan(s):
                continue
            done += 1
            if (1 if s > 0 else -1) != sign:
                flips += 1
        print(f"  clustered f={f:.2f} (~{taken} wells): flip={flips/max(done,1):.3f}")


if __name__ == "__main__":
    main()
