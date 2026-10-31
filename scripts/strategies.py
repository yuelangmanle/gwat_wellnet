#!/usr/bin/env python3
"""Deterministic well-retention strategies vs random removal.

For each stratum, three ways of choosing which wells to keep at a given
network size:
- random (from experiment.py; expectation + MC spread)
- longest-record-first: rank wells by number of annual values (n_annual)
- spatially-spread-first: farthest-point ordering of well coordinates;
  keep the first n wells of the ordering

The two deterministic strategies give a single flip-rate curve each
(sequential deletion = nested subsets); random is the reference band.
Outputs: data/strategy_results.csv
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from experiment import build_matrix, network_slope, DATA  # noqa: E402
from trends import yue_pilon_mk  # noqa: E402

# fine nested-deletion grid for the step curves (strategies.py only)
RETENTION = [round(1.0 - 0.02 * i, 2) for i in range(0, 48)]  # 1.00 -> 0.06


def farthest_point_order(xy):
    """Order of well indices produced by greedy farthest-point sampling."""
    n = len(xy)
    order = np.empty(n, dtype=int)
    selected = np.zeros(n, dtype=bool)
    order[0] = 0
    selected[0] = True
    dist = np.linalg.norm(xy - xy[0], axis=1)
    for i in range(1, n):
        nxt = int(np.argmax(dist))
        order[i] = nxt
        selected[nxt] = True
        d = np.linalg.norm(xy - xy[nxt], axis=1)
        dist = np.minimum(dist, d)
        dist[selected] = -1.0
    return order


def nested_flip_curve(mat, ordered_idx, ref_slope, ref_sign, ref_sig, n_all):
    """Flip indicators along the nested sequence keep=first k wells."""
    out = []
    for f in RETENTION:
        k = max(10, int(np.floor(f * n_all)))
        slope, p = network_slope(mat, ordered_idx[:k])
        if np.isnan(slope):
            continue
        out.append(dict(
            retention=f, n_keep=k,
            flip=int((1 if slope > 0 else -1) != ref_sign),
            sig_loss=int(ref_sig and (p >= 0.05 or (1 if slope > 0 else -1) != ref_sign)),
            slope=slope, abs_bias=abs(slope - ref_slope)))
    return out


def smooth_curve(vals, window=5):
    """Centred rolling mean for a monotone nested-subset curve."""
    s = pd.Series(vals).rolling(window, center=True, min_periods=1).mean()
    return s.values


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

    aquifer = wells["national_aquifer_code"].fillna("UNKNOWN")
    counts = aquifer.value_counts()
    big = counts[counts >= 100].index.tolist()
    strata = [("CONUS", wells.index.values)] + [
        (aq, wells.index[aquifer == aq].values) for aq in big]

    years = np.sort(anoms.year.unique())
    years = years[(years >= 1975) & (years <= 2025)]

    rows = []
    for name, sel in strata:
        mat = build_matrix(anoms, sel, years)
        meta = wells.loc[sel]
        ref_slope, ref_p = network_slope(mat, np.arange(len(sel)))
        if np.isnan(ref_slope):
            continue
        ref_sign = 1 if ref_slope > 0 else -1
        ref_sig = ref_p < 0.05
        n_all = len(sel)

        # deterministic strategy 1: longest record first
        by_len = meta.sort_values("n_annual", ascending=False).index.values
        pos = {w: i for i, w in enumerate(sel)}
        idx_len = np.array([pos[w] for w in by_len])
        for rec in nested_flip_curve(mat, idx_len, ref_slope, ref_sign, ref_sig, n_all):
            rec.update(stratum=name, strategy="longest_first", n_all=n_all,
                       ref_slope=ref_slope, ref_p=ref_p)
            rows.append(rec)

        # deterministic strategy 2: farthest-point spatial spread
        fp_order = farthest_point_order(meta[["x", "y"]].values)
        for rec in nested_flip_curve(mat, fp_order, ref_slope, ref_sign, ref_sig, n_all):
            rec.update(stratum=name, strategy="spatial_spread", n_all=n_all,
                       ref_slope=ref_slope, ref_p=ref_p)
            rows.append(rec)

        # nested removal along random orders (reference band, 50 orders)
        import numpy as _np
        curves = []
        for rep in range(50):
            perm = _np.random.default_rng(1000 + rep).permutation(n_all)
            curves.append([(rec["retention"], rec["n_keep"], rec["flip"], rec["abs_bias"])
                           for rec in nested_flip_curve(mat, perm, ref_slope, ref_sign,
                                                        ref_sig, n_all)])
        # aggregate across orders at matched retention
        from collections import defaultdict
        by_ret = defaultdict(list)
        for cur in curves:
            for ret, nk, flip, bias in cur:
                by_ret[ret].append((nk, flip, bias))
        for ret, vals in sorted(by_ret.items(), key=lambda kv: -kv[0]):
            arr = _np.array(vals)
            rows.append(dict(stratum=name, strategy="nested_random",
                             retention=ret, n_keep=float(_np.mean(arr[:, 0])),
                             flip=float(_np.mean(arr[:, 1])),
                             sig_loss=0, slope=_np.nan, abs_bias=float(_np.mean(arr[:, 2])),
                             stratum_n=len(vals), n_all=n_all,
                             ref_slope=ref_slope, ref_p=ref_p))

    pd.DataFrame(rows).to_csv(os.path.join(DATA, "strategy_results.csv"), index=False)
    print("strategy rows:", len(rows), "->", os.path.join(DATA, "strategy_results.csv"))


if __name__ == "__main__":
    main()
