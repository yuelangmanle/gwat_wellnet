#!/usr/bin/env python3
"""Network-thinning (deletion) experiment.

Question: how reliably does a reduced monitoring network reproduce the trend
conclusion of the full network, and what is the minimum network size that
keeps the conclusion stable?

Design (documented for the manuscript):
- Strata: national aquifer codes with >= MIN_STRAT_WELLS eligible wells,
  plus the CONUS-wide network.
- Reference network trend: mean annual anomaly hydrograph across wells
  (depth-to-water residuals, ft), Theil-Sen slope + Yue-Pilon MK.
- For retention fraction f, sample floor(f*N) wells without replacement
  REPLICATES times, rebuild the network hydrograph over the SAME years,
  recompute slope/p.
- Metrics per (stratum, f): P(sign flip), P(loss of significance given the
  reference is significant), mean absolute slope bias.
- Two removal modes: uniform random (monitoring design) and clustered
  removal by 0.5-degree spatial blocks (budget-driven contraction).
- Outputs: data/experiment_results.csv, data/thresholds.csv.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trends import yue_pilon_mk  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))

MIN_STRAT_WELLS = 100
RETENTION = [0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01]
REPLICATES = 1000
RNG = np.random.default_rng(20261023)


def build_matrix(anoms, wells_idx, years):
    """Wide well x year matrix of residuals (NaN where absent)."""
    sub = anoms[anoms.monitoring_location_id.isin(wells_idx)]
    piv = sub.pivot_table(index="monitoring_location_id", columns="year",
                          values="resid", aggfunc="mean")
    piv = piv.reindex(index=wells_idx, columns=years)
    return piv.values  # rows aligned with wells_idx


def network_slope(mat, rows_idx, min_years=20):
    """nanmean over selected well rows, then TS slope + MK p."""
    if len(rows_idx) == 0:
        return np.nan, np.nan
    series = np.nanmean(mat[rows_idx], axis=0)
    ok = ~np.isnan(series)
    if ok.sum() < min_years:
        return np.nan, np.nan
    slope, p, _, _, _ = yue_pilon_mk(np.arange(len(series))[ok], series[ok])
    return slope, p


def block_labels(wells_meta, cell=0.5):
    """Deterministic 0.5-degree spatial block labels (no dependencies)."""
    return (np.floor(wells_meta["y"].values / cell).astype(int).astype(str)
            + "_" + np.floor(wells_meta["x"].values / cell).astype(int).astype(str))


def run_stratum(name, mat, meta, years, clustered=False, n_rep=REPLICATES):
    n_all = mat.shape[0]
    if n_all < MIN_STRAT_WELLS:
        return []
    if clustered:
        labels = block_labels(meta)
        uniq, counts = np.unique(labels, return_counts=True)
        keep_blocks = set(uniq[counts >= 10])
        mask = np.array([l in keep_blocks for l in labels])
        mat, meta, labels = mat[mask], meta[mask], labels[mask]
        n_all = mat.shape[0]
        if n_all < MIN_STRAT_WELLS:
            return []

    rows_idx_all = np.arange(n_all)
    ref_slope, ref_p = network_slope(mat, rows_idx_all)
    if np.isnan(ref_slope):
        return []
    ref_sign = 1 if ref_slope > 0 else -1
    ref_sig = ref_p < 0.05

    if clustered:
        # remove whole spatial blocks to reach each retention target
        block_groups = {}
        for b in np.unique(labels):
            block_groups[b] = np.where(labels == b)[0]

    out = []
    for f in RETENTION:
        n_keep = max(10, int(np.floor(f * n_all)))
        if n_keep >= n_all:
            continue
        flips = sig_losses = 0
        biases = []
        done = 0
        attempts = 0
        while done < n_rep and attempts < n_rep * 3:
            attempts += 1
            if clustered:
                keys = list(block_groups.keys())
                order = RNG.permutation(len(keys))
                sel, taken = [], 0
                for oi in order:
                    idx = block_groups[keys[oi]]
                    if taken + len(idx) <= max(n_keep * 1.15, n_keep + 10):
                        sel.append(idx)
                        taken += len(idx)
                    if taken >= n_keep * 0.9:
                        break
                keep = np.concatenate(sel) if len(sel) else rows_idx_all
            else:
                keep = RNG.choice(n_all, size=n_keep, replace=False)
            slope, p = network_slope(mat, keep)
            if np.isnan(slope):
                continue
            done += 1
            if (1 if slope > 0 else -1) != ref_sign:
                flips += 1
            if ref_sig and (p >= 0.05 or (1 if slope > 0 else -1) != ref_sign):
                sig_losses += 1
            biases.append(abs(slope - ref_slope))
        if done == 0:
            continue
        out.append(dict(
            stratum=name, n_all=n_all, retention=f, n_keep=n_keep,
            clustered=clustered, ref_slope=ref_slope, ref_p=ref_p,
            ref_sign=ref_sign, n_rep=done,
            p_sign_flip=flips / done,
            p_sig_loss=sig_losses / done,
            mean_abs_bias=float(np.mean(biases))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicates", type=int, default=REPLICATES)
    args = ap.parse_args()

    trends = pd.read_csv(os.path.join(DATA, "well_trends.csv"),
                         dtype={"monitoring_location_id": str})
    wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"),
                        dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    anoms = pd.read_csv(os.path.join(DATA, "annual_anoms.csv"),
                        dtype={"monitoring_location_id": str})

    wells = wells.join(trends.set_index("monitoring_location_id"),
                       how="inner").dropna(subset=["slope_ft_yr"])
    anoms = anoms[anoms.monitoring_location_id.isin(wells.index)]

    years = np.sort(anoms["year"].unique())
    years = years[(years >= 1975) & (years <= 2025)]

    aquifer = wells["national_aquifer_code"].fillna("UNKNOWN")
    counts = aquifer.value_counts()
    big = counts[counts >= MIN_STRAT_WELLS].index.tolist()
    print("strata with >=%d wells: %s" % (MIN_STRAT_WELLS, big), flush=True)

    results = []
    mat_all = build_matrix(anoms, wells.index, years)
    results += run_stratum("CONUS", mat_all, wells, years, clustered=False,
                           n_rep=args.replicates)
    results += run_stratum("CONUS", mat_all, wells, years, clustered=True,
                           n_rep=args.replicates)
    for aq in big:
        sel = wells.index[aquifer == aq]
        mat = build_matrix(anoms, sel, years)
        meta = wells.loc[sel]
        results += run_stratum(aq, mat, meta, years, clustered=False,
                               n_rep=args.replicates)

    res = pd.DataFrame(results)
    res.to_csv(os.path.join(DATA, "experiment_results.csv"), index=False)

    # N_min: smallest retained-network size with P(sign flip) <= 0.10
    th = []
    for (stratum, clus), grp in res.groupby(["stratum", "clustered"]):
        grp = grp.sort_values("n_keep")
        ok = grp[grp.p_sign_flip <= 0.10]
        th.append(dict(stratum=stratum, clustered=clus, n_all=grp.n_all.iloc[0],
                       n_min=int(ok.n_keep.min()) if len(ok) else None,
                       frac_min=(ok.n_keep.min() / grp.n_all.iloc[0]) if len(ok) else None,
                       ref_slope=grp.ref_slope.iloc[0],
                       ref_sig=bool(grp.ref_p.iloc[0] < 0.05)))
    pd.DataFrame(th).to_csv(os.path.join(DATA, "thresholds.csv"), index=False)
    print(res.to_string(max_rows=30))
    print("results ->", os.path.join(DATA, "experiment_results.csv"))


if __name__ == "__main__":
    main()
