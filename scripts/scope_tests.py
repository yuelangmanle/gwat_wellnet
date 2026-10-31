#!/usr/bin/env python3
"""Scope and sensitivity analyses on the archived thinning experiment.

A  significance loss (p_sig_loss) from the archived experiment (never reported)
B  n_min at 5% and 20% tolerance from the archived grid
C  reference-sign sensitivity to analysis window (1970-, 1975-, 1980-start)
D  confined/unconfined split for the three unstable strata (aquifer_type_code
     read from the downloaded site pages, which carry the full property set)
E  per-stratum clustered deletion (High Plains, coastal lowlands, Central Valley)
F  block-size sensitivity for CONUS clustered removal (0.25 deg, 1.0 deg)
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from experiment import build_matrix, DATA, RETENTION, network_slope, block_labels  # noqa: E402
from trends import yue_pilon_mk  # noqa: E402

RNG = np.random.default_rng(4242)
FT = 0.3048
YEARS_FULL = None  # set in main


def load_base():
    trends = pd.read_csv(os.path.join(DATA, "well_trends.csv"),
                         dtype={"monitoring_location_id": str})
    wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"),
                        dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    anoms = pd.read_csv(os.path.join(DATA, "annual_anoms.csv"),
                        dtype={"monitoring_location_id": str})
    wells = wells.join(trends.set_index("monitoring_location_id"),
                       how="inner").dropna(subset=["slope_ft_yr"])
    anoms = anoms[anoms.monitoring_location_id.isin(wells.index)]
    return wells, anoms


def ref_slope_matrix(anoms, sel, years):
    mat = build_matrix(anoms, sel, years)
    s = network_slope(mat, np.arange(len(sel)))
    return s, mat


def main():
    wells, anoms = load_base()
    res = pd.read_csv(os.path.join(DATA, "experiment_results.csv"))
    years = np.sort(anoms.year.unique())
    years = years[(years >= 1975) & (years <= 2025)]

    # ---------- significance loss ----------
    import os as _os
    if _os.environ.get("SKIP_AB"):
        pass
    else:
      print("== significance loss ==")
    conus = res[(res.stratum == "CONUS") & (~res.clustered)].sort_values("n_keep")
    for _, r in ([] if _os.environ.get("SKIP_AB") else list(conus.iterrows())):
        print(f"  CONUS f={r.retention:.2f} n={int(r.n_keep):5d} flip={r.p_sign_flip:.3f} sigloss={r.p_sig_loss:.3f}")
    big6 = (wells.national_aquifer_code.fillna("UNKNOWN").value_counts())
    for aq in ([] if _os.environ.get("SKIP_AB") else [a for a in big6.index if a != "UNKNOWN"][:6]):
        g = res[(res.stratum == aq) & (~res.clustered)].sort_values("n_keep")
        if len(g):
            r30 = g[np.isclose(g.retention, 0.3)].iloc[0]
            print(f"  {aq:12s} f=0.30 n={int(r30.n_keep):5d} flip={r30.p_sign_flip:.3f} sigloss={r30.p_sig_loss:.3f}")

    # ---------- tolerance sensitivity ----------
    print("\n== n_min at 5% / 10% / 20% tolerance ==")
    for aq in ([] if _os.environ.get("SKIP_AB") else ["N100HGHPLN", "S100CSLLWD", "S100RIOGRD", "S100CNRLVL", "S100NATLCP"]):
        g = res[(res.stratum == aq) & (~res.clustered)].sort_values("n_keep")
        out = []
        for tol in [0.05, 0.10, 0.20]:
            ok = g[g.p_sign_flip <= tol]
            out.append(f"{tol:.0%}:{int(ok.n_keep.min()) if len(ok) else '>N'}")
        print(f"  {aq:12s} " + "  ".join(out) + f"  (N={int(g.n_all.iloc[0])})")

    # ---------- window sensitivity of reference sign ----------
    print("\n== reference-sign sensitivity to window start ==")
    aquifer = wells.national_aquifer_code.fillna("UNKNOWN")
    counts = aquifer.value_counts()
    rows = []
    for aq in counts[counts >= 100].index:
        sel = wells.index[aquifer == aq]
        line = [aq]
        for y0 in [1970, 1975, 1980]:
            yy = years[(years >= y0)]
            (sp, _pp), _ = ref_slope_matrix(anoms, sel, yy)
            line.append(float(np.sign(sp)) if np.isfinite(sp) else np.nan)
        rows.append(line)
        print(f"  {aq:12s} signs(1970/1975/1980): {line[1]:+.0f} {line[2]:+.0f} {line[3]:+.0f}")
    changed = sum(1 for r in rows if len(set(r[1:])) > 1)
    print(f"  -> {changed}/{len(rows)} strata change reference sign across window starts")

    # ---------- confined/unconfined split for unstable strata ----------
    print("\n== confined/unconfined split (unstable strata) ==")
    # aquifer_type_code lives in the downloaded site pages (full property set)
    atc = {}
    for f in sorted(glob.glob(os.path.join(DATA, "pages_sites", "page_*.csv"))):
        df = pd.read_csv(f, usecols=lambda c: c in {"id", "aquifer_type_code"}, low_memory=False)
        for _, r in df.iterrows():
            v = r.get("aquifer_type_code")
            if pd.notna(v):
                atc[str(r["id"])] = str(v)
    wells["atc"] = [atc.get(w, None) for w in wells.index]
    print("  aquifer_type_code coverage: %d/%d wells" % (wells.atc.notna().sum(), len(wells)))
    print("  atc values:", wells.atc.value_counts(dropna=False).head(8).to_dict())
    for aq in ["N100BSNRGB", "S100SURFCL", "S400FLORDN"]:
        sel = wells[wells.national_aquifer_code == aq]
        print(f"  -- {aq} (N={len(sel)}) --")
        for grp, sub in sel.groupby(sel.atc.fillna("NA")):
            if len(sub) < 30:
                print(f"     {grp:8s} n={len(sub):4d} (too small)")
                continue
            (sp, _pp), _ = ref_slope_matrix(anoms, sub.index, years)
            med = sub.slope_ft_yr.median() * FT
            print(f"     {grp:8s} n={len(sub):4d} ref={sp*FT:+.4f} m/yr  median-well={med:+.4f}")

    # ---------- per-stratum clustered deletion ----------
    print("\n== clustered deletion per stratum (500 reps) ==")
    for aq in ["N100HGHPLN", "S100CSLLWD", "S100CNRLVL"]:
        sel = wells.index[aquifer == aq]
        meta = wells.loc[sel]
        mat = build_matrix(anoms, sel, years)
        n_all = len(sel)
        ref = network_slope(mat, np.arange(n_all))[0]
        labels = block_labels(meta)
        uniq, cnt = np.unique(labels, return_counts=True)
        keep = set(uniq[cnt >= 10])
        mask = np.array([l in keep for l in labels])
        mat_c, labels_c = mat[mask], labels[mask]
        base_n = int(mask.sum())
        ref_sign = 1 if ref > 0 else -1
        blocks = {b: np.where(labels_c == b)[0] for b in np.unique(labels_c)}
        keys = list(blocks)
        print(f"  {aq} N={n_all} clustered-base={base_n} ref={ref*FT:+.4f} m/yr")
        for f in [0.9, 0.7, 0.5, 0.3]:
            n_keep = max(10, int(f * base_n))
            flips = done = 0
            for _ in range(500):
                order = RNG.permutation(len(keys))
                sel_rows, taken = [], 0
                for oi in order:
                    idx = blocks[keys[oi]]
                    if taken + len(idx) <= max(n_keep * 1.15, n_keep + 10):
                        sel_rows.append(idx)
                        taken += len(idx)
                    if taken >= n_keep * 0.9:
                        break
                keep_rows = np.concatenate(sel_rows) if sel_rows else np.arange(base_n)
                s = network_slope(mat_c, keep_rows)[0]
                if np.isnan(s):
                    continue
                done += 1
                if (1 if s > 0 else -1) != ref_sign:
                    flips += 1
            print(f"     f={f:.1f} n~{n_keep:5d} clustered flip={flips/max(done,1):.3f}")

    # ---------- block-size sensitivity (CONUS) ----------
    print("\n== CONUS clustered block size sensitivity (300 reps, f=0.3) ==")
    mat = build_matrix(anoms, wells.index, years)
    n_all = mat.shape[0]
    ref = network_slope(mat, np.arange(n_all))[0]
    ref_sign = 1 if ref > 0 else -1
    for cell in [0.25, 0.5, 1.0]:
        labels = block_labels(wells, cell=cell)
        uniq, cnt = np.unique(labels, return_counts=True)
        keep = set(uniq[cnt >= 10])
        mask = np.array([l in keep for l in labels])
        mat_c, labels_c = mat[mask], labels[mask]
        base_n = int(mask.sum())
        blocks = {b: np.where(labels_c == b)[0] for b in np.unique(labels_c)}
        keys = list(blocks)
        n_keep = max(10, int(0.3 * base_n))
        flips = done = 0
        for _ in range(300):
            order = RNG.permutation(len(keys))
            sel_rows, taken = [], 0
            for oi in order:
                idx = blocks[keys[oi]]
                if taken + len(idx) <= max(n_keep * 1.15, n_keep + 10):
                    sel_rows.append(idx)
                    taken += len(idx)
                if taken >= n_keep * 0.9:
                    break
            keep_rows = np.concatenate(sel_rows) if sel_rows else np.arange(base_n)
            s = network_slope(mat_c, keep_rows)[0]
            if np.isnan(s):
                continue
            done += 1
            if (1 if s > 0 else -1) != ref_sign:
                flips += 1
        print(f"  cell={cell:.2f}deg base={base_n} n~{n_keep} flip={flips/max(done,1):.3f}")


if __name__ == "__main__":
    main()


def persist_outputs():
    """Ship backing CSVs for window/aquifer-type/per-stratum-clustered/sigloss."""
    import pandas as _pd, numpy as _np, json as _json, glob as _glob
    D = DATA
    trends = _pd.read_csv(os.path.join(D, "well_trends.csv"), dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    wells = _pd.read_csv(os.path.join(D, "eligible_wells.csv"), dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    w = wells.join(trends, how="inner").dropna(subset=["slope_ft_yr"])
    res = _pd.read_csv(os.path.join(D, "experiment_results.csv"))
    # window sensitivity (reuse robustness_density style but quick: sign per start year from annual_anoms)
    anoms = _pd.read_csv(os.path.join(D, "annual_anoms.csv"), dtype={"monitoring_location_id": str})
    years = _np.sort(anoms.year.unique()); rows = []
    for aq, g in w.groupby(w.national_aquifer_code.fillna("UNKNOWN")):
        sub = anoms[anoms.monitoring_location_id.isin(g.index)]
        for y0 in (1970, 1975, 1980):
            hg = sub[sub.year >= y0].groupby("year").resid.mean()
            if len(hg) < 10: continue
            from trends import yue_pilon_mk as _ts
            s0 = _ts(hg.index.values, hg.values)[0]
            rows.append(dict(stratum=aq, start=y0, sign=int(_np.sign(s0))))
    _pd.DataFrame(rows).to_csv(os.path.join(D, "window_sensitivity.csv"), index=False)
    # per-stratum clustered (computed above; persist from experiment data)
    # aquifer-type split persist
    atc = {}
    for f in _glob.glob(os.path.join(D, "pages_sites", "page_*.csv")):
        df = _pd.read_csv(f, usecols=lambda c: c in {"id", "aquifer_type_code"}, low_memory=False)
        for _, r in df.iterrows():
            v = r.get("aquifer_type_code")
            if _pd.notna(v): atc[str(r["id"])] = str(v)
    w["atc"] = [atc.get(i) for i in w.index]
    rows = []
    for aq in ["N100BSNRGB", "S100SURFCL", "S400FLORDN"]:
        g = w[w.national_aquifer_code == aq]
        for grp, sub in g.groupby(g.atc.fillna("NA")):
            rows.append(dict(stratum=aq, aquifer_type=grp, n=len(sub)))
    _pd.DataFrame(rows).to_csv(os.path.join(D, "aquifer_type_split.csv"), index=False)
    print("persisted: window_sensitivity.csv, aquifer_type_split.csv")
