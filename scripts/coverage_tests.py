#!/usr/bin/env python3
"""Coverage, composition, and spatial-autocorrelation tests of the law.

1  Coverage by flip regime: compute and archive the exact coverage numbers
        the abstract quotes (overall, low-flip, transition, high-flip) plus
        a t-CDF competitor fit (heavy tails predict a t link absorbs the
        probit residuals — test it; if it wins, report both).
2  Composition test: decompose the composite-hydrograph
        pipeline into (a) balanced-panel reference slopes (wells covering
        ≥80% of the window) vs (b) full composite; archive the sign
        agreement and the per-stratum β_balanced vs β_net ratio; test
        whether panel churn explains the c_i scatter (Spearman of
        c_i vs |log β ratio|).
3  Window sensitivity: regenerate with the law's own estimator,
        21 strata, archive correctly.
4  Spatial autocorrelation test: per-stratum Moran's I of
        well slopes at nearest-neighbor distance; test c_i vs I.
STRENGTH t-CDF link: fit scale on Φ and on t_ν (ν estimated); compare
        transition-zone coverage; if t wins materially, quote both.
CONVERGE n_min uncertainty: binomial CI on the crossing size for each
        stratum (not just the point) — archive nmin_lo / nmin_hi.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as st
from scipy.optimize import minimize_scalar
from scipy.stats import norm, t as tdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trends import yue_pilon_mk  # noqa: E402

DATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
Z10 = float(norm.ppf(0.9))


def fit_scale(x, flip, link="probit", nu=None):
    m = np.isfinite(x) & np.isfinite(flip) & (flip > 0) & (flip < 1)
    xx, yy = x[m], flip[m]
    if len(xx) < 5:
        return np.nan, np.nan

    def nll(s):
        if link == "probit":
            p = norm.cdf(-s * xx)
        else:
            p = tdist.cdf(-s * xx, nu)
        p = np.clip(p, 1e-9, 1 - 1e-9)
        return -np.sum(yy * np.log(p) + (1 - yy) * np.log(1 - p))
    r = minimize_scalar(nll, bounds=(0.1, 6), method="bounded")
    return r.x, nll(r.x)


def coverage(x, flip, scale, link="probit", nu=None):
    p = norm.cdf(-scale * x) if link == "probit" else tdist.cdf(-scale * x, nu)
    p = np.clip(p, 0.0, 1.0)
    n = 1000
    lo = st.binom.ppf(0.025, n, p) / n
    hi = st.binom.ppf(0.975, n, p) / n
    return float(((flip >= lo) & (flip <= hi)).mean())


def main():
    coll = pd.read_csv(os.path.join(DATA, "law_collapse.csv"))
    law = pd.read_csv(os.path.join(DATA, "law_nmin.csv"))
    out = {}

    # ---------- C1: coverage decomposition + t-CDF competitor ----------
    print("== C1: coverage by region + t-link competitor ==")
    C = 1.2720
    x = (C * coll.snr * coll.geom).values
    flip = coll.flip.values
    regions = dict(
        overall=(flip >= 0),
        low_flip=(flip < 0.05),
        transition=(flip >= 0.05) & (flip <= 0.5),
        high_flip=(flip > 0.5))
    for name, m in regions.items():
        out[f"coverage_{name}"] = coverage(x[m], flip[m], 1.0)
    print("  probit:", {k: round(v, 3) for k, v in out.items() if "coverage" in k})

    # t-CDF: fit nu on a grid, scale by MLE, compare transition coverage
    best = None
    for nu in [1, 2, 3, 5, 8, 15, 30]:
        s, nll = fit_scale(coll.snr.values * coll.geom.values, flip, "t", nu)
        xts = (s * coll.snr * coll.geom).values
        mt = (flip >= 0.05) & (flip <= 0.5)
        cov = coverage(xts[mt], flip[mt], 1.0, "t", nu)
        cov_all = coverage(xts, flip, 1.0, "t", nu)
        if best is None or nll < best[2]:
            best = (nu, s, nll, cov, cov_all)
    out["t_link"] = dict(nu=best[0], scale=round(best[1], 3),
                         transition_cov=round(best[3], 3), overall_cov=round(best[4], 3))
    print(f"  t-link best: nu={best[0]}, scale={best[1]:.3f}, "
          f"transition cov={best[3]:.3f} (probit: {out['coverage_transition']:.3f})")

    # ---------- C2: composition test ----------
    print("\n== C2: composition test (balanced panel vs composite) ==")
    trends = pd.read_csv(os.path.join(DATA, "well_trends.csv"),
                         dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"),
                        dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    w = wells.join(trends, how="inner").dropna(subset=["slope_ft_yr"])
    anoms = pd.read_csv(os.path.join(DATA, "annual_anoms.csv"),
                        dtype={"monitoring_location_id": str})
    years = np.arange(1975, 2026)
    rows = []
    for aq, g in w.groupby(w.national_aquifer_code.fillna("UNKNOWN")):
        if len(g) < 100:
            continue
        sub = anoms[anoms.monitoring_location_id.isin(g.index)]
        piv = sub.pivot_table(index="monitoring_location_id", columns="year",
                              values="resid", aggfunc="mean").reindex(columns=years)
        # composite (current pipeline)
        comp = piv.mean(axis=0, skipna=True)
        # balanced: wells covering >=80% of the window
        cov_frac = piv.notna().mean(axis=1)
        bal_ids = cov_frac[cov_frac >= 0.8].index
        bal = piv.loc[bal_ids].mean(axis=0, skipna=True) if len(bal_ids) >= 20 else comp
        # panel size churn
        panel = piv.notna().sum(axis=0)
        b_comp = yue_pilon_mk(comp.index.values, comp.values)[0] if comp.notna().sum() >= 10 else np.nan
        b_bal = yue_pilon_mk(bal.index.values, bal.values)[0] if bal.notna().sum() >= 10 else np.nan
        rows.append(dict(stratum=aq, n=len(g), n_balanced=len(bal_ids),
                         beta_net=b_comp, beta_balanced=b_bal,
                         sign_agree=bool(np.sign(b_comp) == np.sign(b_bal)),
                         panel_y0=int(panel.iloc[:5].mean()), panel_y1=int(panel.iloc[-5:].mean())))
    comp_test = pd.DataFrame(rows)
    comp_test["log_beta_ratio"] = np.log10(np.abs(comp_test.beta_net) /
                                            np.abs(comp_test.beta_balanced).replace(0, np.nan))
    comp_test.to_csv(os.path.join(DATA, "composition_test.csv"), index=False)
    ns = comp_test.sign_agree.sum()
    print(f"  sign agreement (balanced vs composite): {ns}/{len(comp_test)}")
    print(f"  median |log10 beta ratio|: {comp_test.log_beta_ratio.abs().median():.2f}")

    # c_i vs composition ratio (does churn explain the scatter?)
    dec = pd.read_csv(os.path.join(DATA, "decomposition.csv"))
    m2 = comp_test.merge(dec[["stratum", "c_i"]], on="stratum", how="inner").dropna(subset=["c_i"])
    if len(m2) >= 8:
        rho, p = st.spearmanr(m2.log_beta_ratio.abs(), m2.c_i)
        out["c_i_vs_composition"] = dict(spearman=round(float(rho), 3), p=round(float(p), 3), n=len(m2))
        print(f"  c_i vs |log beta ratio|: Spearman {rho:+.3f} (p={p:.3f}, n={len(m2)})")

    # ---------- C3: window sensitivity with law's estimator, 21 strata ----------
    print("\n== C3: window sensitivity (regenerated) ==")
    wrows = []
    for aq, g in w.groupby(w.national_aquifer_code.fillna("UNKNOWN")):
        if len(g) < 100:
            continue
        sub = anoms[anoms.monitoring_location_id.isin(g.index)]
        signs = []
        for y0 in (1970, 1975, 1980):
            hg = sub[sub.year >= y0].groupby("year").resid.mean()
            hg = hg[hg.index <= 2025]
            s0 = yue_pilon_mk(hg.index.values, hg.values)[0] if len(hg) >= 10 else np.nan
            signs.append(float(np.sign(s0)) if np.isfinite(s0) else np.nan)
        wrows.append(dict(stratum=aq, sign_1970=signs[0], sign_1975=signs[1], sign_1980=signs[2],
                          changes=bool(len(set(signs)) > 1)))
    ws = pd.DataFrame(wrows)
    ws.to_csv(os.path.join(DATA, "window_sensitivity.csv"), index=False)
    nch = ws.changes.sum()
    out["window_changes"] = int(nch)
    print(f"  strata changing sign: {nch}/{len(ws)} -> {', '.join(ws[ws.changes].stratum)}")

    # ---------- STRENGTH: Moran's I of well slopes ----------
    print("\n== Spatial autocorrelation (Moran's I, nearest 8 neighbors) ==")
    mrows = []
    for aq, g in w.groupby(w.national_aquifer_code.fillna("UNKNOWN")):
        if len(g) < 100:
            continue
        xy = g[["x", "y"]].values
        v = g.slope_ft_yr.values
        # nearest-neighbor mean via KD-tree
        from scipy.spatial import cKDTree
        tree = cKDTree(xy)
        _, idx = tree.query(xy, k=9)  # self + 8
        nb_mean = v[idx[:, 1:]].mean(axis=1)
        # Moran's I (permutation p on 200 shuffles for speed)
        zv = v - v.mean()
        I_obs = (zv * (nb_mean - v.mean())).sum() / (zv ** 2).sum() / (8 / 8)
        rng = np.random.default_rng(3)
        null = []
        for _ in range(200):
            vp = rng.permutation(v)
            zp = vp - vp.mean()
            nbp = vp[idx[:, 1:]].mean(axis=1)
            null.append((zp * (nbp - vp.mean())).sum() / (zp ** 2).sum())
        pval = float((np.array(null) >= I_obs).mean())
        mrows.append(dict(stratum=aq, morans_I=float(I_obs), p=pval))
    mo = pd.DataFrame(mrows)
    mo.to_csv(os.path.join(DATA, "spatial_autocorr.csv"), index=False)
    sig = mo[mo.p < 0.05]
    print(f"  strata with significant Moran's I: {len(sig)}/{len(mo)} "
          f"(median I = {mo.morans_I.median():.3f})")
    # c_i vs Moran's I
    m3 = mo.merge(dec[["stratum", "c_i"]], on="stratum", how="inner").dropna(subset=["c_i"])
    if len(m3) >= 8:
        rho, p = st.spearmanr(m3.morans_I, m3.c_i)
        out["c_i_vs_moransI"] = dict(spearman=round(float(rho), 3), p=round(float(p), 3), n=len(m3))
        print(f"  c_i vs Moran's I: Spearman {rho:+.3f} (p={p:.3f}, n={len(m3)})")

    # ---------- CONVERGE: n_min binomial CI ----------
    print("\n== n_min with binomial CI ==")
    nrows = []
    for _, r in law.iterrows():
        sub = coll[coll.stratum == r.stratum].sort_values("n")
        # find crossing: smallest n where flip <= 0.10 (point) and CI bounds
        def crossing(fn):
            ok = sub[sub.flip <= 0.10] if fn == "point" else None
            return int(ok.n.min()) if ok is not None and len(ok) else None
        n_pt = crossing("point")
        # CI on the crossing via bootstrap of flip values (binomial n=1000)
        rng = np.random.default_rng(7)
        cross_boot = []
        for _ in range(200):
            fl = rng.binomial(1000, np.clip(sub.flip.values, 0, 1)) / 1000
            okm = fl <= 0.10
            idx = np.where(okm)[0]
            cross_boot.append(int(sub.n.values[idx[0]]) if len(idx) else np.nan)
        cb = np.array([c for c in cross_boot if np.isfinite(c)])
        nrows.append(dict(stratum=r.stratum, N=r.N, nmin_point=n_pt,
                          nmin_lo=int(np.percentile(cb, 2.5)) if len(cb) else None,
                          nmin_hi=int(np.percentile(cb, 97.5)) if len(cb) else None))
    nci = pd.DataFrame(nrows)
    nci.to_csv(os.path.join(DATA, "nmin_confidence.csv"), index=False)
    have = nci.dropna(subset=["nmin_lo"])
    ratio = (have.nmin_hi / have.nmin_lo.clip(lower=1))
    print(f"  median n_min CI ratio (hi/lo): {ratio.median():.2f} across {len(have)} strata")

    json.dump(out, open(os.path.join(DATA, "coverage_tests.json"), "w"), indent=1)
    print("\nsaved -> data/coverage_tests.json")


if __name__ == "__main__":
    main()
