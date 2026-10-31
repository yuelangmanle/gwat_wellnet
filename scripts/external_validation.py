#!/usr/bin/env python3
"""Validation of the constant and the frequency axis:

A. External validation of c (two strict paths within the archive):
   A1 hold-out: fit c on 11 random strata, predict flip for the other 11,
      score binomial coverage on the held-out half (repeat 200 splits).
   A2 estimator robustness: recompute each stratum's SNR with OLS on annual
      means (instead of Theil-Sen), refit c on the SAME experiments; if
      c_OLS lands near 1.27 the law is pipeline-robust.

B. Frequency axis: restrict every well to one measurement per year (the
   most common month per well), recompute all well slopes, stratum SNR,
   and the law's n_min; compare with the full-measurement baseline.

C. Transition-zone coverage: regress per-stratum c_i on kurtosis /
   skewness / max leverage of the well-slope distribution; if significant,
   quote a heavy-tail correction rule.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as st
from scipy.optimize import minimize_scalar
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
Z = float(norm.ppf(0.9))
C_GLB = json.load(open(os.path.join(DATA, "law_fit.json")))["c"]


def fit_c_binom(x, flip):
    m = np.isfinite(flip) & (flip > 0) & (flip < 1) & np.isfinite(x)
    if m.sum() < 5:
        return np.nan
    xx, yy = x[m], flip[m]

    def nll(cc):
        p = np.clip(norm.cdf(-cc * xx), 1e-9, 1 - 1e-9)
        return -np.sum(yy * np.log(p) + (1 - yy) * np.log(1 - p))
    r = minimize_scalar(nll, bounds=(0.3, 6), method="bounded")
    return r.x


def in_ci(x, flip, c):
    p = norm.cdf(-c * x)
    n = 1000
    lo = st.binom.ppf(0.025, n, np.clip(p, 1e-9, 1)) / n
    hi = st.binom.ppf(0.975, n, np.clip(p, 1e-9, 1)) / n
    return (flip >= lo) & (flip <= hi)


def main():
    coll = pd.read_csv(os.path.join(DATA, "law_collapse.csv"))
    strata = sorted(coll.stratum.unique())

    # ---------- A1: hold-out validation ----------
    print("== A1 hold-out: fit on 11 strata, predict the other 11 (200 splits) ==")
    rng = np.random.default_rng(7)
    covs, cs = [], []
    for rep in range(200):
        perm = rng.permutation(len(strata))
        train = [strata[i] for i in perm[:11]]
        test = [strata[i] for i in perm[11:]]
        tr = coll[coll.stratum.isin(train)]
        c_hat = fit_c_binom((tr.snr * tr.geom).values, tr.flip.values)
        if not np.isfinite(c_hat):
            continue
        te = coll[coll.stratum.isin(test)]
        x = (te.snr * te.geom).values
        ok = in_ci(x, te.flip.values, c_hat)
        covs.append(ok.mean())
        cs.append(c_hat)
    print(f"  held-out binomial coverage: mean {np.mean(covs):.3f} "
          f"(in-sample reference 0.486; perfect calibration 0.95)")
    print(f"  split-level c: median {np.median(cs):.2f}, range "
          f"{np.percentile(cs,5):.2f}-{np.percentile(cs,95):.2f}")
    json.dump(dict(holdout_coverage_mean=float(np.mean(covs)),
                   holdout_c_median=float(np.median(cs))),
              open(os.path.join(DATA, "external_validation.json"), "w"), indent=1)

    # ---------- A2: estimator robustness (OLS on annual means) ----------
    print("\n== A2 estimator robustness: SNR from OLS on annual means ==")
    trends = pd.read_csv(os.path.join(DATA, "well_trends.csv"),
                         dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"),
                        dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    w = wells.join(trends, how="inner")
    anoms = pd.read_csv(os.path.join(DATA, "annual_anoms.csv"),
                        dtype={"monitoring_location_id": str})
    # OLS slope per well on annual anomaly medians (cheap vectorised)
    g = anoms.groupby("monitoring_location_id")
    ols = g.apply(lambda d: np.polyfit(d.year, d.resid, 1)[0]
                  if len(d) >= 5 else np.nan)
    w["slope_ols"] = ols
    rows = []
    for aq, gg in w.groupby(w.national_aquifer_code.fillna("UNKNOWN")):
        if len(gg) < 100:
            continue
        sd_ols = gg.slope_ols.std()
        mean_ts = gg.slope_ft_yr.mean()  # reference unchanged (TS-of-hydrograph)
        # SNR under OLS spread: use |mean of OLS slopes| / SD_OLS as the
        # operator-visible variant, and refit c on existing flips
        snr_ols = abs(gg.slope_ols.mean()) / sd_ols if sd_ols > 0 else np.nan
        rows.append(dict(stratum=aq, N=len(gg), snr_ts=abs(mean_ts) / gg.slope_ft_yr.std(),
                         snr_ols=snr_ols, kurt=float(st.kurtosis(gg.slope_ft_yr, fisher=True)),
                         skew=float(st.skew(gg.slope_ft_yr)),
                         lev=float(((gg.slope_ft_yr - gg.slope_ft_yr.mean()) ** 2).max()
                                   / ((gg.slope_ft_yr - gg.slope_ft_yr.mean()) ** 2).sum())))
    ext = pd.DataFrame(rows)
    # per-stratum c_i (from earlier decomposition-style fit on the archive flips)
    cis = []
    for _, r in ext.iterrows():
        sub = coll[coll.stratum == r.stratum]
        if len(sub) < 4:
            cis.append(np.nan)
            continue
        cis.append(fit_c_binom((sub.snr * sub.geom).values, sub.flip.values))
    ext["c_i"] = cis
    ext.to_csv(os.path.join(DATA, "external_validation.csv"), index=False)

    # c with OLS-based SNR on the same flip data
    coll2 = coll.merge(ext[["stratum", "snr_ols"]], on="stratum", how="inner")
    x = coll2.snr_ols.values * coll2.geom.values
    c_ols = fit_c_binom(x, coll2.flip.values)
    print(f"  c under OLS-spread SNR: {c_ols:.2f}  (Theil-Sen baseline {C_GLB:.2f})")

    # ---------- C: transition-zone/heavy-tail rule ----------
    print("\n== C: c_i vs distribution shape ==")
    ok = ext.dropna(subset=["c_i", "kurt"])
    for var in ["kurt", "skew", "lev"]:
        r, p = st.spearmanr(ok[var], ok.c_i)
        print(f"  Spearman(c_i, {var}) = {r:+.2f}  (p={p:.3f})")
    # heavy-tail rule: kurtosis > 30 -> median c_i of that group vs rest
    ht = ok[ok.kurt > 30]
    rest = ok[ok.kurt <= 30]
    if len(ht) >= 3:
        print(f"  kurt>30 group: n={len(ht)}, median c_i={ht.c_i.median():.2f} "
              f"vs rest {rest.c_i.median():.2f}")
        # coverage gain: global c for low-kurt strata, group median for high-kurt
        low = coll[coll.stratum.isin(rest.stratum)]
        c_low = fit_c_binom((low.snr * low.geom).values, low.flip.values)
        cov_global = in_ci((low.snr * low.geom).values, low.flip.values, C_GLB).mean()
        cov_fit = in_ci((low.snr * low.geom).values, low.flip.values, c_low).mean()
        print(f"  low-kurt strata: c={c_low:.2f}, coverage {cov_global:.2f} (global c) "
              f"-> {cov_fit:.2f} (own c)")
        json.dump(dict(c_ols=float(c_ols), c_low=float(c_low),
                       kurt_split=dict(n_ht=len(ht), c_ht=float(ht.c_i.median()),
                                       c_rest=float(rest.c_i.median()),
                                       cov_global=float(cov_global),
                                       cov_own=float(cov_fit))),
                  open(os.path.join(DATA, "transition_rule.json"), "w"), indent=1)

    # ---------- B: frequency axis (separate stage, needs trends recompute) ----------
    print("\n(B runs in stage 2: see frequency_axis.py)")


if __name__ == "__main__":
    main()
