#!/usr/bin/env python3
"""Methodological robustness of the reference trend and the law.

1  Density-weighted reference: recompute every stratum's reference slope
     with inverse-local-density weights (each well weighted by 1 / number of
     eligible wells in its 0.5-degree block) and report sign agreement with
     the unweighted reference.
2  Bootstrap CI of reference slopes: resample wells with replacement
     (1,000x), recompute the network slope, report the 2.5/97.5 percentiles
     and whether 0 lies inside — quantifies "the aquifer-scale trend is not
     yet estimable" for unstable strata.
5  Leave-one-out validation of the SNR logistic law: fit the logistic on
     20 strata, predict flip@30% for the held-out stratum, and count
     risk-class misclassifications (threshold 10%).
Outputs: data/robustness_density.csv, data/robustness_bootstrap.csv,
data/robustness_loo.csv (+ printed summary).
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from experiment import build_matrix, DATA  # noqa: E402
from trends import yue_pilon_mk  # noqa: E402

RNG = np.random.default_rng(20261023)


def net_slope(mat, rows, weights=None):
    sub = mat[rows]
    if weights is not None:
        w = weights[rows][:, None]
        series = np.nansum(sub * w, axis=0) / np.where(np.nansum((~np.isnan(sub)) * w, axis=0) == 0,
                                                       np.nan, np.nansum((~np.isnan(sub)) * w, axis=0))
    else:
        series = np.nanmean(sub, axis=0)
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

    # ---- density weighting ----
    dens_rows = []
    boot_rows = []
    for name, sel in strata:
        meta = wells.loc[sel]
        mat = build_matrix(anoms, sel, years)
        n = len(sel)
        s_unw = net_slope(mat, np.arange(n))
        # inverse-density weights: 1 / wells within the same 0.5-degree block
        bx = np.floor(meta.x.values / 0.5).astype(int)
        by = np.floor(meta.y.values / 0.5).astype(int)
        keys = list(zip(by, bx))
        cnt = pd.Series(keys).map(pd.Series(keys).value_counts()).values
        wgt = 1.0 / cnt
        wgt = wgt / wgt.sum() * n  # normalize to effective n
        s_den = net_slope(mat, np.arange(n), weights=wgt)
        dens_rows.append(dict(stratum=name, n=n, slope_unw=s_unw, slope_den=s_den,
                              sign_agree=bool(np.sign(s_unw) == np.sign(s_den))))
        # ---- bootstrap CI of reference ----
        slopes = []
        for _ in range(1000):
            rows = RNG.choice(n, size=n, replace=True)
            slopes.append(net_slope(mat, rows))
        sl = np.array([x for x in slopes if np.isfinite(x)])
        boot_rows.append(dict(stratum=name, n_boot=len(sl),
                              lo=np.percentile(sl, 2.5), hi=np.percentile(sl, 97.5),
                              ci_excludes_zero=not (np.percentile(sl, 2.5) <= 0 <= np.percentile(sl, 97.5))))
    dens = pd.DataFrame(dens_rows)
    dens.to_csv(os.path.join(DATA, "robustness_density.csv"), index=False)
    boot = pd.DataFrame(boot_rows)
    boot.to_csv(os.path.join(DATA, "robustness_bootstrap.csv"), index=False)

    # ---- LOO of the SNR law ----
    res = pd.read_csv(os.path.join(DATA, "experiment_results.csv"))
    spread = wells.groupby(wells.national_aquifer_code.fillna("UNKNOWN")).slope_ft_yr.apply(
        lambda s: float(np.median(np.abs(s - np.median(s)))))
    pts = []
    for aq, g in res[(~res.clustered) & (res.stratum != "CONUS")].groupby("stratum"):
        g30 = g[np.isclose(g.retention, 0.3)]
        if len(g30):
            sd = spread.get(aq, np.nan)
            if np.isfinite(sd) and sd > 0:
                pts.append((aq, abs(g30.ref_slope.iloc[0]) / sd, g30.p_sign_flip.iloc[0]))
    pts = pd.DataFrame(pts, columns=["aq", "snr", "flip"])

    def logistic(x, b, x0):
        return 1.0 / (1.0 + np.exp(b * (x - x0)))

    loo = []
    for i in range(len(pts)):
        tr = pts.drop(i)
        te = pts.iloc[i]
        try:
            popt, _ = curve_fit(logistic, tr.snr, tr.flip, p0=[4, 0.2], maxfev=5000)
            pred = float(logistic(te.snr, *popt))
        except Exception:
            pred = np.nan
        loo.append(dict(aq=te.aq, snr=te.snr, observed=te.flip, predicted=pred,
                        obs_class=int(te.flip > 0.10), pred_class=int(pred > 0.10)))
    loo = pd.DataFrame(loo)
    loo.to_csv(os.path.join(DATA, "robustness_loo.csv"), index=False)

    print("== density weighting ==")
    print("sign agreement: %d/%d strata" % (dens.sign_agree.sum(), len(dens)))
    print(dens[~dens.sign_agree].to_string(index=False))
    print("\n== bootstrap CI of reference ==")
    for _, r in boot.iterrows():
        flag = "" if r.ci_excludes_zero else "  <-- CI includes zero"
        print(f"{r.stratum:12s} n_boot={r.n_boot:4d} CI [{r.lo:+.4f}, {r.hi:+.4f}] ft/yr{flag}")
    print("\n== LOO of logistic ==")
    acc = (loo.obs_class == loo.pred_class).mean()
    print("LOO risk-class agreement: %d/%d (%.0f%%)" % (
        int((loo.obs_class == loo.pred_class).sum()), len(loo), 100 * acc))
    print("mean |predicted - observed| flip: %.3f" % loo.predicted.sub(loo.observed).abs().mean())
    # full-fit parameters for the caption
    popt, _ = curve_fit(logistic, pts.snr, pts.flip, p0=[4, 0.2], maxfev=5000)
    print("\nfull logistic: flip = 1/(1+exp(%.2f*(SNR-%.2f)))" % tuple(popt))


if __name__ == "__main__":
    main()
