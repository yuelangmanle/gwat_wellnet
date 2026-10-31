#!/usr/bin/env python3
"""A finite-population sampling law for trend-sign flips.

Theory.  The network-trend estimator behaves like the mean of the retained
wells' slopes.  Thinning samples n of N wells WITHOUT replacement, so the
mean of the retained slopes has variance  S^2 (N - n) / (n N)  with S the
finite-population spread of well slopes.  Hence

    P(flip) = Phi( -c * SNR * sqrt( n N / (N - n)) ),
    SNR     = |beta_net| / S,   c = calibration constant (~1).

Inversion at tolerance alpha (z = Phi^-1(1 - alpha)) gives the closed-form
minimum reliable size

    n_min = N z^2 / ( z^2 + (c SNR)^2 N ).

Outputs: data/law_collapse.csv, data/law_nmin.csv, data/law_fit.json.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))
ALPHA = 0.10
Z = float(norm.ppf(1 - ALPHA))


def load():
    res = pd.read_csv(os.path.join(DATA, "experiment_results.csv"))
    wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"),
                        dtype={"monitoring_location_id": str})
    trends = pd.read_csv(os.path.join(DATA, "well_trends.csv"),
                         dtype={"monitoring_location_id": str})
    return res, wells.merge(trends, on="monitoring_location_id", how="inner")


def well_slopes(w, aq):
    if aq == "CONUS":
        return w.slope_ft_yr.values
    if aq == "UNKNOWN":
        return w[w.national_aquifer_code.isna()].slope_ft_yr.values
    return w[w.national_aquifer_code == aq].slope_ft_yr.values


def build_frame(res, w):
    """One row per (stratum, retention, random mode) with SNR and geometry."""
    rows = []
    for aq, g in res[~res.clustered].groupby("stratum"):
        s = well_slopes(w, aq)
        if len(s) < 50:
            continue
        S = float(np.std(s, ddof=0))
        ref = float(g.ref_slope.iloc[0])
        snr = abs(ref) / S
        N = int(g.n_all.iloc[0])
        # retentions mapping to the same n_keep (small strata) are one size:
        # average the replicate flip estimates, sum the replicate counts
        gg = g.groupby("n_keep").agg(flip=("p_sign_flip", "mean"),
                                     n_rep=("n_rep", "sum")).reset_index()
        for _, r in gg.iterrows():
            n = int(r.n_keep)
            rows.append(dict(
                stratum=aq, N=N, n=n, retention=n / N,
                snr=snr, geom=np.sqrt(n * N / (N - n)),
                flip=r.flip, n_rep=int(r.n_rep),
                n_wells_with_slope=len(s)))
    return pd.DataFrame(rows)


def pred_flip(c, df):
    return norm.cdf(-c * df.snr * df.geom)


def fit_c(df, exclude=None):
    """Binomial MLE of the global calibration constant c."""
    d = df if exclude is None else df[df.stratum != exclude]
    ok = (d.flip > 0) & (d.flip < 1) & (d.n < d.N)  # n=N rows carry no signal
    d = d[ok]
    k = d.flip * d.n_rep

    def nll(c):
        p = np.clip(norm.cdf(-c * d.snr * d.geom), 1e-9, 1 - 1e-9)
        return -float(np.sum(k * np.log(p) + (d.n_rep - k) * np.log(1 - p)))

    fit = minimize_scalar(nll, bounds=(0.5, 3.0), method="bounded",
                          options={"xatol": 1e-4})
    return float(fit.x), float(fit.fun)


def predict_nmin(c, snr, N):
    lam = (c * snr) ** 2 * N
    return Z ** 2 * N / (Z ** 2 + lam)


def main():
    res, w = load()
    df = build_frame(res, w)
    c, nll = fit_c(df)

    # bootstrap CI over strata (cluster bootstrap)
    rng = np.random.default_rng(11)
    strata = df.stratum.unique()
    boots = []
    for _ in range(400):
        pick = rng.choice(strata, size=len(strata), replace=True)
        parts = [df[df.stratum == s] for s in pick]
        cb, _ = fit_c(pd.concat(parts, ignore_index=True))
        boots.append(cb)
    ci = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))

    # goodness of fit: observed vs predicted, binomial bands
    df["pred"] = pred_flip(c, df)
    se = np.sqrt(np.clip(df.flip * (1 - df.flip), 2.5e-4, None) / df.n_rep)
    df["in_ci"] = (df.flip - df.pred).abs() <= 1.96 * se
    in_ci = float(df.in_ci.mean())

    # leave-one-stratum-out: refit c, score held-out flip predictions
    loo_ll, loo_in = [], []
    for s in strata:
        c_loo, _ = fit_c(df, exclude=s)
        d = df[df.stratum == s]
        p = np.clip(pred_flip(c_loo, d), 1e-9, 1 - 1e-9)
        k = d.flip * d.n_rep
        loo_ll.append(-float(np.sum(k * np.log(p) + (d.n_rep - k) * np.log(1 - p))))
        loo_in.append(float(((d.flip - p).abs() <= 1.96 * np.sqrt(
            np.clip(d.flip * (1 - d.flip), 2.5e-4, None) / d.n_rep)).mean()))
    loo_mean_in = float(np.mean(loo_in))

    # simulated vs predicted n_min
    nmin_rows = []
    for aq, g in res[~res.clustered].groupby("stratum"):
        d = df[df.stratum == aq]
        if not len(d):
            continue
        ok = g[g.p_sign_flip <= ALPHA]
        N = int(g.n_all.iloc[0])
        nmin_sim = float(ok.n_keep.min()) if len(ok) else np.nan
        snr = float(d.snr.iloc[0])
        nmin_pred = predict_nmin(c, snr, N)
        nmin_rows.append(dict(
            stratum=aq, N=N, snr=snr, S_ft=float(np.std(well_slopes(w, aq), ddof=0)),
            ref_slope_ft=float(g.ref_slope.iloc[0]),
            nmin_sim=nmin_sim, censored=bool(np.isnan(nmin_sim)),
            nmin_pred=float(nmin_pred), ret_pred=float(nmin_pred / N),
            c_loo=predict_nmin(fit_c(df, exclude=aq)[0], snr, N)))
    nmin = pd.DataFrame(nmin_rows)
    uncens = nmin[~nmin.censored]
    # log retention error for uncensored, direction check for censored
    err = np.log10(uncens.nmin_sim / uncens.N) - np.log10(uncens.nmin_pred / uncens.N)
    med_err = float(np.median(np.abs(err)))
    cens_ok = int(((nmin[nmin.censored].ret_pred > 0.86)).sum())
    cens_n = int(nmin.censored.sum())

    # clustered inversion (CONUS): effective n under clustered deletion
    cl = res[(res.stratum == "CONUS") & res.clustered].sort_values("n_keep")
    snr_conus = float(nmin.loc[nmin.stratum == "CONUS", "snr"].iloc[0])
    Nc = int(cl.n_all.max())
    lam = (c * snr_conus) ** 2 * Nc
    inv = []
    for _, r in cl.iterrows():
        if not (0.001 < r.p_sign_flip < 0.999):
            continue
        z = float(norm.ppf(1 - r.p_sign_flip))
        neff = Nc * z ** 2 / (z ** 2 + lam)
        inv.append(dict(retention=r.retention, n_keep=int(r.n_keep),
                        flip=r.p_sign_flip, n_eff=neff,
                        deff=r.n_keep / max(neff, 1e-9)))
    inv = pd.DataFrame(inv)
    band = inv[(inv.retention >= 0.10) & (inv.retention <= 0.30)]
    neff_cap = (float(band.n_eff.min()), float(band.n_eff.max()))
    deff_30 = float(inv.loc[inv.retention == 0.30, "deff"].iloc[0]) \
        if (inv.retention == 0.30).any() else float("nan")

    out = dict(
        c=c, ci=ci, z=Z, alpha=ALPHA,
        in_ci_band=in_ci, loo_mean_in_ci=loo_mean_in,
        n_points=int(len(df)), n_strata=int(len(strata)),
        nmin_med_log10_err=med_err,
        nmin_within_factor2=float(((uncens.nmin_sim / uncens.nmin_pred).between(0.5, 2)).mean()),
        censored_pred_above_86pct=f"{cens_ok}/{cens_n}",
        conus=snr_conus,
        clustered=dict(N_clustered=Nc, neff_cap_10_30=neff_cap, deff_at_30pct=deff_30,
                       max_deff=float(inv.deff.max())),
    )
    df.to_csv(os.path.join(DATA, "law_collapse.csv"), index=False)
    nmin.to_csv(os.path.join(DATA, "law_nmin.csv"), index=False)
    inv.to_csv(os.path.join(DATA, "law_clustered.csv"), index=False)
    with open(os.path.join(DATA, "law_fit.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
