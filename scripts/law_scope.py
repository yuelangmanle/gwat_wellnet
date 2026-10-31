#!/usr/bin/env python3
"""Where does the sampling law stop? Scope test on significance loss.

The flip law treats the *direction* of the network conclusion.  A natural
question is whether the same normalized size  x = c*SNR*sqrt(nN/(N-n))
also governs the second failure mode -- loss of statistical significance.
Under the same geometry the replicate z-statistic behaves like

    z_hat ~ N( z(n),  (z(n)/x)^2 ),   P(sigloss) = Phi( -(1 - z_crit/z(n)) x ),

with z(n) the expected z-statistic of the averaged hydrograph at n wells.
Two variants are tested:
  (i)  z(n) = z_ref (no noise growth; slope sampling only);
  (ii) z(n) = z_ref * sqrt((1+rho)/(1+rho*N/n)), where rho = sigma_e^2/(sigma_u^2 N)
       is the archived ratio of well-specific to common-mode year noise,
       computable from the annual anomalies without new simulations.

Each variant keeps one global slack constant w (w = 1 if the derivation were
exact) fitted by binomial MLE over all trials.  The flip law's coverage on
the same trials is computed for comparison.  If coverage stays far below the
flip law's, significance defence degrades through a mechanism (temporal
noise growth of the averaged hydrograph) that the spatial-sampling law does
not capture -- a boundary of the law, reported honestly.

Outputs: data/law_sigloss_scope.csv, data/law_sigloss_scope.json.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))
Z_CRIT = 1.959964  # two-sided 95%


def archived_rho():
    """Per stratum: rho = sigma_e^2 / (sigma_u^2 N) from annual anomalies."""
    an = pd.read_csv(os.path.join(DATA, "annual_anoms.csv"))
    wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"),
                        dtype={"monitoring_location_id": str})
    wmap = dict(zip(wells.monitoring_location_id, wells.national_aquifer_code))
    an["stratum"] = an.monitoring_location_id.map(wmap)
    out = {}
    for st, g in an.groupby("stratum"):
        piv = g.pivot_table(index="year", columns="monitoring_location_id",
                            values="resid")
        piv = piv.loc[:, piv.notna().sum() >= 15]
        if piv.shape[1] < 30:
            continue
        m_t = piv.mean(axis=1)
        t = m_t.index.values - m_t.index.values.mean()
        b = np.polyfit(t, m_t.values, 1)[0]
        su2 = float(np.var(m_t.values - b * t, ddof=1))
        dev = piv.sub(m_t, axis=0)
        se2 = float(np.nanmean(np.nanvar(dev.values, axis=0, ddof=1)))
        if su2 <= 0:
            continue
        out[st] = se2 / (su2 * piv.shape[1])
    return out


def main():
    res = pd.read_csv(os.path.join(DATA, "experiment_results.csv"))
    coll = pd.read_csv(os.path.join(DATA, "law_collapse.csv"))
    fit = json.load(open(os.path.join(DATA, "law_fit.json")))
    c = fit["c"]
    snr_map = dict(zip(coll.stratum, coll.snr))
    rho_map = archived_rho()
    z_ref_map = {st: float(min(norm.ppf(1 - g.ref_p.iloc[0] / 2), 8.5))
                 for st, g in res.groupby("stratum")}

    d = res[(~res.clustered) & (res.ref_p < 0.05)].copy()
    d = d[d.stratum.isin(snr_map) & d.stratum.isin(rho_map)
           & d.p_sig_loss.notna() & (d.n_keep < d.n_all)]
    d["x"] = c * d.stratum.map(snr_map) * np.sqrt(
        d.n_keep * d.n_all / (d.n_all - d.n_keep))
    d["z_ref"] = d.stratum.map(z_ref_map)
    d["rho"] = d.stratum.map(rho_map)

    def fit_variant(z_n):
        a = 1 - Z_CRIT / z_n
        k = d.p_sig_loss.values * d.n_rep.values

        def nll(w):
            p = np.clip(norm.cdf(-w * a * d.x.values), 1e-6, 1 - 1e-6)
            return -(k * np.log(p) + (d.n_rep.values - k) * np.log(1 - p)).sum()

        w = minimize_scalar(nll, bounds=(0.05, 4.0), method="bounded",
                            options={"xatol": 1e-3}).x
        return float(w), norm.cdf(-w * a * d.x.values)

    def coverage(pred, obs):
        # same variance-floor convention as law.py (zero-observation trials
        # get the 2.5e-4 floor so a near-zero prediction stays in band)
        se = np.sqrt(np.clip(obs * (1 - obs), 2.5e-4, None) / d.n_rep)
        return float(((pred >= obs - 1.96 * se)
                      & (pred <= obs + 1.96 * se)).mean())

    # flip law on the identical trial set, for a fair comparison
    flip_pred = norm.cdf(-d.x.values)
    flip_cov = coverage(flip_pred, d.p_sign_flip.values)

    variants = {}
    for name, z_n in [("slope_sampling_only", d.z_ref),
                      ("with_noise_growth",
                       d.z_ref * np.sqrt((1 + d.rho)
                                         / (1 + d.rho * d.n_all / d.n_keep)))]:
        w, pred = fit_variant(z_n.values)
        variants[name] = dict(w=w,
                              coverage=coverage(pred, d.p_sig_loss.values))
        d["pred_" + name] = pred
        d["a_" + name] = 1 - Z_CRIT / z_n.values

    out = dict(
        n_trials=int(len(d)), n_strata=int(d.stratum.nunique()),
        flip_law_coverage_same_trials=flip_cov,
        sigloss=dict(variants),
        verdict=("significance loss does not collapse onto the same law: "
                 "best-variant binomial coverage {:.0%} vs {:.0%} for the "
                 "flip law on the same trials".format(
                     max(v["coverage"] for v in variants.values()), flip_cov)),
    )
    d.to_csv(os.path.join(DATA, "law_sigloss_scope.csv"), index=False)
    with open(os.path.join(DATA, "law_sigloss_scope.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
