#!/usr/bin/env python3
"""Synthetic-population validation of the sampling law (out-of-domain
testbed) + forward Kish check of the clustered arm.

Part 1 — synthetic testbed: populations of wells with KNOWN slope structure
(normal / skewed / offsetting mixture), passed through the paper's own pipeline
(annual series -> average across wells -> Theil-Sen -> sign flip vs reference).
If the law is sound, noiseless populations give apparent c ~ 1; the question is
what skewness, offsetting composition, and temporal noise do to c.

Part 2 — forward two-stage check: predict the national clustered flip curve
from n_eff = n / Kish-DEFF (ICC 0.32, m_bar 41.5, cv 1.04 -> DEFF 29.06)
plugged into the calibrated law, and compare against the archived clustered
observations (no free parameters).
"""
import json
import os

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

rng = np.random.default_rng(42)
T = 51                     # 1975-2025 annual series
YEARS = np.arange(T, dtype=float)
N = 400                    # wells per population
SNR_POP = 0.15             # |mean slope| / SD of slopes (weak-signal stratum)
FRACT = [0.05, 0.08, 0.12, 0.18, 0.28, 0.45]
REPS = 500
NOISE_F = 0.3              # per-well slope SE = 0.3 * S (moderate temporal noise)


def theil_sen(y, x=YEARS):
    """O(T^2) pairwise slope median."""
    n = len(x)
    d = (y[:, None] - y[None, :]) / (x[:, None] - x[None, :])
    iu = np.triu_indices(n, k=1)
    return float(np.median(d[iu]))


def draw_slopes(shape, n, mean, sd):
    if shape == "normal":
        b = rng.normal(0.0, 1.0, n)
    elif shape == "skewed":
        b = rng.lognormal(0.0, 1.0, n)
    elif shape == "offsetting":
        comp = rng.random(n) < 0.6
        b = np.where(comp, rng.normal(1.0, 0.365, n), rng.normal(-0.9, 0.365, n))
    return (b - b.mean()) / b.std() * sd + mean


def run_population(shape, noise):
    slopes = draw_slopes(shape, N, SNR_POP, 1.0)          # S = 1 by construction
    intercepts = rng.normal(0.0, 5.0, N)
    sigma = (NOISE_F / 0.00951) if noise else 0.0          # annual-noise SD
    # annual series matrix (N x T)
    eps = rng.normal(0.0, sigma, (N, T)) if sigma else 0.0
    Y = intercepts[:, None] + slopes[:, None] * YEARS[None, :] + eps

    beta_net = theil_sen(Y.mean(axis=0))                   # paper's reference
    S_obs = float(np.std(slopes, ddof=0))                  # observed SD of well slopes
    snr = abs(beta_net) / S_obs

    rows = []
    for f in FRACT:
        n = max(10, int(round(f * N)))
        flips = 0
        for _ in range(REPS):
            idx = rng.choice(N, n, replace=False)
            sub = theil_sen(Y[idx].mean(axis=0))
            flips += (np.sign(sub) != np.sign(beta_net))
        rows.append({"n": n, "ret": n / N, "flip": flips / REPS,
                     "k": int(flips), "R": REPS})
    return beta_net, snr, rows


def fit_c(snr, rows):
    """Binomial MLE over informative trials (0 < k < R), as in law.py."""
    info = [r for r in rows if 0 < r["k"] < r["R"]]

    def nll(c):
        s = 0.0
        for r in info:
            x = c * snr * np.sqrt(r["n"] * N / (N - r["n"]))
            p = min(max(norm.cdf(-x), 1e-9), 1 - 1e-9)
            s += -(r["k"] * np.log(p) + (r["R"] - r["k"]) * np.log(1 - p))
        return s

    res = minimize_scalar(nll, bounds=(0.2, 6.0), method="bounded")
    # crude parametric CI from the informative-trial Fisher info
    return float(res.x), len(info)


def main():
    out = {"N": N, "T": T, "snr_pop": SNR_POP, "reps": REPS, "cells": {}}
    print(f"{'cell':28s} {'c_app':>6s} {'inf':>4s}  curve")
    for shape in ["normal", "skewed", "offsetting"]:
        for noise in [0.0, NOISE_F]:
            beta, snr, rows = run_population(shape, noise)
            c_hat, n_inf = fit_c(snr, rows)
            key = f"{shape}{'_noise' if noise else ''}"
            out["cells"][key] = {"beta_net": beta, "snr": snr, "c_app": c_hat,
                                 "informative": n_inf, "trials": rows}
            curve = " ".join(f"{r['ret']:.2f}:{r['flip']:.2f}" for r in rows)
            print(f"{key:28s} {c_hat:6.2f} {n_inf:4d}  {curve}")

    # ---- Part 2: forward Kish prediction of the national clustered curve ----
    import pandas as pd
    res = pd.read_csv(os.path.join(DATA, "experiment_results.csv"))
    cl = res[(res.stratum == "CONUS") & (res.clustered)]
    lc = pd.read_csv(os.path.join(DATA, "law_collapse.csv"))
    snr_conus = float(lc[lc.stratum == "CONUS"].snr.iloc[0])
    Nnat = 25754
    DEFF = (1 + 40.5 * 0.32) * (1 + 1.04 ** 2)     # Kish: 29.06, no free params
    c_conus = 1.27
    fwd = []
    for _, r in cl.iterrows():
        n = int(r.n_keep)
        if n < 20:
            continue
        n_eff = n / DEFF
        x = c_conus * snr_conus * np.sqrt(n_eff * Nnat / (Nnat - n_eff))
        p_pred = float(norm.cdf(-x))
        k, R = r.p_sign_flip, 1.0                   # p is a mean over replicates
        fwd.append({"ret": float(r.retention), "n": n, "obs": float(r.p_sign_flip),
                    "pred": p_pred,
                    "lo": float(max(0.0, r.p_sign_flip - 1.96 * np.sqrt(
                        max(r.p_sign_flip, 1e-4) * (1 - r.p_sign_flip) / 1000))),
                    "hi": float(min(1.0, r.p_sign_flip + 1.96 * np.sqrt(
                        max(r.p_sign_flip, 1e-4) * (1 - r.p_sign_flip) / 1000)))})
    inb = sum(1 for f in fwd if f["lo"] <= f["pred"] <= f["hi"])
    out["kish_forward"] = {"deff": DEFF, "snr": snr_conus, "c": c_conus,
                           "points": fwd, "within_bands": f"{inb}/{len(fwd)}"}
    print(f"\nKish forward (no free params): within binomial bands {inb}/{len(fwd)}")
    for f in fwd:
        print(f"  ret={f['ret']:.2f} n={f['n']:6d} obs={f['obs']:.3f} pred={f['pred']:.3f}")

    with open(os.path.join(DATA, "sim_validation.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote data/sim_validation.json")


if __name__ == "__main__":
    main()
