#!/usr/bin/env python3
"""Decomposition experiment: benchmark versus pipeline estimator.

Benchmark: draw SRS subsets directly from each stratum's OBSERVED well-slope
population and track the sign of the subset MEAN against the full-population
mean. This is the exact design-based law with c = 1:
    P_mean(flip) = Phi(-SNR_mean * sqrt(nN/(N-n))),  SNR_mean = |mean(b)|/SD(b).
Comparing the actual experiment's curve (Theil-Sen of mean hydrograph) to this
benchmark turns the per-stratum calibration into a measured functional gap:
    c_i = x_actual / x_benchmark  at matched flip probability,
where x = SNR*sqrt(nN/(N-n)) evaluated at each trial.

Outputs: data/decomposition.csv (per stratum: c_i measured, mean-benchmark
coverage, functional gap), data/decomposition_summary.json.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

DATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
RNG = np.random.default_rng(31337)
RET = [0.9, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01]


def main():
    trends = pd.read_csv(os.path.join(DATA, "well_trends.csv"),
                         dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"),
                        dtype={"monitoring_location_id": str}).set_index("monitoring_location_id")
    w = wells.join(trends, how="inner").dropna(subset=["slope_ft_yr"])
    coll = pd.read_csv(os.path.join(DATA, "law_collapse.csv"))

    rows = []
    for aq, g in w.groupby(w.national_aquifer_code.fillna("UNKNOWN")):
        b = g.slope_ft_yr.values
        N = len(b)
        if N < 100:
            continue
        mu, sd = b.mean(), b.std(ddof=0)
        if sd == 0:
            continue
        snr_mean = abs(mu) / sd
        sign = 1 if mu > 0 else -1
        # benchmark: exact-by-simulation mean-functional flips (1,000 reps)
        bench = {}
        for f in RET:
            n = max(10, int(f * N))
            if n >= N:
                continue
            flips = 0
            for _ in range(1000):
                sub = RNG.choice(b, size=n, replace=False)
                if (1 if sub.mean() > 0 else -1) != sign:
                    flips += 1
            bench[n] = flips / 1000
        # actual experiment curves for this stratum
        act = coll[(coll.stratum == aq)]
        if not len(act):
            continue
        # c_i: ratio of x-values at matched flip ~0.10 (or nearest)
        def x_at(curve, target=0.10):
            xs = np.sqrt(curve.n.values * N / (N - np.minimum(curve.n.values, N - 1)))
            flips = curve.flip.values
            ok = np.isfinite(flips) & (flips > 0)
            if ok.sum() < 3:
                return np.nan
            # probit of flip vs snr*x -> slope = c_i
            from scipy.optimize import minimize_scalar
            xx = curve.snr.values * xs
            m = ok
            def nll(cc):
                p = np.clip(norm.cdf(-cc * xx[m]), 1e-9, 1 - 1e-9)
                return -np.sum(flips[m] * np.log(p) + (1 - flips[m]) * np.log(1 - p))
            r = minimize_scalar(nll, bounds=(0.2, 8), method="bounded")
            return r.x
        c_actual = x_at(act)
        # benchmark c (should be ~1)
        bd = pd.DataFrame({"n": list(bench.keys()), "flip": list(bench.values())})
        bd["snr"] = snr_mean
        c_bench = x_at(bd)
        rows.append(dict(stratum=aq, N=N, mean_slope=mu, sd=sd,
                         snr_mean=snr_mean, c_i=c_actual, c_bench=c_bench,
                         functional_gain=(c_bench / c_actual if c_actual and c_bench else np.nan)))
        print(f"{aq:12s} N={N:5d} SNR_mean={snr_mean:6.3f}  c_i(actual)={c_actual:5.2f}  "
              f"c(benchmark)={c_bench:5.2f}  bench-c>1 means robust functional helps: "
              f"{'yes' if c_bench and c_bench > 1.15 else 'no'}")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DATA, "decomposition.csv"), index=False)
    summ = dict(n_strata=len(df),
                c_i_range=[round(float(df.c_i.min()), 2), round(float(df.c_i.max()), 2)],
                c_bench_median=round(float(df.c_bench.median()), 3),
                functional_gain_median=round(float(df.functional_gain.median()), 3),
                corr_c_i_glb_absSpearman=None)
    json.dump(summ, open(os.path.join(DATA, "decomposition_summary.json"), "w"), indent=1)
    print("\nsummary:", summ)


if __name__ == "__main__":
    main()
