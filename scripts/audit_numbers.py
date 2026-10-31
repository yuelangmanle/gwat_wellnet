#!/usr/bin/env python3
"""Submission audit: recompute every number claimed in the manuscript from
the data files and report agreement. Also verifies strata counts, state
coverage, and unit conversions."""
import json
import os
import re

import numpy as np
import pandas as pd

DATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
MS = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "manuscript", "manuscript_draft.md"))
FT = 0.3048

res = pd.read_csv(os.path.join(DATA, "experiment_results.csv"))
th = pd.read_csv(os.path.join(DATA, "thresholds.csv"))
trends = pd.read_csv(os.path.join(DATA, "well_trends.csv"), dtype={"monitoring_location_id": str})
wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"), dtype={"monitoring_location_id": str})
qa = json.load(open(os.path.join(DATA, "qa_report.json")))
sens = pd.read_csv(os.path.join(DATA, "sensitivity_aggregation.csv"))
strat = pd.read_csv(os.path.join(DATA, "strategy_results.csv"))
w = wells.merge(trends, on="monitoring_location_id")

if not os.path.exists(MS):
    raise SystemExit("manuscript_draft.md not found - the number audit "
                     "checks the manuscript text against the data, so it "
                     "needs the manuscript source in place.")
ms = open(MS).read()
checks = []


def chk(name, claimed, actual, tol=0.0):
    if isinstance(claimed, float) or isinstance(actual, float):
        ok = abs(float(claimed) - float(actual)) <= tol
    else:
        ok = claimed == actual
    checks.append((name, claimed, actual, "PASS" if ok else "FAIL"))
    return ok


# --- data funnel ---
chk("raw measurements", 9456561, qa["raw_measurement_rows"])
chk("wells with any record (843,302 claim)", 843302, qa["wells_with_any_record"])
chk("well-day records", 8814831, qa["well_day_records"])
chk("eligible wells", 29717, qa["eligible_wells"])
chk("pumping+dry dropped (234,752)", 234752, qa["drop_pumping_level"])
chk("outside CONUS dropped", 103049, qa["drop_outside_CONUS_well_days"])
chk("eligibility fraction 3.5%", 3.5, round(100 * 29717 / 843302, 1), 0.05)

# --- well trends ---
chk("n wells with trends", 29717, len(trends))
chk("falling wells 18,513", 18513, int((trends.slope_ft_yr > 0).sum()))
chk("rising wells 11,134 (stable 70 excluded)", 11134, int((trends.slope_ft_yr < 0).sum()))
chk("stable wells 70", 70, int((trends.slope_ft_yr == 0).sum()))
chk("falling + rising + stable partition = 29,717", 29717,
    int((trends.slope_ft_yr > 0).sum() + (trends.slope_ft_yr < 0).sum() + (trends.slope_ft_yr == 0).sum()))
chk("share falling 62.3%", 62.3, round(100 * (trends.slope_ft_yr > 0).mean(), 1), 0.05)
chk("significant 70.1%", 70.1, round(100 * (trends.p_value < 0.05).mean(), 1), 0.05)
chk("AR-corrected 60.2%", 60.2, round(100 * trends.ar_corrected.mean(), 1), 0.05)
chk("median slope m/yr 0.010", 0.010, round(trends.slope_ft_yr.median() * FT, 3), 0.0005)
q1, q3 = trends.slope_ft_yr.quantile([0.25, 0.75])
chk("IQR low m/yr -0.010", -0.010, round(q1 * FT, 3), 0.0005)
chk("IQR high m/yr 0.083", 0.083, round(q3 * FT, 3), 0.0005)
chk("min slope m/yr -9.8", -9.8, round(trends.slope_ft_yr.min() * FT, 1), 0.05)
chk("max slope m/yr 5.4", 5.4, round(trends.slope_ft_yr.max() * FT, 1), 0.05)
beyond = (trends.slope_ft_yr.abs() > 5).mean() * 100
chk("beyond +/-1.5 m/yr share (0.3%)", 0.3, round(beyond, 1), 0.05)

# --- strata counts ---
aq = w.national_aquifer_code.fillna("UNKNOWN")
counts = aq.value_counts()
strata_all = sorted(counts[counts >= 100].index.tolist())
assigned = [s for s in strata_all if s != "UNKNOWN"]
chk("aquifer strata >=100 wells (incl UNKNOWN, ms says 21)", 21, len(strata_all))
chk("assigned strata (excl UNKNOWN)", 20, len(assigned))
chk("UNKNOWN wells 6,932", 6932, int(counts.get("UNKNOWN", 0)))
chk("HGHPLN wells 8,160", 8160, int(counts.get("N100HGHPLN", 0)))
chk("NATLCP wells 1,959", 1959, int(counts.get("S100NATLCP", 0)))
chk("CSLLWD wells 1,145", 1145, int(counts.get("S100CSLLWD", 0)))
# states
states = w.state_code.dropna().astype(float).astype(int)
conus_states = set(states.unique())
chk("states+DC covered", 49, len(conus_states))

# --- reference slopes (m/yr) ---
def ref_slope(code):
    g = res[(res.stratum == code) & (~res.clustered)]
    return float(g.ref_slope.iloc[0]) if len(g) else np.nan
chk("CSLLWD ref m/yr -0.17", -0.17, round(ref_slope("S100CSLLWD") * FT, 2), 0.005)
chk("CACSTL ref m/yr 0.07", 0.07, round(ref_slope("N100CACSTL") * FT, 2), 0.005)
chk("BSNRGB ref m/yr 0.002", 0.002, round(ref_slope("N100BSNRGB") * FT, 3), 0.0005)
# CNRLVL fastest individual decline
cnr = trends[trends.monitoring_location_id.isin(w[w.national_aquifer_code == "S100CNRLVL"].monitoring_location_id)]
chk("CNRLVL fastest falling ~1.4 m/yr", 1.4, round(float((cnr.slope_ft_yr * FT).max()), 1), 0.05)

# --- CONUS experiment ---
c = res[(res.stratum == "CONUS") & (~res.clustered)]
chk("CONUS flip@5% = 0.000", 0.000, float(c[np.isclose(c.retention, 0.05)].p_sign_flip.iloc[0]), 0.001)
chk("CONUS flip@1% = 0.075", 0.075, float(c[np.isclose(c.retention, 0.01)].p_sign_flip.iloc[0]), 0.002)
chk("CONUS n@1% = 297", 297, int(c[np.isclose(c.retention, 0.01)].n_keep.iloc[0]))
chk("CONUS n@5% = 1,485", 1485, int(c[np.isclose(c.retention, 0.05)].n_keep.iloc[0]))
cl = res[(res.stratum == "CONUS") & (res.clustered)]
chk("CONUS clustered flip@30% = 0.10", 0.10, round(float(cl[np.isclose(cl.retention, 0.30)].p_sign_flip.iloc[0]), 2), 0.005)
chk("CONUS clustered n@30% = 7,726", 7726, int(cl[np.isclose(cl.retention, 0.30)].n_keep.iloc[0]))
chk("CONUS clustered n_all 25,754", 25754, int(cl.n_all.iloc[0]))

# --- thresholds ---
def nmin(code):
    g = th[(th.stratum == code) & (~th.clustered)]
    return float(g.n_min.iloc[0]) if len(g) and pd.notna(g.n_min.iloc[0]) else None
for code, claim in [("S100CSLLWD", 11), ("N100MSRVVL", 12), ("N100HGHPLN", 81),
                    ("S100NATLCP", 39), ("N100CACSTL", 25), ("N100GLCIAL", 71),
                    ("S500INTRMD", 10), ("N300COPLTS", 10), ("S100RIOGRD", 519),
                    ("S100CNRLVL", 146)]:
    chk(f"n_min {code} = {claim}", claim, nmin(code))
chk("n_min RIOGRD frac ~70%", 0.699, round(nmin("S100RIOGRD") / 742, 3), 0.01)
chk("n_min CNRLVL frac ~80%", 0.798, round(nmin("S100CNRLVL") / 183, 3), 0.01)
never = sorted(th[(~th.clustered) & th.n_min.isna()].stratum.tolist())
chk("strata with no n_min (3: BSNRGB, SURFCL, FLORDN)", 3, len(never))

# --- sampling law ---
lawfit = json.load(open(os.path.join(DATA, "law_fit.json")))
nminlaw = pd.read_csv(os.path.join(DATA, "law_nmin.csv"))
collaw = pd.read_csv(os.path.join(DATA, "law_collapse.csv"))
chk("law c = 1.27", 1.27, round(lawfit["c"], 2), 0.005)
chk("law CI 1.04-1.55", True, (round(lawfit["ci"][0], 2), round(lawfit["ci"][1], 2)) == (1.04, 1.55))
chk("collapse points 216", 216, len(collaw))
chk("17/19 within factor 2", 17, int(((nminlaw[~nminlaw.censored].nmin_sim /
                                       nminlaw[~nminlaw.censored].nmin_pred).between(0.5, 2)).sum()))
uncens = nminlaw[~nminlaw.censored]
chk("uncensored strata 19", 19, len(uncens))

# --- grid-floor censoring and criterion granularity (as disclosed in text) ---
fl = nminlaw.dropna(subset=["nmin_sim"]).copy()
fl["floor"] = fl.N.apply(lambda n: max(10, int(0.01 * n)))
onf = fl[fl.nmin_sim == fl.floor]
claimed_floor = {"CONUS", "N100HGHPLN", "N100MSRVVL", "N300COPLTS", "S100CSLLWD",
                 "S500INTRMD", "UNKNOWN"}
chk("floor-censored minima 7 strata", 7, len(onf))
chk("derived floor set matches manuscript list", True,
    set(onf.stratum) == claimed_floor)
chk("no minimum below the 1% floor", 0, int((fl.nmin_sim < fl.floor).sum()))
nconf = pd.read_csv(os.path.join(DATA, "nmin_confidence.csv")).merge(
    nminlaw[["stratum", "nmin_pred"]], on="stratum")
chk("pred within sim binomial CI 4/19", 4,
    int(((nconf.nmin_pred >= nconf.nmin_lo) & (nconf.nmin_pred <= nconf.nmin_hi)).sum()))
lc = pd.read_csv(os.path.join(DATA, "law_collapse.csv"))
chk("trials with observed flip>0.5", 7, int((lc.flip > 0.5).sum()))
chk("band coverage at flip>0.5 is 0/7", 0, int(lc[lc.flip > 0.5].in_ci.sum()))
chk("window-start sign changes 4/21", 4,
    int(pd.read_csv(os.path.join(DATA, "window_sensitivity.csv")).changes.sum()))
grp = wells.groupby(["x", "y"]).size()
chk("co-located wells 3064", 3064, int(grp[grp > 1].sum()))
aan = pd.read_csv(os.path.join(DATA, "annual_anoms.csv"),
                  usecols=["monitoring_location_id", "year"])
_wy = aan[(aan.year >= 1975) & (aan.year <= 2025)].groupby(
    "monitoring_location_id").year.nunique().reindex(
    wells.monitoring_location_id, fill_value=0)
chk("median window-year coverage 29 of 51", 29, int(_wy.median()))
chk("wells covering all 51 window years 3.5%", 3.5,
    round(100 * (_wy == 51).mean(), 1), 0.05)
_g = sorted(res.retention.unique())[::-1]
_rr = []
for _st in nminlaw.stratum:
    _N = int(nminlaw[nminlaw.stratum == _st].N.iloc[0])
    _nk = [max(10, int(np.floor(r * _N))) for r in _g]
    _rr += [a / b for a, b in zip(_nk[:-1], _nk[1:]) if b > 0]
chk("adjacent grid sizes differ 1.5x median", 1.5, round(float(np.median(_rr)), 1), 0.05)
cen = nminlaw[nminlaw.censored].set_index("stratum").ret_pred
chk("censored predictions 94/91/88%", True,
    (round(cen["N100BSNRGB"], 2), round(cen["S100SURFCL"], 2), round(cen["S400FLORDN"], 2)) == (0.94, 0.91, 0.88))
cam = collaw[(collaw.stratum == "S300CAMORD") & (collaw.n == 30)].iloc[0]
chk("CAMORD obs flip@30% 0.29", 0.29, round(cam.flip, 2), 0.005)
chk("CAMORD law pred 0.25", 0.25, round(cam.pred, 2), 0.005)
chk("clustered n_eff cap 184-281", True,
    round(lawfit["clustered"]["neff_cap_10_30"][0]) == 184 and
    round(lawfit["clustered"]["neff_cap_10_30"][1]) == 281)
chk("HGHPLN law n_min 82 vs sim 81", True,
    int(nminlaw.set_index('stratum').loc['N100HGHPLN', 'nmin_pred']) == 82)
# manuscript formula present
chk("ms has closed-form n_min", True, "n_min = Nz²/(z²+(c·SNR)²N)" in ms)

# --- aggregation sensitivity ---
chk("sensitivity rows (21 aq + CONUS)", 22, len(sens))
chk("sign agreement count 20", 20, int(sens.sign_agree.sum()))
dis = sorted(sens[~sens.sign_agree].stratum.tolist())
chk("disagreeing = FLORDN, RIOGRD", sorted(["S400FLORDN", "S100RIOGRD"]), sorted(dis))

# --- strategies ---
g = strat[strat.stratum == "N100HGHPLN"]
lf = g[g.strategy == "longest_first"]
nr = g[g.strategy == "nested_random"]
r10lf = float(lf[np.isclose(lf.retention, 0.10, atol=0.006)].abs_bias.iloc[0])
r10nr = float(nr[np.isclose(nr.retention, 0.10, atol=0.006)].abs_bias.iloc[0])
chk("HP longest bias@10% m/yr 0.017", 0.017, round(r10lf * FT, 3), 0.001)
chk("HP band bias@10% m/yr 0.003", 0.003, round(r10nr * FT, 3), 0.0007)
chk("ratio 5.4x", 5.4, round(r10lf / r10nr, 1), 0.1)

# --- CV spatial-spread bias vs reference slope (3.1x) ---
cv = strat[(strat.stratum == "S100CNRLVL") & (strat.strategy == "spatial_spread")]
cvmax = float(cv.abs_bias.max())
cvref = float(cv.ref_slope.iloc[0])
chk("CV spatial max bias ratio 3.1x", 3.1, round(cvmax / abs(cvref), 1), 0.1)

# --- abstract count claims ---
chk("abstract names three unstable systems", True,
    "(Basin and Range, surficial, Floridan)" in ms)
chk("abstract says 88-94%", True, "88\u201394%" in ms)
ctj = json.load(open(os.path.join(DATA, "coverage_tests.json")))
chk("low-flip band coverage 74%", 0.74, round(ctj["coverage_low_flip"], 2), 0.005)
chk("transition band coverage 23%", 0.23, round(ctj["coverage_transition"], 2), 0.005)
chk("ms cites 74% low-flip and 23% transition", True,
    "74% of trials" in ms and "23% coverage" in ms)

# --- synthetic testbed + Kish forward (sim_validation.json) ---
sv = json.load(open(os.path.join(DATA, "sim_validation.json")))
cn = [v["c_app"] for k, v in sv["cells"].items() if not k.endswith("noise")]
cy = [v["c_app"] for k, v in sv["cells"].items() if k.endswith("noise")]
chk("testbed noiseless c 1.01-1.04", "1.01\u20131.04", f"{min(cn):.2f}\u2013{max(cn):.2f}", 0.005)
chk("testbed noise c 0.91-0.95", "0.91\u20130.95", f"{min(cy):.2f}\u2013{max(cy):.2f}", 0.005)
kf = [p for p in sv["kish_forward"]["points"] if abs(p["ret"] - 0.2) < 0.01 or abs(p["ret"] - 0.3) < 0.01]
dev = max(abs(p["obs"] - p["pred"]) for p in kf)
chk("Kish forward 20-30% band within ~1pt", 0.01, round(dev, 3), 0.005)

dec = pd.read_csv(os.path.join(DATA, "decomposition.csv"))
civ = dec.c_i.dropna()
chk("c_i range 0.55-4.05", "0.55\u20134.05", f"{civ.min():.2f}\u2013{civ.max():.2f}", 0.005)

print(f"{'check':58s} {'claim':>12s} {'actual':>12s}  verdict")
fails = 0
for name, claimed, actual, verdict in checks:
    if verdict == "FAIL":
        fails += 1
    print(f"{name:58s} {str(claimed):>12s} {str(actual):>12s}  {verdict}")
print(f"\n{len(checks)-fails}/{len(checks)} PASS, {fails} FAIL")
