#!/usr/bin/env python3
"""Build the analysis dataset from downloaded OGC API pages.

- Loads all field-measurement pages (parameter 72019, depth to water below
  land surface, ft) and all monitoring-location pages.
- Joins site metadata (state, national aquifer, well depth) on
  monitoring_location_id.
- Applies QA filters and selects wells eligible for trend analysis
  (>= MIN_YEARS distinct calendar years, >= MIN_MEAS measurements).
- Outputs: data/well_measurements.csv (eligible long-format series),
  data/eligible_wells.csv (site attributes), data/qa_report.json.
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))

MIN_MEAS = 30          # minimum discrete measurements per well
MIN_YEARS = 20         # minimum distinct calendar years per well
MIN_VALUE, MAX_VALUE = 0.0, 5000.0   # plausible depth-to-water range (ft)

MEAS_COLS = ["x", "y", "monitoring_location_id", "value", "time",
             "qualifier", "approval_status", "year", "month"]
SITE_COLS = ["id", "x", "y", "state_code", "state_name", "county_name",
             "national_aquifer_code", "aquifer_code", "site_type_code",
             "well_constructed_depth", "altitude", "construction_date"]


def load_pages(tag, cols):
    """Load pages from one dir (tag='sites') or the union of shard dirs
    (tag='gw' matches pages_gw_*)."""
    dirs = sorted(glob.glob(os.path.join(DATA, "pages_%s_*" % tag))) or \
        [os.path.join(DATA, "pages_%s" % tag)]
    files = []
    for d in dirs:
        files += sorted(glob.glob(os.path.join(d, "page_*.csv")))
    if not files:
        sys.exit("no pages found for %s" % tag)
    frames = []
    for i, f in enumerate(files):
        try:
            df = pd.read_csv(f, usecols=lambda c: c in set(cols), low_memory=False)
            frames.append(df)
        except Exception as e:  # noqa: BLE001
            print("WARN: skipping %s (%s)" % (f, e))
        if (i + 1) % 20 == 0:
            print("  %s: %d/%d pages" % (tag, i + 1, len(files)), flush=True)
    out = pd.concat(frames, ignore_index=True)
    return out


def main():
    qa = {}

    meas = load_pages("gw", MEAS_COLS)
    qa["raw_measurement_rows"] = len(meas)
    qa["wells_with_any_record"] = int(
        meas["monitoring_location_id"].astype(str).nunique())
    print("raw measurement rows:", len(meas),
          "at", qa["wells_with_any_record"], "wells", flush=True)

    meas["monitoring_location_id"] = meas["monitoring_location_id"].astype(str)
    meas["time"] = pd.to_datetime(meas["time"], errors="coerce")
    meas["value"] = pd.to_numeric(meas["value"], errors="coerce")
    meas["year"] = meas["time"].dt.year
    meas["month"] = meas["time"].dt.month

    # --- QA filters -------------------------------------------------------
    n0 = len(meas)
    meas = meas.dropna(subset=["time", "value"])
    qa["drop_invalid_time_or_value"] = n0 - len(meas)

    n1 = len(meas)
    in_range = meas["value"].between(MIN_VALUE, MAX_VALUE)
    qa["drop_value_out_of_range"] = int((~in_range).sum())
    meas = meas[in_range]

    # approved/validated records only, if approval status present
    if "approval_status" in meas.columns:
        n2 = len(meas)
        meas = meas[meas["approval_status"].fillna("").str.contains("Approved", case=False)]
        qa["drop_not_approved"] = n2 - len(meas)

    # drop non-ambient (pumping/flowing) levels and dry-well flags
    qcol = "qualifier" if "qualifier" in meas.columns else None
    if qcol:
        flag = meas[qcol].fillna("").str.lower()
        bad = flag.str.contains("pump", na=False) | flag.str.contains("dry", na=False)
        qa["drop_pumping_level"] = int(bad.sum())
        meas = meas[~bad]

    # one record per well-day (median of same-day repeats)
    meas["date"] = meas["time"].dt.date
    meas = (meas.groupby(["monitoring_location_id", "date"], as_index=False)
                .agg(value=("value", "median"), month=("month", "first"),
                     year=("year", "first"), x=("x", "first"), y=("y", "first")))
    meas["time"] = pd.to_datetime(meas["date"])
    qa["well_day_records"] = len(meas)

    # CONUS scope: drop Alaska, Hawaii, Puerto Rico & other territories
    conus = meas["y"].between(24, 50) & meas["x"].between(-125, -66)
    qa["drop_outside_CONUS_well_days"] = int((~conus).sum())
    meas = meas[conus]

    # --- site metadata join ----------------------------------------------
    sites = load_pages("sites", SITE_COLS)
    sites["id"] = sites["id"].astype(str)
    sites = sites.drop_duplicates(subset="id")
    qa["site_rows"] = len(sites)

    meas = meas.merge(sites[SITE_COLS].rename(columns={"id": "monitoring_location_id"}),
                      on="monitoring_location_id", how="left")

    # --- eligible wells ---------------------------------------------------
    g = meas.groupby("monitoring_location_id")
    stats = g.agg(n_meas=("value", "size"), n_years=("year", "nunique"),
                  y0=("year", "min"), y1=("year", "max"),
                  n_months=("month", "nunique"))
    elig = stats[(stats["n_meas"] >= MIN_MEAS) & (stats["n_years"] >= MIN_YEARS)].copy()
    elig["span_years"] = elig["y1"] - elig["y0"] + 1
    qa["eligible_wells"] = len(elig)
    qa["eligibility_rule"] = dict(min_measurements=MIN_MEAS, min_years=MIN_YEARS)

    meas = meas[meas["monitoring_location_id"].isin(elig.index)].copy()

    out_m = os.path.join(DATA, "well_measurements.csv")
    meas.to_csv(out_m, index=False)
    site_meta = sites.set_index("id").drop(columns=["id"], errors="ignore")
    elig.join(site_meta, how="left").to_csv(os.path.join(DATA, "eligible_wells.csv"))
    with open(os.path.join(DATA, "qa_report.json"), "w") as fh:
        json.dump(qa, fh, indent=2)
    print(json.dumps(qa, indent=2))
    print("wrote", out_m)


if __name__ == "__main__":
    main()
