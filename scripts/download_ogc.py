#!/usr/bin/env python3
"""Bulk-download USGS water data OGC API collections via cursor paging.

Each page is fetched as JSON (the cursor lives in the same response's next
link), flattened to CSV locally, and the cursor checkpointed after every page
so the job can resume after a crash. Retries transient server-side timeouts
with exponential backoff.
"""
import argparse
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request

BASE = "https://api.waterdata.usgs.gov/ogcapi/v1/collections"
PAGE_DIR = os.path.dirname(os.path.abspath(__file__))


def fetch_json(url, timeout=300, max_tries=30):
    import urllib.error
    last = None
    for attempt in range(1, max_tries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "gwat-wellnet-research/0.1"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if r.getcode() == 200:
                    body = r.read()
                    if body.strip():
                        return json.loads(body)
                    raise IOError("empty body")
                raise IOError("HTTP %s" % r.getcode())
        except urllib.error.HTTPError as e:
            last = e
            if e.code == 429:
                retry_after = e.headers.get("retry-after") if e.headers else None
                wait = int(retry_after) + 15 if retry_after else 600
                sys.stderr.write("429 rate-limited; honoring retry-after %ss\n" % wait)
                time.sleep(wait)
            else:
                wait = min(2 ** attempt, 90)
                sys.stderr.write("attempt %d failed (%s); sleep %ds\n" % (attempt, e, wait))
                time.sleep(wait)
        except Exception as e:  # noqa: BLE001
            last = e
            wait = min(2 ** attempt, 90)
            sys.stderr.write("attempt %d failed (%s); sleep %ds\n" % (attempt, e, wait))
            time.sleep(wait)
    raise IOError("giving up on %s: %s" % (url, last))


def flatten(doc):
    """Yield dicts for each feature: x, y + all properties."""
    rows = []
    for feat in doc.get("features", []):
        props = dict(feat.get("properties") or {})
        geom = feat.get("geometry") or {}
        coords = geom.get("coordinates") or [None, None]
        if "x" not in props:
            props["x"] = coords[0] if len(coords) > 0 else None
        if "y" not in props:
            props["y"] = coords[1] if len(coords) > 0 else None
        rows.append(props)
    return rows


def write_page(path, rows):
    if not rows:
        return 0
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return len(rows)


def download(collection, params, name, limit=50000, sleep_s=4.0):
    outdir = os.path.abspath(os.path.join(PAGE_DIR, "..", "data", "pages_%s" % name))
    os.makedirs(outdir, exist_ok=True)
    state_path = os.path.join(outdir, "_cursor.json")

    cursor, page, total_rows, done = None, 0, 0, False
    if os.path.exists(state_path):
        with open(state_path) as fh:
            st = json.load(fh)
        cursor, page, done = st.get("cursor"), st.get("page", 0), st.get("done", False)
        for fn in sorted(os.listdir(outdir)):
            if fn.endswith(".csv"):
                with open(os.path.join(outdir, fn), errors="ignore") as fh:
                    total_rows += max(0, sum(1 for _ in fh) - 1)
        if done:
            print("%s already complete (%d rows, %d pages)" % (name, total_rows, page), flush=True)
            return
        if cursor is None:
            page = 0  # no cursor saved -> start fresh
        print("resuming %s at page %d, cursor %s" % (name, page, cursor), flush=True)

    while True:
        q = dict(params)
        q["f"] = "json"
        q["limit"] = str(limit)
        url = "%s/%s/items?%s" % (BASE, collection, urllib.parse.urlencode(q))
        if cursor:
            url += "&cursor=" + urllib.parse.quote(cursor)
        doc = fetch_json(url)

        rows = flatten(doc)
        page += 1
        page_file = os.path.join(outdir, "page_%06d.csv" % page)
        n = write_page(page_file, rows)
        total_rows += n
        print("page %d: %d rows (cum %d)" % (page, n, total_rows), flush=True)

        nxt = None
        for l in doc.get("links", []):
            if l.get("rel") == "next":
                nxt = l.get("href")
        if not nxt or n == 0:
            with open(state_path, "w") as fh:
                json.dump({"cursor": None, "page": page, "done": True}, fh)
            print("%s complete: %d rows across %d pages" % (name, total_rows, page), flush=True)
            return
        qs = urllib.parse.parse_qs(urllib.parse.urlparse(nxt).query)
        cursor = qs.get("cursor", [None])[0]
        with open(state_path, "w") as fh:
            json.dump({"cursor": cursor, "page": page}, fh)
        time.sleep(sleep_s)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--collection", required=True)
    ap.add_argument("--params", default="", help="key=value&key2=value2 extra filters")
    ap.add_argument("--name", required=True, help="short name for page directory")
    ap.add_argument("--limit", type=int, default=50000)
    args = ap.parse_args()
    extra = {}
    if args.params:
        for kv in args.params.split("&"):
            k, _, v = kv.partition("=")
            extra[k] = v
    download(args.collection, extra, args.name, limit=args.limit)
