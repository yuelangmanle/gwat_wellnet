#!/usr/bin/env python3
"""Main figures (self-explanatory design).

Figure contracts (core conclusion -> panels), each panel carrying its own
finding title so the figure reads without the caption:
- Fig 1 (two-column): a where the eligible wells are (trend-sign map with
  state anchors); b minimum reliable size per stratum: simulated dot vs
  law diamond against the full network, regime colours, censored systems
  labelled with the law's predicted fraction.
- Fig 2 (two-column): a six largest aquifer strata: flip risk vs size with
  risk zones and crossing beads; b the CONUS network: random vs clustered
  loss with the law curve and the effective-size reading.
- Fig 3 (two-column): a every trial collapses onto the sampling law
  (regime-coloured, theory curve emphasised); b the design chart: named
  networks placed at (SNR, N), isopleths of required retention.
- figS (Supporting Information): retention strategy moves bias, not sign risk.
"""
import json
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))
FIGS = os.path.abspath(os.path.join(HERE, "..", "figures"))
os.makedirs(FIGS, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    # math ($n$, $\Phi$...) from the same face as the surrounding text,
    # so exports embed one font family instead of Arial + DejaVu mix
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "font.size": 7,
    "axes.labelsize": 7.5,
    "axes.titlesize": 8,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 1.0,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "legend.frameon": False,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
})

TWO_COL = 17.15 / 2.54   # 6.752 in
ONE_COL = 8.25 / 2.54    # 3.248 in
# Groundwater figure guideline: lines >= 1 pt, gray fills 20-70 % black,
# no type below 5 pt, TIF/EPS >= 300 dpi at 8.25 cm (one column) or
# 17.15 cm (two columns).
MIN_LW = 1.0
GRAY_LIGHT = "#C4C8CE"   # ~22 % black: lightest gray allowed

# restrained palette: neutral family + signal family + accent
C_NEUTRAL = "#8A8F98"
C_SIGNAL = "#2166AC"     # rising water level
C_SIGNAL2 = "#B2182B"    # declining water level
C_ACCENT = "#E08214"
C_DARK = "#252525"
C_SAFE = "#228833"
C_WARN = "#E08214"
C_RISK = "#B2182B"

FIT = json.load(open(os.path.join(DATA, "law_fit.json")))
C_LAW = FIT["c"]
Z = FIT["z"]

_aq_names_path = os.path.join(DATA, "aquifer_names.json")
AQUIFER_NAMES = json.load(open(_aq_names_path)) if os.path.exists(_aq_names_path) else {}
SHORT = {
    "N100HGHPLN": "High Plains",
    "S100NATLCP": "N. Atlantic\nCoastal Plain",
    "N100BSNRGB": "Basin & Range",
    "N100GLCIAL": "glacial aquifers",
    "N100MSRVVL": "MS River\nalluvium",
    "S100CSLLWD": "coastal lowlands",
    "S100CNRLVL": "Central Valley",
    "S100RIOGRD": "Rio Grande",
    "S300CAMORD": "Cambrian–Ordov.",
    "N600CMBPLV": "Columbia Plateau",
    "S400RSWLBS": "Roswell Basin",
    "N400BISCYN": "Biscayne",
    "S100MSEMBM": "MS embayment",
    "S100SURFCL": "surficial system",
    "S500INTRMD": "Intermediate",
    "N300COPLTS": "Colorado Plateaus",
    "N100CACSTL": "California coastal",
    "N100ALLUVL": "alluvial aquifers",
    "S400FLORDN": "Floridan",
    "N9999OTHER": "other aquifers",
    "UNKNOWN": "unassigned",
}
QUAL = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#AA3377", "#66CCEE"]


def aq_name(code):
    return SHORT.get(code) or AQUIFER_NAMES.get(code, code)


def regime_color(ret):
    """Shared regime palette across figures: required fraction of network."""
    return C_SAFE if ret <= 0.05 else (C_WARN if ret <= 0.5 else C_RISK)


def enforce_min_linewidth(fig, min_lw=MIN_LW):
    """Raise every drawn stroke (spines, ticks, lines, patch edges, arrows,
    legend handles) to the journal minimum; zero-width strokes stay off."""
    def bump(artist):
        lw = artist.get_linewidth()
        if lw is not None and 0 < lw < min_lw:
            artist.set_linewidth(min_lw)
    for ax in fig.get_axes():
        for sp in ax.spines.values():
            bump(sp)
        ax.tick_params(which="both", width=min_lw)
        for ln in ax.get_lines():
            bump(ln)
        for p in ax.patches:
            bump(p)
        for t in ax.texts:
            ap = getattr(t, "arrow_patch", None)
            if ap is not None:
                bump(ap)
        leg = ax.get_legend()
        if leg is not None:
            handles = getattr(leg, "legend_handles", None) or \
                getattr(leg, "legendHandles", [])
            for h in handles:
                if isinstance(h, Line2D):
                    bump(h)


def normalize_raster_width(path, width_in, compression=None):
    """Tight bounding boxes grow or shrink the canvas by a few millimetres;
    reset the resolution tag so the file's physical width is exactly the
    column width the journal specifies (pixels untouched, dpi stays >= 300)."""
    from PIL import Image
    im = Image.open(path)
    if compression and im.mode == "RGBA":
        # production TIFFs: flatten the (opaque) alpha channel onto white
        flat = Image.new("RGB", im.size, "white")
        flat.paste(im, mask=im.getchannel("A"))
        im = flat
    px_w, px_h = im.size
    dpi = px_w / width_in
    kw = {"dpi": (dpi, dpi)}
    if compression:
        kw["compression"] = compression
    im.save(path, **kw)
    assert dpi >= 300, (path, dpi)
    return px_w, px_h, dpi


def save_pub(fig, name, width_in):
    fig.set_size_inches(width_in, fig.get_size_inches()[1])
    enforce_min_linewidth(fig)
    base = os.path.join(FIGS, name)
    fig.savefig(base + ".svg", bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".png", dpi=450, bbox_inches="tight")
    fig.savefig(base + ".tiff", dpi=450, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    normalize_raster_width(base + ".png", width_in)
    w, h, dpi = normalize_raster_width(base + ".tiff", width_in, "tiff_lzw")
    print("saved %s  %dx%d px  %.2f x %.2f cm  %.0f dpi"
          % (base, w, h, width_in * 2.54, h / dpi * 2.54, dpi))


def load():
    res = pd.read_csv(os.path.join(DATA, "experiment_results.csv"))
    strat = pd.read_csv(os.path.join(DATA, "strategy_results.csv"))
    wells = pd.read_csv(os.path.join(DATA, "eligible_wells.csv"),
                        dtype={"monitoring_location_id": str})
    trends = pd.read_csv(os.path.join(DATA, "well_trends.csv"),
                         dtype={"monitoring_location_id": str})
    wells = wells.merge(trends, on="monitoring_location_id", how="inner")
    law = pd.read_csv(os.path.join(DATA, "law_nmin.csv"))
    coll = pd.read_csv(os.path.join(DATA, "law_collapse.csv"))
    return res, strat, wells, law, coll


def risk_zones(ax, ymax, key=""):
    """Shaded flip-risk bands with right-edge labels."""
    ax.axhspan(0, 0.10, color=C_SAFE, alpha=0.06, lw=0, zorder=0)
    if ymax > 0.12:
        ax.axhspan(0.10, 0.30, color=C_WARN, alpha=0.07, lw=0, zorder=0)
        ax.axhspan(0.30, ymax, color=C_RISK, alpha=0.05, lw=0, zorder=0)
    ax.axhline(0.10, color=C_NEUTRAL, lw=1.0, ls=(0, (3, 3)), zorder=1)


# ------------------------------------------------- collision-free placement
def _segments_cross_rect(x0, y0, x1, y1, xs, ys):
    """True if any polyline segment's extent overlaps the rect (data coords)."""
    sx0, sx1 = np.minimum(xs[:-1], xs[1:]), np.maximum(xs[:-1], xs[1:])
    sy0, sy1 = np.minimum(ys[:-1], ys[1:]), np.maximum(ys[:-1], ys[1:])
    over_x = (sx0 < x1) & (sx1 > x0)
    over_y = (sy0 < y1) & (sy1 > y0)
    return bool((over_x & over_y).any())


def place_clear(ax, txt, candidates, fontsize, color, weight="normal",
                pad=1.5, **kw):
    """Place txt at the first candidate (x, y, ha, va) whose rendered box
    crosses no Line2D segment in ax; fall back to fewest crossings."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lines = []
    for ln in ax.get_lines():
        if ln.get_linestyle() in (None, "None", "", " "):
            continue
        xs, ys = np.asarray(ln.get_xdata(), float), np.asarray(ln.get_ydata(), float)
        m = np.isfinite(xs) & np.isfinite(ys)
        if m.sum() < 2:
            continue
        lines.append((xs[m], ys[m]))
    inv = ax.transData.inverted()
    best, best_n = None, None
    for cand in candidates:
        x, y, ha, va = cand
        t = ax.text(x, y, txt, fontsize=fontsize, color=color, ha=ha, va=va,
                    fontweight=weight, **kw)
        fig.canvas.draw()
        bb = t.get_window_extent(renderer)
        (dx0, dy0), (dx1, dy1) = inv.transform(
            [(bb.x0 - pad, bb.y0 - pad), (bb.x1 + pad, bb.y1 + pad)])
        n = sum(_segments_cross_rect(min(dx0, dx1), min(dy0, dy1),
                                     max(dx0, dx1), max(dy0, dy1), xs, ys)
                for xs, ys in lines)
        # also avoid other texts already placed in this axes
        for ot in ax.texts:
            if ot is t:
                continue
            obb = ot.get_window_extent(renderer)
            obb = obb.expanded(1.06, 1.25)
            if obb.overlaps(bb):
                n += 1
                break
        t.remove()
        if n == 0:
            ax.text(x, y, txt, fontsize=fontsize, color=color, ha=ha, va=va,
                    fontweight=weight, **kw)
            return True
        if best_n is None or n < best_n:
            best, best_n = cand, n
    x, y, ha, va = best
    ax.text(x, y, txt, fontsize=fontsize, color=color, ha=ha, va=va,
            fontweight=weight, **kw)
    return False


# ---------------------------------------------------------------- fig 1
def fig1(res):
    fig = plt.figure(figsize=(TWO_COL, 1.98))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1],
                          wspace=0.40, left=0.085, right=0.985, top=0.86,
                          bottom=0.185)
    # ---- panel a: six largest aquifer strata
    axb = fig.add_subplot(gs[0, 0])
    big = (res[(~res.clustered) & (~res.stratum.isin(["CONUS", "UNKNOWN"]))]
           .groupby("stratum").n_all.first().sort_values(ascending=False).head(6))
    for i, aq in enumerate(big.index):
        g = res[(res.stratum == aq) & (~res.clustered)].sort_values("n_keep")
        axb.plot(g.n_keep, g.p_sign_flip, "-", color=QUAL[i], lw=1.1, zorder=3)
    risk_zones(axb, 0.42)
    axb.set_ylim(0, 0.42)
    axb.set_xlim(8, 12000)
    axb.set_xscale("log")
    # crossing beads with numeric tags placed into the safe zone; names in a key
    beads = {"S100CSLLWD": "11", "N100MSRVVL": "12", "S100NATLCP": "39",
             "N100GLCIAL": "71", "N100HGHPLN": "81"}
    for aq, tag in beads.items():
        g = res[(res.stratum == aq) & (~res.clustered)].sort_values("n_keep")
        ok = g[g.p_sign_flip <= 0.10]
        if not len(ok):
            continue
        x = float(ok.n_keep.min())
        i = list(big.index).index(aq)
        axb.plot([x], [0.10], "o", ms=4.5, mfc=QUAL[i], mec="white",
                 mew=0.7, zorder=6)
        place_clear(axb, tag, [(x * 1.6, 0.072, "center", "center"),
                               (x * 1.6, 0.042, "center", "center"),
                               (x * 2.6, 0.030, "center", "center"),
                               (x * 1.15, 0.070, "center", "center"),
                               (x * 3.4, 0.055, "center", "center")],
                    fontsize=6.0, color=QUAL[i], weight="bold",
                    bbox=dict(boxstyle="square,pad=0.12", facecolor="white",
                              edgecolor="none", alpha=0.85))
    key_lines = [("coastal lowlands — 11 wells", "S100CSLLWD"),
                 ("MS alluvium — 12", "N100MSRVVL"),
                 ("N. Atlantic CP — 39", "S100NATLCP"),
                 ("glacial — 71", "N100GLCIAL"),
                 ("High Plains — 81", "N100HGHPLN"),
                 ("Basin & Range — never <10%", "N100BSNRGB")]
    for k, (txt, aqk) in enumerate(key_lines):
        i = list(big.index).index(aqk)
        axb.text(0.98, 0.985 - k * 0.082, txt, transform=axb.transAxes,
                 fontsize=5.4, color=QUAL[i], ha="right", va="top", zorder=7,
                 bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                           edgecolor="none", alpha=0.82))
    axb.set_xlabel("Wells retained in aquifer network, $n$")
    axb.set_ylabel("P(trend sign flips)")
    axb.text(11500, 0.028, "safe  <10%", fontsize=5.8, color=C_SAFE, ha="right")
    axb.text(11500, 0.20, "10–30%", fontsize=5.8, color=C_WARN, ha="right")
    axb.text(9.3, 0.412, ">30%", fontsize=5.8, color=C_RISK, ha="left", va="top")
    axb.set_title("Sign-reversal probability across six major aquifers",
                  fontsize=7, pad=3)
    axb.text(-0.155, 1.16, "a", transform=axb.transAxes, fontsize=9,
             fontweight="bold")

    # ---- panel b: CONUS random vs clustered, with the law and effective size
    axc = fig.add_subplot(gs[0, 1])
    conus = res[(res.stratum == "CONUS") & (~res.clustered)].sort_values("n_keep")
    clus = res[(res.stratum == "CONUS") & res.clustered].sort_values("n_keep")
    risk_zones(axc, 0.42)
    snr = law_snr_conus()
    nn = np.logspace(2.35, 4.45, 240)
    nn = nn[nn < 29717 * 0.999]
    axc.plot(nn, norm.cdf(-C_LAW * snr * np.sqrt(nn * 29717 / (29717 - nn))),
             color=C_NEUTRAL, lw=1.0, ls=(0, (2, 2)), zorder=2,
             label="sampling law (random)")
    p = conus.p_sign_flip.values
    nrep = conus.n_rep.values
    half = 1.96 * np.sqrt(np.clip(p * (1 - p), 1e-4, None) / nrep)
    axc.fill_between(conus.n_keep, np.clip(p - half, 0, 1), np.clip(p + half, 0, 1),
                     color=C_NEUTRAL, alpha=0.35, lw=0, zorder=2)
    axc.plot(conus.n_keep, p, "-o", color=C_DARK, ms=2.4, lw=1.2, zorder=4,
             label="random loss")
    axc.plot(clus.n_keep, clus.p_sign_flip, "-s", color=C_ACCENT, ms=2.6,
             lw=1.2, zorder=4, label="clustered loss\n(budget-cut pattern)")
    axc.set_ylim(0, 0.42)
    axc.set_xlim(200, 33000)
    axc.set_xscale("log")
    axc.set_xlabel("Wells retained in national (CONUS) network, $n$")
    axc.set_ylabel("P(trend sign flips)")
    axc.plot([297], [0.10], "o", ms=4.5, mfc=C_DARK, mec="white", mew=0.7, zorder=6)
    axc.plot([7726], [0.10], "s", ms=4.5, mfc=C_ACCENT, mec="white", mew=0.7, zorder=6)
    axc.annotate("297", (297, 0.10), textcoords="offset points",
                 xytext=(0, 11), fontsize=5.6, color=C_DARK, ha="center")
    axc.annotate("7,726", (7726, 0.10), textcoords="offset points",
                 xytext=(0, 11), fontsize=5.6, color=C_ACCENT, ha="center")
    axc.annotate("", xy=(300, 0.108), xytext=(7726, 0.128),
                 arrowprops=dict(arrowstyle="->", color=C_ACCENT, lw=1.0,
                                 connectionstyle="arc3,rad=-0.22"))
    axc.text(420, 0.295, "7,726 clustered wells carry the\nsame flip risk as ~180–280 random\nwells (surviving 0.5° blocks)",
             fontsize=5.8, color=C_DARK, ha="left", zorder=7,
             bbox=dict(boxstyle="square,pad=0.2", facecolor="white",
                       edgecolor="none", alpha=0.85))
    axc.legend(loc="upper right", fontsize=5.6, handlelength=1.5,
               frameon=True, framealpha=0.9, edgecolor="none")
    axc.set_title("National network (CONUS): random vs. clustered loss",
                  fontsize=7, pad=3)
    axc.text(-0.155, 1.16, "b", transform=axc.transAxes, fontsize=9,
             fontweight="bold")

    save_pub(fig, "fig2_flip_curves", TWO_COL)


# ---------------------------------------------------------------- fig 2
STATE_ABBR = {
    "Washington": "WA", "California": "CA", "Montana": "MT", "Minnesota": "MN",
    "Wisconsin": "WI", "Nebraska": "NE", "Texas": "TX", "Arkansas": "AR",
    "Georgia": "GA", "Florida": "FL", "New York": "NY", "Idaho": "ID",
}


def state_anchors(ax):
    """Two-letter state labels at rough centroids, behind the well points."""
    gj = json.load(open(os.path.join(DATA, "us-states.json")))
    for feat in gj["features"]:
        name = feat["properties"]["name"]
        if name not in STATE_ABBR:
            continue
        geom = feat["geometry"]
        rings = ([geom["coordinates"]] if geom["type"] == "Polygon"
                 else geom["coordinates"])
        ring = max((r for poly in rings for r in poly), key=len)
        ring = np.asarray(ring)
        cx, cy = ring[:, 0].mean(), ring[:, 1].mean()
        if not (-124 < cx < -67 and 25 < cy < 49.5):
            continue  # centroid would clip at the map frame
        ax.text(cx, cy, STATE_ABBR[name],
                fontsize=5.2, color="#7A8089", ha="center", va="center",
                zorder=1.5, fontweight="bold")


def fig2(res, wells, law):
    fig = plt.figure(figsize=(TWO_COL, 2.36))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.52,
                          left=0.055, right=0.985, top=0.88, bottom=0.15)

    ax = fig.add_subplot(gs[0, 0])
    gj = json.load(open(os.path.join(DATA, "us-states.json")))
    for feat in gj["features"]:
        geom = feat["geometry"]
        polys = (geom["coordinates"] if geom["type"] == "Polygon"
                 else [p for poly in geom["coordinates"] for p in poly])
        for poly in polys:
            ring = np.asarray(poly)
            if ring[:, 0].max() < -128 or ring[:, 0].min() > -64:
                continue
            if ring[:, 1].max() < 23 or ring[:, 1].min() > 51:
                continue
            ax.plot(ring[:, 0], ring[:, 1], color=GRAY_LIGHT, lw=1.0, zorder=1)
    state_anchors(ax)
    rng = np.random.default_rng(7)
    jit = rng.normal(0, 0.08, len(wells))
    dec = wells.slope_ft_yr > 0
    ris = wells.slope_ft_yr < 0
    ax.scatter(wells.x[dec] + jit[dec], wells.y[dec] + jit[dec], s=0.35,
               c=C_SIGNAL2, alpha=0.55, lw=0, rasterized=True, zorder=2)
    ax.scatter(wells.x[ris] + jit[ris], wells.y[ris] + jit[ris], s=0.35,
               c=C_SIGNAL, alpha=0.55, lw=0, rasterized=True, zorder=2)
    ax.scatter(wells.x[~dec & ~ris] + jit[~dec & ~ris], wells.y[~dec & ~ris]
               + jit[~dec & ~ris], s=0.35, c=C_NEUTRAL, alpha=0.55, lw=0,
               rasterized=True, zorder=2)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_aspect(1.0 / np.cos(np.radians(39)), adjustable="datalim")
    ax.set_xlim(-126, -65)
    ax.set_ylim(23, 51)
    ax.autoscale_view()
    ax.set_title("Long-term monitoring wells (n = 29,717; 1950–2026)",
                 fontsize=7, pad=3)
    handles = [Line2D([], [], marker="o", ls="", ms=2.8, color=C_SIGNAL2,
                      label="falling: %d wells" % dec.sum()),
               Line2D([], [], marker="o", ls="", ms=2.8, color=C_SIGNAL,
                      label="rising: %d wells" % ris.sum()),
               Line2D([], [], marker="o", ls="", ms=2.8, color=C_NEUTRAL,
                      label="stable: %d" % (~dec & ~ris).sum())]
    ax.legend(handles=handles, loc="lower left", fontsize=5.8, handletextpad=0.3)
    ax.text(-0.10, 1.09, "a", transform=ax.transAxes, fontsize=9,
            fontweight="bold")

    # ---- panel b: minimum reliable size, simulation vs law (dumbbell rows)
    ax = fig.add_subplot(gs[0, 1])
    base = (res[(~res.clustered) & (~res.stratum.isin(["CONUS", "UNKNOWN"]))]
            .groupby("stratum").n_all.first().sort_values(ascending=False))
    th_codes = list(base.head(9).index) + [c for c in ["N100BSNRGB", "S100SURFCL",
                                                       "S400FLORDN"]
                                           if c not in base.head(9).index]
    d = law[law.stratum.isin(th_codes)].copy()
    d["ret"] = d.nmin_sim / d.N
    d = d.sort_values("ret", na_position="last").reset_index(drop=True)
    y = np.arange(len(d))
    xmax = 13000
    for i, r in d.iterrows():
        ret = r.ret if np.isfinite(r.ret) else 1.0
        col = regime_color(ret)
        # full-network reference
        ax.plot([min(r.N, xmax), xmax], [i, i], color=GRAY_LIGHT, lw=1.6, zorder=1)
        ax.plot(r.N, i, "o", ms=3.8, mfc="white", mec=C_NEUTRAL, mew=0.9, zorder=3)
        if np.isfinite(r.ret):
            ax.plot([r.nmin_sim, r.N], [i, i], color=col, lw=1.2, zorder=2)
            ax.plot(r.nmin_sim, i, "o", ms=3.8, color=col, zorder=4)
            if r.nmin_sim < 15:   # near left edge: label to the right, not above
                ax.annotate("%d" % r.nmin_sim, (r.nmin_sim, i),
                            textcoords="offset points", xytext=(8, -2.5),
                            ha="left", va="center", fontsize=5.2, color=col)
            else:
                ax.annotate("%d" % r.nmin_sim, (r.nmin_sim, i),
                            textcoords="offset points", xytext=(-1, 4),
                            ha="center", fontsize=5.2, color=col)
        else:
            ax.plot(0.9 * r.N, i, ">", ms=3.8, color=C_RISK, zorder=4)
            ax.annotate("never ≤10%% (law: %d%%)" % round(100 * r.ret_pred),
                        (0.9 * r.N, i), textcoords="offset points",
                        xytext=(2, 4), ha="right", fontsize=5.2, color=C_RISK)
        ax.plot(r.nmin_pred, i, "D", ms=3.0, mfc=C_DARK, mec=C_DARK, zorder=5,
                alpha=0.9)
    ax.set_xscale("log")
    ax.set_xlim(6, xmax)
    ax.set_ylim(-0.8, len(d) - 0.2)
    ax.set_yticks(y)
    ax.set_yticklabels([aq_name(a).replace("\n", " ") for a in d.stratum],
                       fontsize=5.6)
    ax.set_xlabel("Wells needed for ≤10% flip risk (log scale)")
    ax.set_title("Minimum network size for ≤10% sign-reversal risk", fontsize=7,
                 pad=3, loc="left")
    ax.text(-0.34, 1.09, "b", transform=ax.transAxes, fontsize=9,
            fontweight="bold")

    save_pub(fig, "fig1_map_thresholds", TWO_COL)


# ---------------------------------------------------------------- fig 3
DESIGN_LABELS = [
    # stratum, label, (dx, dy) in points, ha
    ("N100HGHPLN", "High Plains\n8,160 wells · needs ~1%", (7, -3), "left"),
    ("S100CSLLWD", "coastal lowlands\n1,145 · ~1%", (-9, -14), "right"),
    ("CONUS", "CONUS\n29,717 · ~1%", (-8, -12), "right"),
    ("N100BSNRGB", "Basin & Range\n1,557 · needs 94%", (8, 0), "left"),
    ("S400FLORDN", "Floridan\n819 · 88%", (8, -2), "left"),
    ("S100SURFCL", "surficial\n195 · 91%", (8, 2), "left"),
]


def fig3(law, coll):
    fig, axes = plt.subplots(1, 2, figsize=(TWO_COL, 2.40))
    fig.subplots_adjust(wspace=0.36, left=0.075, right=0.925, top=0.86,
                        bottom=0.215)
    ret_map = dict(zip(law.stratum, law.ret_pred))

    # ---- panel a: the collapse
    ax = axes[0]
    x = C_LAW * coll.snr * coll.geom
    cols = [regime_color(ret_map.get(s, 0.5)) for s in coll.stratum]
    size = 6 + 26 * np.sqrt(coll.N / coll.N.max())
    ax.scatter(x, coll.flip, s=size, c=cols, alpha=0.55, lw=0, zorder=3)
    xs = np.logspace(-1.05, 1.62, 300)
    ax.plot(xs, norm.cdf(-xs), color=C_DARK, lw=1.6, zorder=4)
    ax.axhline(0.10, color=C_NEUTRAL, lw=1.0, ls=(0, (3, 3)), zorder=2)
    ax.axvline(Z, color=C_RISK, lw=1.0, ls=(0, (3, 3)), zorder=2)
    ax.set_xscale("log")
    ax.set_xlim(0.09, 42)
    ax.set_ylim(0, 0.8)
    ax.set_xlabel("Normalized size  $c\\,\\mathrm{SNR}\\,(nN/(N{-}n))^{1/2}$")
    ax.set_ylabel("P(trend sign flips)")
    ax.text(40, 0.30, "theory: $\\Phi(-x)$\nzero fitted shape parameter",
            fontsize=6.0, color=C_DARK, ha="right", va="top", zorder=6,
            bbox=dict(boxstyle="square,pad=0.2", facecolor="white",
                      edgecolor="none", alpha=0.85))
    ax.text(Z * 1.06, 0.70, "10% risk", fontsize=5.8, color=C_RISK,
            rotation=90, va="top", ha="left")
    conus = coll[coll.stratum == "CONUS"]
    xc = C_LAW * conus.snr.iloc[0] * conus.geom.values
    ax.plot(xc[np.argmin(np.abs(conus.flip.values - 0.075))],
            0.075, "*", ms=7, color=C_SIGNAL2, mec="white", mew=0.5, zorder=6)
    ax.annotate("CONUS", (xc[np.argmin(np.abs(conus.flip.values - 0.075))], 0.075),
                textcoords="offset points", xytext=(6, -7), fontsize=5.6,
                color=C_SIGNAL2, zorder=7,
                bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.85))
    ax.set_title("Collapse of 216 thinning trials onto the sampling law", fontsize=7, pad=3)
    # regime key in the empty lower-left (points sit high-left / low-right)
    reg = [Line2D([], [], marker="o", ls="", ms=3.6, color=C_SAFE,
                  label="net needs ≤5%"),
           Line2D([], [], marker="o", ls="", ms=3.6, color=C_WARN,
                  label="5–50%"),
           Line2D([], [], marker="o", ls="", ms=3.6, color=C_RISK,
                  label=">50%")]
    leg = ax.legend(handles=reg, loc="upper right", fontsize=5.6,
                    handletextpad=0.3, labelspacing=0.35, frameon=False,
                    borderaxespad=0.4)
    leg.set_title("marker area scales with $N$", prop={"size": 5.4})
    leg.get_title().set_color(C_DARK)
    ax.text(-0.155, 1.12, "a", transform=ax.transAxes, fontsize=9,
            fontweight="bold")

    # ---- panel b: the design chart
    ax = axes[1]
    snr_grid = np.logspace(-2.45, -0.42, 250)
    # label each isopleth at a staggered interior point (target N per line)
    label_N = {0.01: 60000, 0.05: 45000, 0.10: 25000, 0.30: 8000, 0.90: 4000}
    for q in [0.01, 0.05, 0.10, 0.30, 0.90]:
        acol = C_RISK if q == 0.90 else ("#7E848D" if q == 0.10 else "#9AA0A8")
        Nq = Z ** 2 * (1 - q) / (q * C_LAW ** 2) / snr_grid ** 2
        m = (Nq > 70) & (Nq < 90000)
        ax.plot(snr_grid[m], Nq[m], lw=1.0 if q != 0.90 else 1.2,
                color=acol, ls="-" if q != 0.90 else (0, (4, 2)), zorder=2)
        Nlab = min(label_N[q], 60000)
        snr_lab = np.sqrt(Z ** 2 * (1 - q) / (q * C_LAW ** 2) / Nlab)
        ax.annotate("%d%%" % (q * 100), (snr_lab, Nlab),
                    textcoords="offset points", xytext=(3, 0),
                    fontsize=5.6, color=acol, ha="left", va="center")
    # not-separable region below the 90% line
    N90 = Z ** 2 * 0.1 / (0.9 * C_LAW ** 2) / snr_grid ** 2
    ax.fill_between(snr_grid, 60, N90, color=C_RISK, alpha=0.05, lw=0, zorder=1)
    ax.text(0.0045, 68, "needs ≥90% of wells:\ntrend not separable from noise",
            fontsize=5.6, color=C_RISK, va="bottom", zorder=7,
            bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                      edgecolor="none", alpha=0.85))
    ax.text(0.28, 30000, "high SNR:\n≤5% of wells suffice",
            fontsize=5.6, color=C_SAFE, va="top", ha="center", zorder=7,
            bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                      edgecolor="none", alpha=0.85))

    d = law[law.stratum != "UNKNOWN"]
    pos = {s: (r.snr, r.N) for s, r in d.set_index("stratum").iterrows()}
    for _, r in d.iterrows():
        ret = r.ret_pred
        col = regime_color(ret)
        star = r.stratum == "CONUS"
        ax.plot(r.snr, r.N, "*" if star else "o",
                ms=8 if star else 4.2, color=col, mec="white", mew=0.5,
                zorder=5)
    for s, lab, (dx, dy), ha in DESIGN_LABELS:
        if s not in pos:
            continue
        sx, sy = pos[s]
        ax.annotate(lab, (sx, sy), textcoords="offset points",
                    xytext=(dx, dy), fontsize=5.2, color=C_DARK, ha=ha,
                    va="center", zorder=7,
                    bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                              edgecolor="none", alpha=0.85))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.0042, 0.62)
    ax.set_ylim(60, 70000)
    ax.set_xlabel("Trend signal-to-noise  SNR =\n|network trend| / SD(well trends)")
    ax.set_ylabel("Network size $N$ (wells)")
    ax.set_title("Network design chart: required retention fraction",
                 fontsize=7, pad=3)
    ax.text(-0.155, 1.12, "b", transform=ax.transAxes, fontsize=9,
            fontweight="bold")

    save_pub(fig, "fig3_signal_noise", TWO_COL)


# ---------------------------------------------------------------- fig 4
def fig4(res, strat):
    fig, axes = plt.subplots(1, 2, figsize=(TWO_COL, 2.24))
    fig.subplots_adjust(wspace=0.30, left=0.085, right=0.99, top=0.84,
                        bottom=0.20)

    for k, (st, title) in enumerate(
            [("N100HGHPLN",
              "High Plains (strong signal): no order flips the sign;\norder moves the bias"),
             ("S100CNRLVL",
              "Central Valley (weak signal): spread holds the sign,\n"
              "not the magnitude")]):
        ax = axes[k]
        h = strat[strat.stratum == st]
        nr = h[(h.strategy == "nested_random") & (h.abs_bias > 0)].sort_values("n_keep")
        lf = h[(h.strategy == "longest_first") & (h.abs_bias > 0)].sort_values("n_keep")
        sp = h[(h.strategy == "spatial_spread") & (h.abs_bias > 0)].sort_values("n_keep")
        if len(nr):
            nrep = nr.stratum_n.iloc[0]
            half = 1.96 * np.sqrt(np.clip(nr.abs_bias, 1e-6, None) / max(nrep, 1))
            ax.fill_between(nr.n_keep, np.clip(nr.abs_bias - half, 1e-4, None) * 0.3048,
                            (nr.abs_bias + half) * 0.3048, color=C_NEUTRAL,
                            alpha=0.45, lw=0, zorder=2)
            ax.plot(nr.n_keep, nr.abs_bias * 0.3048, color=C_DARK, lw=1.0, zorder=3)
        ax.plot(lf.n_keep, lf.abs_bias * 0.3048, "-", color=C_SIGNAL2, lw=1.3, zorder=4)
        ax.plot(sp.n_keep, sp.abs_bias * 0.3048, "--", color=C_ACCENT, lw=1.3, zorder=4)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(3e-4, 0.4)
        ax.set_title(title, fontsize=6.5, pad=3)
        ax.set_xlabel("Wells retained, $n$")
        if k == 0:
            ax.set_ylabel("Slope bias |$\\hat{\\beta}-\\beta|$ (m yr$^{-1}$)")
            place_clear(ax, "random orders\n(mean ±95%)",
                        [(470, 0.018, "left", "center"), (450, 0.09, "left", "center"),
                         (900, 0.30, "left", "center"), (3000, 0.16, "left", "center")],
                        fontsize=5.4, color=C_DARK)
            place_clear(ax, "longest records first",
                        [(600, 0.30, "left", "center"), (700, 0.22, "left", "center"),
                         (450, 0.38, "left", "center")], fontsize=5.6, color=C_SIGNAL2)
            place_clear(ax, "spatial spread first",
                        [(550, 0.0016, "left", "center"), (700, 0.0009, "left", "center"),
                         (450, 0.004, "left", "center")], fontsize=5.6, color=C_ACCENT)
            ax.annotate("5.4× the random-order bias\nat 816 wells",
                        xy=(816, 0.0166), xytext=(1800, 0.06),
                        fontsize=5.6, color=C_SIGNAL2,
                        arrowprops=dict(arrowstyle="-", color=C_SIGNAL2, lw=1.0))
        else:
            ax.text(190, 0.0012, "random and longest-first\nflip below ~150 wells",
                    fontsize=5.6, color=C_DARK)
            # identify the series in panel b too (colours shared with panel a)
            ax.plot([], [], color=C_DARK, lw=1.0, label="random orders (±95%)")
            ax.plot([], [], color=C_SIGNAL2, lw=1.3, label="longest records first")
            ax.plot([], [], color=C_ACCENT, lw=1.3, ls=(0, (4, 2)),
                    label="spatial spread first")
            ax.legend(loc="upper right", fontsize=5.4, handlelength=1.6,
                      frameon=True, framealpha=0.92, edgecolor="none",
                      borderaxespad=0.4)
        ax.text(-0.17, 1.15, "ab"[k], transform=ax.transAxes, fontsize=9,
                fontweight="bold")

    save_pub(fig, "figS_strategies", TWO_COL)


def law_snr_conus():
    coll = pd.read_csv(os.path.join(DATA, "law_collapse.csv"))
    return float(coll[coll.stratum == "CONUS"].snr.iloc[0])


if __name__ == "__main__":
    res, strat, wells, law, coll = load()
    fig1(res)
    fig2(res, wells, law)
    fig3(law, coll)
    fig4(res, strat)
    print("all figures done ->", FIGS)
