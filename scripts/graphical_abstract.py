#!/usr/bin/env python3
"""Graphical abstract (law version): 5.5x5.5 in, 300 dpi,
TIFF + PDF. Two panels: the law collapse (a) and the design chart (b),
with c and z defined on-artwork. Fresh canvas: no subtitle ghosting."""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

mpl.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "pdf.fonttype": 42, "svg.fonttype": "none",
                     "mathtext.fontset": "custom",
                     "mathtext.rm": "Arial",
                     "mathtext.it": "Arial:italic",
                     "mathtext.bf": "Arial:bold",
                     "axes.spines.right": False, "axes.spines.top": False})
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))
FIGS = os.path.abspath(os.path.join(HERE, "..", "figures"))
C_SIG, C_DARK, C_NEU, C_RISK = "#2166AC", "#252525", "#8A8F98", "#B2182B"
C = 1.27
Z = float(norm.ppf(0.9))

coll = pd.read_csv(os.path.join(DATA, "law_collapse.csv"))
law = pd.read_csv(os.path.join(DATA, "law_nmin.csv"))

fig, axes = plt.subplots(1, 2, figsize=(5.5, 5.5),
                         gridspec_kw=dict(wspace=0.44, width_ratios=[1.1, 1.0]))
ax = axes[0]
x = C * coll.snr * coll.geom
ax.scatter(x, coll.flip, s=5 + 20 * np.sqrt(coll.N / coll.N.max()),
           color=C_SIG, alpha=0.4, lw=0)
xs = np.logspace(-1.0, 1.6, 300)
ax.plot(xs, norm.cdf(-xs), color=C_DARK, lw=1.8)
ax.axhline(0.10, color=C_NEU, lw=1.0, ls=(0, (4, 3)))
ax.set_xscale("log")
ax.set_xlim(0.11, 40)
ax.set_ylim(0, 1)
ax.set_xlabel("normalized size c\u00b7SNR\u00b7\u221a(nN/(N\u2212n))", fontsize=9)
ax.set_ylabel("P(trend sign flips)", fontsize=9)
ax.tick_params(labelsize=8)
# axes-fraction placement: fully inside the plot, clear of tick labels
ax.text(0.40, 0.80, "216 thinning trials\nfollow Φ(\u2212x)", fontsize=9,
        color=C_DARK, ha="center", va="center", transform=ax.transAxes,
        zorder=7,
        bbox=dict(boxstyle="square,pad=0.25", facecolor="white",
                  edgecolor="none", alpha=0.9))
ax.set_title("Finite-population sampling law", fontsize=10, fontweight="bold", color=C_DARK)

ax = axes[1]
snr = np.logspace(-2.2, -0.5, 220)
for q, col in [(0.05, "#8FB7D9"), (0.30, "#5E8FBE"), (0.90, C_RISK)]:
    Nq = Z ** 2 * (1 - q) / (q * C ** 2) / snr ** 2
    m = (Nq > 80) & (Nq < 6e4)
    ax.plot(snr[m], Nq[m], lw=1.6 if q != 0.9 else 2.0,
            color=col, ls="-" if q != 0.9 else (0, (5, 3)))
    # label each isopleth at its own height, off the bottom edge and clear
    # of the neighbours (all three curves end at the N floor and would
    # otherwise pile up on the axis)
    anchor_N = {0.05: 2.0e4, 0.30: 2.0e3, 0.90: 1.5e3}[q]
    xa = np.sqrt(Z ** 2 * (1 - q) / (q * C ** 2 * anchor_N))
    # 30% label flips to the empty left side of its curve; the other two
    # sit on their curve's clear side already
    dx, halign = ((-6, "right") if q == 0.30 else (4, "left"))
    ax.annotate("keep %d%%" % (q * 100), (xa, anchor_N),
                textcoords="offset points", xytext=(dx, 16), ha=halign,
                fontsize=8, color=col, zorder=7,
                bbox=dict(boxstyle="square,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.85))
d = law[law.stratum != "UNKNOWN"]
for _, r in d.iterrows():
    ret = r.ret_pred
    colr = "#4A9B4E" if ret <= 0.05 else ("#E0A030" if ret <= 0.5 else C_RISK)
    ax.plot(r.snr, r.N, "o", ms=3.5, color=colr, mec="white", mew=0.4)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0.006, 0.9)
ax.set_ylim(80, 6e4)
ax.set_xlabel("network SNR", fontsize=9)
ax.set_ylabel("network size N (wells)", fontsize=9)
ax.tick_params(labelsize=8)
ax.set_title("Minimum network design chart", fontsize=10, fontweight="bold", color=C_DARK)
fig.subplots_adjust(left=0.14, right=0.97, top=0.86, bottom=0.15)
fig.suptitle("How sparse is too sparse?", fontsize=14, fontweight="bold",
             color=C_DARK, y=0.985)
fig.text(0.5, 0.04,
         "n$_{\\mathrm{min}}$ = Nz$^2$/(z$^2$+(c\u00b7SNR)$^2$N);   c = 1.27;   "
         "z = normal quantile at the risk tolerance",
         fontsize=8, color=C_DARK, ha="center")
base = os.path.join(FIGS, "graphical_abstract")
fig.savefig(base + ".pdf", bbox_inches="tight")
fig.savefig(base + ".tiff", dpi=300, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(base + ".png", dpi=300, bbox_inches="tight")

# Journal spec is exactly 5.5 in square at 300 dpi (1650 x 1650 px); the
# tight bounding box trims whitespace, so centre the artwork on a white
# canvas of that size instead of stretching it.
from PIL import Image  # noqa: E402
SIDE = int(round(5.5 * 300))
for ext, kw in ((".tiff", {"compression": "tiff_lzw"}), (".png", {})):
    im = Image.open(base + ext).convert("RGB")
    assert im.width <= SIDE and im.height <= SIDE, im.size
    canvas = Image.new("RGB", (SIDE, SIDE), "white")
    canvas.paste(im, ((SIDE - im.width) // 2, (SIDE - im.height) // 2))
    canvas.save(base + ext, dpi=(300, 300), **kw)
print("GA regenerated ->", base, "%dx%d px @300 dpi" % (SIDE, SIDE))
