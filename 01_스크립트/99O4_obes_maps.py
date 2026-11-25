#!/usr/bin/env python3
# District maps.
import os, unicodedata
import numpy as np, pandas as pd, geopandas as gpd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
NFC = lambda s: unicodedata.normalize("NFC", str(s))
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")); AN = os.path.join(ROOT, "04_분석"); RES = os.path.join(AN, "OBES", "results_std3")
SHP = os.environ.get("ST_SHP", "final.shp")
OUTD = os.path.join(AN, "OBES", "figs"); os.makedirs(OUTD, exist_ok=True)
PANEL_LABEL_PT = 60; plt.rcParams.update({"font.family": "Arial", "font.size": 16})
G = gpd.read_file(SHP); G["region"] = G.region.map(NFC)
FO = pd.read_csv(os.path.join(ROOT, "02_패널", "OBES_fh_outcomes.csv"), encoding="utf-8-sig"); FO["fbd_region"] = FO.fbd_region.map(NFC)
sr = FO.groupby("fbd_region").SR.mean().rename("v").reset_index(); ms = FO[FO.year.isin([2018, 2019])].groupby("fbd_region").meas_SR.mean().rename("v").reset_index()
def ld(f):
    r = pd.read_csv(os.path.join(RES, f)); r["fbd_region"] = r.fbd_region.map(NFC); return r.rename(columns={"OR": "v"})
P = [("A", sr, "seq", "Self-reported obesity,\nage-sex-standardized %,\n2011–2024 mean"), ("B", ld("S3_spatial_OR.csv"), "div", "Spatial odds ratio,\nself-reported obesity"),
     ("C", ms, "seq", "Measured obesity,\nage-sex-standardized %,\n2018–2019 mean"), ("D", ld("S14_measured_spatialOR.csv"), "div", "Spatial odds ratio,\nmeasured obesity")]
fig, axes = plt.subplots(2, 2, figsize=(16, 15))
for ax, (lab, r, kind, title) in zip(axes.flat, P):
    cols = ["fbd_region", "v"] + (["class"] if "class" in r else [])
    g = G.merge(r[cols], left_on="region", right_on="fbd_region", how="left", validate="one_to_one")
    miss = g.v.isna()
    if kind == "div":
        x = np.abs(np.log(g.v.dropna())).max(); norm = TwoSlopeNorm(vmin=np.exp(-x), vcenter=1.0, vmax=np.exp(x)); cmap = "RdBu_r"
    else:
        norm = plt.Normalize(g.v.quantile(0.02), g.v.quantile(0.98)); cmap = "YlOrRd"
    g[~miss].plot(column="v", cmap=cmap, norm=norm, ax=ax, edgecolor="#9a9a9a", linewidth=0.15)
    if miss.any(): g[miss].plot(ax=ax, color="#d9d9d9", edgecolor="#9a9a9a", linewidth=0.15, hatch="///")
    if "class" in g:
        g[g["class"] == "hotspot"].boundary.plot(ax=ax, color="black", linewidth=0.9)
        g[g["class"] == "coldspot"].boundary.plot(ax=ax, color="black", linewidth=0.9, linestyle=(0, (2, 1.5)))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); cb = fig.colorbar(sm, ax=ax, shrink=0.45, pad=0.01, location="right"); cb.set_label(title, fontsize=15); cb.ax.tick_params(labelsize=13)
    ax.set_axis_off()
    ax.annotate("", xy=(0.86, 0.36), xytext=(0.86, 0.10), xycoords="axes fraction", arrowprops=dict(facecolor="black", width=27, headwidth=81, headlength=70))
    ax.text(0.86, 0.40, "N", transform=ax.transAxes, ha="center", va="bottom", fontsize=48, fontweight="bold")
plt.subplots_adjust(left=0.05, right=0.97, top=0.95, bottom=0.05, wspace=0.05, hspace=0.12); fig.canvas.draw()
for (lab, *_), ax in zip(P, axes.flat):
    b = ax.get_position(); fig.text(b.x0 - 0.03, b.y1 + 0.01, lab, fontsize=PANEL_LABEL_PT, fontweight="bold", ha="left", va="bottom")
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
fig.legend(handles=[Line2D([0], [0], color="black", lw=1.5, label="Hotspot, posterior probability > 0.95"), Line2D([0], [0], color="black", lw=1.5, ls=(0, (2, 1.5)), label="Coldspot, posterior probability < 0.05"),
                    Patch(facecolor="#d9d9d9", edgecolor="#9a9a9a", hatch="///", label="No measured data (Sejong)")], loc="lower center", ncol=3, fontsize=14, frameon=False, bbox_to_anchor=(0.5, 0.0))
for ext in ["png", "tiff"]: fig.savefig(os.path.join(OUTD, f"Fig_maps_OBES.{ext}"), dpi=300 if ext == "png" else 600, pil_kwargs={"compression": "tiff_lzw"} if ext == "tiff" else None)
print("[99O4] 저장", OUTD)
