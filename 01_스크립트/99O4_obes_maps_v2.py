#!/usr/bin/env python3
# Figure 1 district maps at print size (170 mm wide, 8 pt text, scale bar, 600 dpi). Inputs: OBES_SHP (or ST_SHP), OBES_FO, OBES_RES, OBES_OUT.
import os, sys, unicodedata
import numpy as np, pandas as pd, geopandas as gpd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle, FancyArrow
from PIL import Image
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, ".."))
SHP = os.environ.get("OBES_SHP", os.environ.get("ST_SHP", "final.shp"))
FOF = os.environ.get("OBES_FO", os.path.join(ROOT, "02_패널", "OBES_fh_outcomes.csv"))
RES = os.environ.get("OBES_RES", os.path.join(ROOT, "04_분석", "OBES", "results_std3"))
OUTD = os.environ.get("OBES_OUT", os.path.join(ROOT, "04_분석", "OBES", "MS", "out", "Figures")); os.makedirs(OUTD, exist_ok=True)
MM = 1 / 25.4; FW = 170 * MM
FS = 8; PANEL_PT = 12
plt.rcParams.update({"font.family": "Arial", "font.size": FS, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "xtick.major.size": 2.5,
                     "hatch.linewidth": 0.5})
G = gpd.read_file(SHP); G["region"] = G.region.map(NFC)
FO = pd.read_csv(FOF, encoding="utf-8-sig"); FO["fbd_region"] = FO.fbd_region.map(NFC)
sr = FO.groupby("fbd_region").SR.mean().rename("v").reset_index(); ms = FO[FO.year.isin([2018, 2019])].groupby("fbd_region").meas_SR.mean().rename("v").reset_index()
def ld(f):
    r = pd.read_csv(os.path.join(RES, f)); r["fbd_region"] = r.fbd_region.map(NFC); return r.rename(columns={"OR": "v"})
P = [("A", sr, "seq", "Self-reported obesity, age-sex-standardized %,\n2011–2024 mean"), ("B", ld("S3_spatial_OR.csv"), "div", "Spatial odds ratio,\nself-reported obesity"),
     ("C", ms, "seq", "Measured obesity, age-sex-standardized %,\n2018–2019 mean"), ("D", ld("S14_measured_spatialOR.csv"), "div", "Spatial odds ratio,\nmeasured obesity")]

x0, x1, y0, y1 = G.total_bounds; padx = (x1 - x0) * 0.02; pady = (y1 - y0) * 0.02
XL, YL = (x0 - padx, x1 + padx), (y0 - pady, y1 + pady); ASP = (YL[1] - YL[0]) / (XL[1] - XL[0])
LM, GAPX = 0.10, 0.25; MW = (FW - 2 * LM - GAPX) / 2; MH = MW * ASP
TOP, LET, CBG, CBH, CBLAB, ROWG, LEGH = 0.06, 0.20, 0.06, 0.09, 0.50, 0.12, 0.42
PH = LET + MH + CBG + CBH + CBLAB
FH = TOP + 2 * PH + ROWG + LEGH
fig = plt.figure(figsize=(FW, FH), facecolor="white")
def A(l, b, w, h): return fig.add_axes([l / FW, b / FH, w / FW, h / FH])
for i, (lab, r, kind, title) in enumerate(P):
    col, row = i % 2, i // 2
    left = LM + col * (MW + GAPX); ptop = FH - TOP - row * (PH + ROWG)
    ax = A(left, ptop - LET - MH, MW, MH)
    cw = MW * 0.70; cax = A(left + (MW - cw) / 2, ptop - LET - MH - CBG - CBH, cw, CBH)
    cols = ["fbd_region", "v"] + (["class"] if "class" in r else [])
    g = G.merge(r[cols], left_on="region", right_on="fbd_region", how="left", validate="one_to_one")
    miss = g.v.isna()
    if kind == "div":
        x = np.abs(np.log(g.v.dropna())).max(); norm = TwoSlopeNorm(vmin=np.exp(-x), vcenter=1.0, vmax=np.exp(x)); cmap = "RdBu_r"
    else:
        norm = plt.Normalize(g.v.quantile(0.02), g.v.quantile(0.98)); cmap = "YlOrRd"
    g[~miss].plot(column="v", cmap=cmap, norm=norm, ax=ax, edgecolor="#9a9a9a", linewidth=0.08)
    if miss.any(): g[miss].plot(ax=ax, color="#d9d9d9", edgecolor="#6f6f6f", linewidth=0.08, hatch="//////")
    if "class" in g:
        g[g["class"] == "hotspot"].boundary.plot(ax=ax, color="black", linewidth=0.55)
        g[g["class"] == "coldspot"].boundary.plot(ax=ax, color="black", linewidth=0.55, linestyle=(0, (2, 1.2)))
    ax.set_xlim(XL); ax.set_ylim(YL); ax.set_aspect("equal"); ax.set_axis_off()
    ax.add_patch(FancyArrow(0.94, 0.855, 0, 0.075, width=0.012, head_width=0.040, head_length=0.035, length_includes_head=True,
                            transform=ax.transAxes, color="black", lw=0))
    ax.text(0.94, 0.94, "N", transform=ax.transAxes, ha="center", va="bottom", fontsize=FS, fontweight="bold")
    sx, sy, bh = XL[1] - 0.09 * (XL[1] - XL[0]) - 200_000, YL[0] + 0.06 * (YL[1] - YL[0]), 0.012 * (YL[1] - YL[0])
    for k, fc in enumerate(["black", "white"]):
        ax.add_patch(Rectangle((sx + k * 100_000, sy), 100_000, bh, facecolor=fc, edgecolor="black", lw=0.5))
    for d in (0, 100, 200):
        ax.text(sx + d * 1000, sy + bh * 1.6, str(d), ha="center", va="bottom", fontsize=FS)
    ax.text(sx + 200_000 + 0.03 * (XL[1] - XL[0]), sy + bh / 2, "km", ha="left", va="center", fontsize=FS)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label(title, fontsize=FS, labelpad=2); cb.ax.tick_params(labelsize=FS, pad=1.5); cb.outline.set_linewidth(0.5)
    if kind == "div":
        t = [v for v in np.arange(0.80, 1.2501, 0.05) if norm.vmin <= v <= norm.vmax]; cb.set_ticks(t); cb.set_ticklabels([f"{v:.2f}" for v in t])
    fig.text((left - 0.05) / FW, (ptop - 0.02) / FH, lab, fontsize=PANEL_PT, fontweight="bold", ha="left", va="top")
fig.legend(handles=[Line2D([0], [0], color="black", lw=0.9, label="Hotspot, posterior probability > 0.95"),
                    Line2D([0], [0], color="black", lw=0.9, ls=(0, (2, 1.2)), label="Coldspot, posterior probability < 0.05"),
                    Patch(facecolor="#d9d9d9", edgecolor="#6f6f6f", hatch="//////", lw=0.5, label="No measured data (Sejong)")],
           loc="lower center", ncol=3, fontsize=FS, frameon=False, bbox_to_anchor=(0.5, 0.04 / FH), handlelength=2.2, columnspacing=1.2, handletextpad=0.5)
tmp = os.path.join(OUTD, "_Figure1_rgba.tiff")
fig.savefig(tmp, dpi=600, facecolor="white")
im = Image.open(tmp).convert("RGB"); im.save(os.path.join(OUTD, "Figure1.tiff"), compression="tiff_lzw", dpi=(600, 600)); os.remove(tmp)
fig.savefig(os.path.join(OUTD, "Figure1.png"), dpi=300, facecolor="white")
print("saved", OUTD, im.size, f"{FW / MM:.0f} x {FH / MM:.0f} mm")
