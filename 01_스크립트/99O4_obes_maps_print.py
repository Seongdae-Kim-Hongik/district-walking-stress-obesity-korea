#!/usr/bin/env python3
# Figure 1 district maps at print size (self-reported and measured obesity, spatial odds ratios), joined by district code. Environment: OBES_SHP or ST_SHP, OBES_FO, OBES_RES, OBES_PANEL, OBES_OUT.
import os, sys, unicodedata
import numpy as np, pandas as pd, geopandas as gpd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle, Polygon
from PIL import Image
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, ".."))
SHP = os.environ.get("OBES_SHP", os.environ.get("ST_SHP", "data/boundary/final.shp"))
FOF = os.environ.get("OBES_FO", os.path.join(ROOT, "02_패널", "OBES_fh_outcomes.csv"))
RES = os.environ.get("OBES_RES", os.path.join(ROOT, "04_분석", "OBES", "results_std3"))
PAN = os.environ.get("OBES_PANEL", os.path.join(ROOT, "04_분석", "OBES", "panel_OBES_std3.csv"))
OUTD = os.environ.get("OBES_OUT", os.path.join(ROOT, "04_분석", "OBES", "MS", "out", "Figures")); os.makedirs(OUTD, exist_ok=True)
MM = 1 / 25.4; FW = 170 * MM
FS = 8; PANEL_PT = 12
plt.rcParams.update({"font.family": "Arial", "font.size": FS, "axes.linewidth": 0.6, "hatch.linewidth": 0.5})
G = gpd.read_file(SHP)
if G.crs is None: G = G.set_crs(epsg=5179)
assert G.crs.to_epsg() == 5179; x0_, y0_, x1_, y1_ = G.total_bounds; assert 7.0e5 < x0_ < x1_ < 1.45e6 and 1.40e6 < y0_ < y1_ < 2.10e6
G["sigungu_cd"] = G.sigungu_cd.astype(int)
CW = pd.read_csv(PAN, encoding="utf-8-sig", usecols=["fbd_region", "sigungu_cd"]).drop_duplicates(); CW["fbd_region"] = CW.fbd_region.map(NFC); CW["sigungu_cd"] = CW.sigungu_cd.astype(int)
assert CW.fbd_region.is_unique and CW.sigungu_cd.is_unique and len(CW) == 229
assert set(CW.sigungu_cd) == set(G.sigungu_cd), "지도 코드 집합 ≠ 분석 코드 집합"
FO = pd.read_csv(FOF, encoding="utf-8-sig"); FO["fbd_region"] = FO.fbd_region.map(NFC)
sr = FO.groupby("fbd_region").SR.mean().rename("v").reset_index(); ms = FO[FO.year.isin([2018, 2019])].groupby("fbd_region").meas_SR.mean().rename("v").reset_index()
def ld(f):
    r = pd.read_csv(os.path.join(RES, f)); r["fbd_region"] = r.fbd_region.map(NFC); return r.rename(columns={"OR": "v"})
def coded(r):
    r = r.merge(CW, on="fbd_region", how="left", validate="one_to_one"); assert r.sigungu_cd.notna().all(); return r
P = [("A", sr, "seq", "Self-reported obesity, age-sex-standardized %, 2011–2024 mean"), ("B", ld("S3_spatial_OR.csv"), "div", "Spatial odds ratio, self-reported obesity"),
     ("C", ms, "seq", "Measured obesity, age-sex-standardized %, 2018–2019 mean"), ("D", ld("S14_measured_spatialOR.csv"), "div", "Spatial odds ratio, measured obesity")]
BRK = {"seq": [-np.inf, 28, 30, 32, 35, np.inf], "div": [-np.inf, 0.90, 0.95, 1.00, 1.05, 1.10, np.inf]}
LAB = {"seq": ["<28", "28–<30", "30–<32", "32–<35", "≥35"], "div": ["<0.90", "0.90–<0.95", "0.95–<1.00", "1.00–<1.05", "1.05–<1.10", "≥1.10"]}
CMP = {"seq": plt.get_cmap("YlOrRd")(np.linspace(0.12, 0.92, 5)), "div": plt.get_cmap("RdBu_r")(np.linspace(0.08, 0.92, 6))}
VER = []
x0, x1, y0, y1 = G.total_bounds[[0, 2, 1, 3]]; padx = (x1 - x0) * 0.02; pady = (y1 - y0) * 0.02
XL, YL = (x0 - padx, x1 + padx), (y0 - pady, y1 + pady); ASP = (YL[1] - YL[0]) / (XL[1] - XL[0])
LM, GAPX = 0.10, 0.25; MW = (FW - 2 * LM - GAPX) / 2; MH = MW * ASP
TOP, LET, CBG, CBH, CBLAB, ROWG, LEGH = 0.06, 0.20, 0.04, 0.42, 0.0, 0.12, 0.40
PH = LET + MH + CBG + CBH + CBLAB
FH = TOP + 2 * PH + ROWG + LEGH
fig = plt.figure(figsize=(FW, FH), facecolor="white")
def A(l, b, w, h): return fig.add_axes([l / FW, b / FH, w / FW, h / FH])
K = 3
for i, (lab, r, kind, title) in enumerate(P):
    col, row = i % 2, i // 2
    left = LM + col * (MW + GAPX); ptop = FH - TOP - row * (PH + ROWG)
    ax = A(left, ptop - LET - MH, MW, MH)
    lax = A(left, ptop - LET - MH - CBG - CBH, MW, CBH); lax.set_axis_off()
    cols = ["sigungu_cd", "v"] + (["class"] if "class" in r else [])
    g = G.merge(coded(r)[cols], on="sigungu_cd", how="left", validate="one_to_one")
    VER.append((lab, len(g), int(g.v.notna().sum()), set(G.sigungu_cd) == set(coded(r).sigungu_cd) or lab == "C"))
    miss = g.v.isna()
    cmap = ListedColormap(CMP[kind]); norm = BoundaryNorm(BRK[kind], cmap.N)
    g[~miss].plot(column="v", cmap=cmap, norm=norm, ax=ax, edgecolor="#9a9a9a", linewidth=0.08)
    if miss.any(): g[miss].plot(ax=ax, color="#d9d9d9", edgecolor="#6f6f6f", linewidth=0.08, hatch="//////")
    if "class" in g:
        g[g["class"] == "hotspot"].boundary.plot(ax=ax, color="black", linewidth=0.55)
        g[g["class"] == "coldspot"].boundary.plot(ax=ax, color="black", linewidth=0.55, linestyle=(0, (2, 1.2)))
    ax.set_xlim(XL); ax.set_ylim(YL); ax.set_aspect("equal"); ax.set_axis_off()
    ax_, ay_, ah, aw = 0.92, 0.70, 0.075 * K * 0.78, 0.012 * K * 1.4
    ax.add_patch(Polygon([(ax_ - aw / 2, ay_), (ax_ + aw / 2, ay_), (ax_, ay_ + ah)], closed=True, transform=ax.transAxes, facecolor="black", edgecolor="black", lw=0.4))
    ax.text(ax_, ay_ + ah + 0.01, "N", transform=ax.transAxes, ha="center", va="bottom", fontsize=FS * K, fontweight="bold")
    sx, sy, bh = XL[0] + 0.04 * (XL[1] - XL[0]), YL[0] + 0.40 * (YL[1] - YL[0]), 0.012 * (YL[1] - YL[0])
    for k, fc in enumerate(["black", "white"]):
        ax.add_patch(Rectangle((sx + k * 50_000, sy), 50_000, bh, facecolor=fc, edgecolor="black", lw=0.5))
    for d in (0, 50, 100):
        ax.text(sx + d * 1000, sy + bh * 1.6, str(d), ha="center", va="bottom", fontsize=FS)
    ax.text(sx + 50_000, sy - bh * 0.8, "km", ha="center", va="top", fontsize=FS)
    hs = [Patch(facecolor=c, edgecolor="#555", lw=0.4, label=t) for c, t in zip(CMP[kind], LAB[kind])]
    lax.legend(handles=hs, loc="upper center", ncol=3, fontsize=FS, frameon=False, title=title, title_fontsize=FS, handlelength=1.4, handleheight=1.0, columnspacing=0.9, handletextpad=0.4, borderaxespad=0)
    fig.text((left - 0.05) / FW, (ptop - 0.02) / FH, lab, fontsize=PANEL_PT, fontweight="bold", ha="left", va="top")
fig.legend(handles=[Line2D([0], [0], color="black", lw=0.9, label="Hotspot, posterior probability > .95"),
                    Line2D([0], [0], color="black", lw=0.9, ls=(0, (2, 1.2)), label="Coldspot, posterior probability < .05"),
                    Patch(facecolor="#d9d9d9", edgecolor="#6f6f6f", hatch="//////", lw=0.5, label="No measured data (Sejong)")],
           loc="lower center", ncol=3, fontsize=FS, frameon=False, bbox_to_anchor=(0.5, 0.04 / FH), handlelength=2.2, columnspacing=1.2, handletextpad=0.5)
print("[99O4] 코드 결합 검증(패널, 지도 행, 값 있는 행, 코드 집합 일치):", VER)
assert all(n == 229 for _, n, _, _ in VER) and [k for _, _, k, _ in VER][:3] == [229, 229, 228] and all(ok for *_, ok in VER)
tmp = os.path.join(OUTD, "_Figure1_rgba.tiff")
fig.savefig(tmp, dpi=600, facecolor="white")
im = Image.open(tmp).convert("RGB"); im.save(os.path.join(OUTD, "Figure1.tiff"), compression="tiff_lzw", dpi=(600, 600)); os.remove(tmp)
fig.savefig(os.path.join(OUTD, "Figure1.png"), dpi=300, facecolor="white")
print("[99O4] 저장", OUTD, im.size, f"{FW / MM:.0f} x {FH / MM:.0f} mm")
