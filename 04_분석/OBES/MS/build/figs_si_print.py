# Supplementary Figures S1-S4 at print size from the additional-analysis tables (extra2/): window refits (spatial OR and exceedance probability), window grid,
# observed, fitted and residual maps, and all sensitivity specifications. Maps are joined by district code with a one-to-one and set-equality check.
# Environment: OBES_EXTRA2 (results of extra2/, default 04_분석/OBES/results_extra2), OBES_PANEL (analytic panel), OBES_SHP or ST_SHP (boundary file), OBES_NUMBERS, OBES_SI_OUT.
import os, unicodedata
import numpy as np, pandas as pd, geopandas as gpd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import figs_common as C
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.environ.get("OBES_SI_OUT", os.path.join(os.path.dirname(HERE), "out", "SI_figures")); os.makedirs(OUT, exist_ok=True)
X2 = os.environ.get("OBES_EXTRA2", os.path.normpath(os.path.join(HERE, "..", "..", "results_extra2"))); NFC = lambda s: unicodedata.normalize("NFC", s)
SHP = os.environ.get("OBES_SHP", os.environ.get("ST_SHP", "data/boundary/final.shp"))
MM = 1 / 25.4; FW = 170 * MM; FS = 8; PANEL_PT = 12
plt.rcParams.update({"font.family": "Arial", "font.size": FS, "axes.linewidth": 0.6})
G = gpd.read_file(SHP)
if G.crs is None: G = G.set_crs(epsg=5179)
G["sigungu_cd"] = G.sigungu_cd.astype(int)
CW = pd.read_csv(os.environ.get("OBES_PANEL", os.path.normpath(os.path.join(HERE, "..", "..", "panel_OBES_std3.csv"))), encoding="utf-8-sig", usecols=["fbd_region", "sigungu_cd"]).drop_duplicates(); CW["fbd_region"] = CW.fbd_region.map(NFC); CW["sigungu_cd"] = CW.sigungu_cd.astype(int)
assert CW.fbd_region.is_unique and CW.sigungu_cd.is_unique and len(CW) == 229 and set(CW.sigungu_cd) == set(G.sigungu_cd)
SIDO = G.assign(sido=G.sigungu_cd // 1000).dissolve("sido")
def cjoin(d):
    d = d.copy()
    if "sigungu_cd" not in d: d["fbd_region"] = d.fbd_region.map(NFC); d = d.merge(CW, on="fbd_region", how="left", validate="one_to_one")
    d["sigungu_cd"] = d.sigungu_cd.astype(int); assert set(d.sigungu_cd) == set(G.sigungu_cd) and d.sigungu_cd.is_unique
    return G.merge(d.drop(columns=[c for c in ("fbd_region",) if c in d]), on="sigungu_cd", how="left", validate="one_to_one")
def base(ax):
    SIDO.boundary.plot(ax=ax, color="#555555", linewidth=0.25); ax.set_axis_off(); ax.set_aspect("equal")
def letter(ax, s, t=""):
    ax.text(0.0, 1.0, s, transform=ax.transAxes, fontsize=PANEL_PT, fontweight="bold", va="bottom", ha="left")
    if t: ax.text(0.09, 1.0, t, transform=ax.transAxes, fontsize=FS, va="bottom", ha="left")
def classes(ax, g, col, brk, cols):
    g = g.copy(); g["k"] = pd.cut(g[col], brk, right=False, labels=False)
    for k, c in enumerate(cols): g[g.k == k].plot(color=c, ax=ax, edgecolor="#a0a0a0", linewidth=0.05)
def legend(fig, rect, cols, labs, title):
    ax = fig.add_axes(rect); ax.set_axis_off()
    ax.legend(handles=[Patch(facecolor=c, edgecolor="#777777", linewidth=0.4, label=l) for c, l in zip(cols, labs)], loc="center", ncol=len(labs), frameon=False, fontsize=FS, title=title, title_fontsize=FS, handlelength=1.4, columnspacing=1.0)
def save(fig, n):
    p = os.path.join(OUT, f"FigureS{n}.png"); fig.savefig(p, dpi=300, facecolor="white"); plt.close(fig); return p
DIV_B = [-np.inf, 0.90, 0.95, 1.00, 1.05, 1.10, np.inf]; DIV_L = ["<0.90", "0.90–<0.95", "0.95–<1.00", "1.00–<1.05", "1.05–<1.10", "≥1.10"]; DIV_C = plt.get_cmap("RdBu_r")(np.linspace(0.08, 0.92, 6))
SEQ_B = [-np.inf, 28, 30, 32, 35, np.inf]; SEQ_L = ["<28", "28–<30", "30–<32", "32–<35", "≥35"]; SEQ_C = plt.get_cmap("YlOrRd")(np.linspace(0.12, 0.92, 5))
WIN = [("W1_2011-15", "2011–2015"), ("W2_2016-18", "2016–2018"), ("W3_2020-24", "2020–2024")]

W = pd.read_csv(os.path.join(X2, "W85_window_spatialOR.csv"))
PB = [0, .05, .5, .8, .95, 1.0001]; PL = ["<.05", ".05–<.50", ".50–<.80", ".80–<.95", "≥.95"]; PC = ["#2166AC", "#D1E5F0", "#FEE0B6", "#FDAE61", "#B2182B"]
fig, axs = plt.subplots(2, 3, figsize=(FW, 150 * MM)); plt.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.12, wspace=0.02, hspace=0.32)
for j, (w, wl) in enumerate(WIN):
    g = cjoin(W[W.window == w][["fbd_region", "OR", "P_gt1", "class"]])
    classes(axs[0, j], g, "OR", DIV_B, DIV_C)
    g[g["class"] == "hotspot"].boundary.plot(ax=axs[0, j], color="black", linewidth=0.45); g[g["class"] == "coldspot"].boundary.plot(ax=axs[0, j], color="black", linewidth=0.45, linestyle=(0, (2, 1.2)))
    base(axs[0, j]); letter(axs[0, j], "ABC"[j], wl)
    classes(axs[1, j], g, "P_gt1", PB, PC); base(axs[1, j]); letter(axs[1, j], "DEF"[j], wl)
legend(fig, [0.05, 0.535, 0.9, 0.05], DIV_C, DIV_L, "Spatial odds ratio, refitted within window (solid outline, hotspot; dashed, coldspot)")
legend(fig, [0.05, 0.01, 0.9, 0.06], PC, PL, "Posterior probability that the spatial odds ratio exceeds 1")
P1 = save(fig, 1)

D = pd.read_csv(os.path.join(X2, "tables", "wingrid_window_grid_table.csv"))
ORD = [v for g_, vs in C.GROUPS for v in vs]; assert set(ORD) == set(D["cov"])
fig = plt.figure(figsize=(FW, 160 * MM)); ax = fig.add_axes([0.37, 0.02, 0.62, 0.90]); nr = len(ORD)
for i, v in enumerate(ORD):
    y = nr - 1 - i
    for j, (w, wl) in enumerate(WIN):
        r = D[(D["cov"] == v) & (D.window == w)].iloc[0]; cr = bool(r.cred)
        col = ("#d6e4f2" if r.OR < 1 else "#f6d8d4") if cr else "white"
        ax.add_patch(plt.Rectangle((j, y), 1, 1, facecolor=col, edgecolor="#bbbbbb", linewidth=0.4))
        ax.text(j + 0.5, y + 0.5, f"{C.r3(r.OR)} ({C.r3(r.lo)}–{C.r3(r.hi)})", ha="center", va="center", fontsize=FS, fontweight="bold" if cr else "normal")
        if bool(r.interaction_vs_W1_cred): ax.plot(j + 0.93, y + 0.78, marker="^", color="black", ms=3.5)
    if i == len(ORD) - 5: ax.axhline(y + 1, color="black", linewidth=0.8)
ax.set_xlim(0, 3); ax.set_ylim(0, nr); ax.set_xticks([0.5, 1.5, 2.5]); ax.xaxis.tick_top()
_n = D.drop_duplicates("window").set_index("window")
_nob = pd.read_csv(os.path.join(X2, "W85_window_fit.csv")).set_index("window").n_obs
ax.set_xticklabels([f"{wl}\n(n = {int(_nob[w]):,} district-years)" for w, wl in WIN], fontsize=FS)
ax.set_yticks(np.arange(nr) + 0.5); ax.set_yticklabels([C.LABEL[v] for v in ORD][::-1], fontsize=FS); ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
P2 = save(fig, 2)

E = pd.read_csv(os.path.join(X2, "tables", "obsfit_district_obs_fit_resid.csv"))
g = cjoin(E[["sigungu_cd", "obs", "fit", "res"]])
fig, axs = plt.subplots(1, 3, figsize=(FW, 95 * MM)); plt.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.30, wspace=0.03)
classes(axs[0], g, "obs", SEQ_B, SEQ_C); base(axs[0]); letter(axs[0], "A", "Observed")
classes(axs[1], g, "fit", SEQ_B, SEQ_C); base(axs[1]); letter(axs[1], "B", "Fitted")
RB = [-np.inf, -0.05, -0.02, 0.02, 0.05, np.inf]; RL = ["<−0.05", "−0.05–<−0.02", "−0.02–<0.02", "0.02–<0.05", "≥0.05"]; RC = plt.get_cmap("PuOr_r")(np.linspace(0.1, 0.9, 5))
classes(axs[2], g, "res", RB, RC); base(axs[2]); letter(axs[2], "C", "Residual")
legend(fig, [0.0, 0.15, 1.0, 0.09], SEQ_C, SEQ_L, "Panels A and B: self-reported obesity, %, 2011–2024 mean")
legend(fig, [0.0, 0.03, 1.0, 0.09], RC, RL, "Panel C: observed − fitted, logit")
P3 = save(fig, 3)

L = pd.read_csv(os.path.join(X2, "E2_sensitivity_long_all.csv")); R6 = pd.read_csv(os.path.join(X2, "E2_retention_robust6.csv")).set_index("cov")
ROB = [C.W, C.ED, C.ST, C.DR, C.RE, C.SM]; SHN = {C.W: "Walking", C.ED: "Tertiary education", C.ST: "Perceived stress", C.DR: "High-risk drinking", C.RE: "General restaurants", C.SM: "Current smoking"}
fig = plt.figure(figsize=(FW, 100 * MM)); ax = fig.add_axes([0.20, 0.13, 0.62, 0.77]); rng = np.random.default_rng(20261125)
for i, v in enumerate(ROB):
    d = L[L["cov"] == v]; y = len(ROB) - 1 - i; assert len(d) == int(R6.loc[v, "n_specs"])
    same = (d.OR - 1) * (C.M[v]["OR"] - 1) > 0; keep = d.cred.astype(bool) & same; assert int(keep.sum()) == int(R6.loc[v, "kept"])
    jy = y + rng.uniform(-0.25, 0.25, len(d))
    ax.scatter(d.OR[keep], jy[keep.values], s=7, color="#1f4e79", linewidths=0, zorder=2)
    ax.scatter(d.OR[~keep], jy[~keep.values], s=12, facecolors="white", edgecolors="#b2182b", linewidths=0.7, zorder=3)
    ax.plot([C.M[v]["OR"]], [y], marker="D", ms=5, color="black", zorder=4)
    ax.text(1.045, y, f"{int(R6.loc[v, 'kept'])}/{int(R6.loc[v, 'n_specs'])}", va="center", fontsize=FS, transform=ax.transData, clip_on=False)
ax.axvline(1, color="grey", linewidth=0.6); ax.set_yticks(range(len(ROB))); ax.set_yticklabels([SHN[v] for v in ROB][::-1]); ax.set_xlim(0.92, 1.04)
ax.set_xlabel("Odds ratio per SD in each specification"); ax.text(1.045, len(ROB) - 0.4, "Credible,\nn/specifications", fontsize=FS, va="bottom", clip_on=False)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#1f4e79", ms=4, label="95% CrI excludes 1, same direction"), Line2D([], [], marker="o", ls="", mfc="white", mec="#b2182b", ms=4, label="Interval includes 1"),
                   Line2D([], [], marker="D", ls="", color="black", ms=4, label="Main model")], loc="lower left", fontsize=FS, frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
P4 = save(fig, 4)
print("SI figures", [os.path.basename(p) for p in (P1, P2, P3, P4)])
