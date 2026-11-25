# Figure 2 forest plot at print size (170 mm wide, 8 pt text, OR and 95% CrI column, 600 dpi). Figure 1 is drawn by 01_스크립트/99O4_obes_maps_v2.py.
import os
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
import content_obes as C
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(os.path.dirname(HERE), "out", "Figures"); os.makedirs(OUT, exist_ok=True)
N = C.N; r3 = C.r3
MM = 1 / 25.4; FW, FH = 170 * MM, 228 * MM; FS = 8
plt.rcParams.update({"font.family": "Arial", "font.size": FS, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2.5, "ytick.major.size": 0})
VARS = [v for g, vs in C.GROUPS[:4] for v in vs]
SER = [("Main model", lambda v: N["main"][v], "#1f4e79", "o"), ("Knorr-Held type IV", lambda v: N["kh"][f"KH_IV|{v}"], "#2e8b57", "s"),
       ("Between districts", lambda v: N["bw"][f"{v}|between"], "#b8860b", "^"), ("Within districts", lambda v: N["bw"][f"{v}|within"], "#8b1a1a", "v"),
       ("Measured obesity, 2018-2019", lambda v: N["meas_OR"][v], "#6a3d9a", "D")]
K = len(SER); nV = len(VARS)
ypos = {(i, k): (nV - 1 - i) * K + (K - 1 - k) for i in range(nV) for k in range(K)}
fig = plt.figure(figsize=(FW, FH), facecolor="white")
L, PW, TW = 52 * MM, 70 * MM, 44 * MM
B, T = 23 * MM, 6 * MM
ax = fig.add_axes([L / FW, B / FH, PW / FW, (FH - B - T) / FH])
tx = fig.add_axes([(L + PW + 2 * MM) / FW, B / FH, TW / FW, (FH - B - T) / FH], sharey=ax); tx.set_axis_off()
lo_lim, hi_lim = 0.885, 1.10
for i in range(nV):
    if i % 2 == 0:
        for a in (ax, tx):
            a.axhspan(ypos[(i, K - 1)] - 0.5, ypos[(i, 0)] + 0.5, color="#f0f0f0", lw=0, zorder=0)
for k, (lab, f, col, mk) in enumerate(SER):
    for i, v in enumerate(VARS):
        d = f(v); y = ypos[(i, k)]
        ax.plot([d["lo"], d["hi"]], [y, y], color=col, lw=0.9, zorder=2)
        ax.plot(d["OR"], y, marker=mk, color=col, mfc=col if d["cred"] else "white", mec=col, mew=0.8, ms=3.6, lw=0, zorder=3)
        tx.text(0.0, y, f"{r3(d['OR'])} ({r3(d['lo'])}–{r3(d['hi'])})", color=col, fontsize=FS, va="center", ha="left", transform=tx.get_yaxis_transform())
ax.axvline(1, color="#555", lw=0.6, ls="--", zorder=1)
ax.set_xscale("log"); ticks = [0.90, 0.95, 1.0, 1.05, 1.10]; ax.set_xticks(ticks); ax.set_xticklabels([f"{t:.2f}" for t in ticks]); ax.minorticks_off()
ax.set_xlim(lo_lim, hi_lim)
ax.set_yticks([ypos[(i, 2)] for i in range(nV)]); ax.set_yticklabels([C.LABEL[v] + (" †" if v in N["robust"] else "") for v in VARS], fontsize=FS)
ax.set_ylim(-0.6, nV * K - 0.4)
ax.set_xlabel("Odds ratio per standard deviation (log scale)", fontsize=FS, labelpad=3); ax.tick_params(axis="x", labelsize=FS, pad=2)
for s in ["top", "right"]: ax.spines[s].set_visible(False)
tx.text(0.0, 1.0, "OR (95% CrI)", transform=tx.transAxes, fontsize=FS, fontweight="bold", ha="left", va="bottom")
fig.legend(handles=[Line2D([0], [0], color=col, marker=mk, mfc=col, mec=col, ms=3.6, lw=0.9, label=lab) for lab, _, col, mk in SER],
           loc="lower center", bbox_to_anchor=(0.5, 1 * MM / FH), ncol=3, frameon=False, fontsize=FS, handlelength=2.0, columnspacing=1.5, handletextpad=0.5)
tmp = os.path.join(OUT, "_Figure2_rgba.tiff")
fig.savefig(tmp, dpi=600, facecolor="white")
im = Image.open(tmp).convert("RGB"); im.save(os.path.join(OUT, "Figure2.tiff"), compression="tiff_lzw", dpi=(600, 600)); os.remove(tmp)
fig.savefig(os.path.join(OUT, "Figure2.png"), dpi=300, facecolor="white")
print("figs v2 ok", im.size, f"{FW / MM:.0f} x {FH / MM:.0f} mm")
