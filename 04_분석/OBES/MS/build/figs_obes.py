# SUPERSEDED in v1.2.0: earlier figure layout, replaced by figs_obes_print.py (Figure 2) and 01_스크립트/99O4_obes_maps_print.py (Figure 1).
# Figures.
import os, json, shutil
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import content_obes as C
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(os.path.dirname(HERE), "out", "Figures"); os.makedirs(OUT, exist_ok=True)
AN = os.path.normpath(os.path.join(HERE, "..", "..")); N = C.N
from PIL import Image
im = Image.open(os.path.join(AN, "figs", "Fig_maps_OBES.tiff")); im.save(os.path.join(OUT, "Figure1.tiff"), compression="tiff_lzw", dpi=(600, 600))
plt.rcParams.update({"font.family": "Arial", "font.size": 11})
VARS = [v for g, vs in C.GROUPS[:4] for v in vs]
SER = [("Main model", lambda v: N["main"][v], "#1f4e79", "o"), ("Knorr-Held type IV", lambda v: N["kh"][f"KH_IV|{v}"] | {"lo": None}, "#2e8b57", "s"),
       ("Between districts", lambda v: N["bw"][f"{v}|between"], "#b8860b", "^"), ("Within districts", lambda v: N["bw"][f"{v}|within"], "#8b1a1a", "v"),
       ("Measured obesity, 2018-2019", lambda v: N["meas_OR"][v], "#6a3d9a", "D")]
import pandas as pd
kh = pd.read_csv(os.path.join(AN, "results_std3", "S8_KH.csv")); kh = kh[kh.model == "KH_IV"].set_index("cov")
fig, ax = plt.subplots(figsize=(8.5, 9.5)); y0 = np.arange(len(VARS))[::-1]; off = np.linspace(-0.32, 0.32, len(SER))
for (lab, f, col, mk), dy in zip(SER, off):
    for v, y in zip(VARS, y0):
        d = {"OR": kh.loc[v, "OR"], "lo": kh.loc[v, "lo"], "hi": kh.loc[v, "hi"], "cred": bool(kh.loc[v, "cred"])} if lab.startswith("Knorr") else f(v)
        ax.plot([d["lo"], d["hi"]], [y + dy, y + dy], color=col, lw=1.4)
        ax.plot(d["OR"], y + dy, marker=mk, color=col, mfc=col if d["cred"] else "white", ms=6.5, label=lab if v == VARS[0] else None)
ax.axvline(1, color="#555", lw=0.8, ls="--")
ax.set_yticks(y0); ax.set_yticklabels([C.LABEL[v] + (" †" if v in N["robust"] else "") for v in VARS]); ax.set_xscale("log")
ticks = [0.90, 0.95, 1.0, 1.05, 1.10]; ax.set_xticks(ticks); ax.set_xticklabels([f"{t:.2f}" for t in ticks]); ax.minorticks_off()
ax.set_xlabel("Odds ratio per standard deviation (log scale)"); ax.set_ylim(-0.7, len(VARS) - 0.3)
for s in ["top", "right"]: ax.spines[s].set_visible(False)
ax.legend(loc="lower center", bbox_to_anchor=(0.42, -0.2), ncol=2, frameon=False, fontsize=10)
plt.tight_layout()
for ext, dpi in (("png", 300), ("tiff", 600)): fig.savefig(os.path.join(OUT, f"Figure2.{ext}"), dpi=dpi, pil_kwargs={"compression": "tiff_lzw"} if ext == "tiff" else None, bbox_inches="tight")
shutil.copy(os.path.join(AN, "figs", "Fig_maps_OBES.png"), os.path.join(OUT, "Figure1.png"))
print("figs ok")
