#!/usr/bin/env python3
# Post hoc tables without new fits: pointwise WAIC difference SE, S17 lowest-WAIC reference, decision boundaries, S11 hotspot Jaccard, self-report vs measured spatial OR bootstrap,
# outcome floor counts, S7 specification log, fit-setting counts. Usage: python3 05_posthoc.py <results_std3> <results_extra> <panel csv> <outcomes csv>
import sys, os, glob, json
import numpy as np, pandas as pd
from decimal import Decimal as _D, ROUND_HALF_UP as _RHU
def rhu(x, nd):
    try: return float(_D(repr(float(x))).quantize(_D(1).scaleb(-nd), rounding=_RHU))
    except Exception: return x
R2, EXR, PANEL, FOUT = sys.argv[1:5]
WD = os.path.join(R2, "waic_local"); WDX = os.path.join(EXR, "waic_local")
def lw(path):
    d = pd.read_csv(path); d["k"] = d.year.astype(str) + "|" + d.fbd_region; return d.set_index("k")
W = {os.path.basename(p)[:-4]: lw(p) for p in glob.glob(os.path.join(WD, "*.csv")) + glob.glob(os.path.join(WDX, "*.csv"))}
def cmp(a, b):
    A, B = W[a], W[b]; k = A.index.intersection(B.index); assert len(k) == len(A) == len(B), (a, b, len(k), len(A), len(B))
    dif = (A.loc[k, "local_waic"] - B.loc[k, "local_waic"]).values; n = len(dif)
    return dict(model=a, baseline=b, n=n, WAIC_model=rhu(A.local_waic.sum(), 1), WAIC_baseline=rhu(B.local_waic.sum(), 1), dWAIC=rhu(dif.sum(), 1), SE_dWAIC=rhu(np.sqrt(n * dif.var(ddof=1)), 1))
def judge(dw, se):
    a = "≤−6 clearly better" if dw <= -6 else ("−6~−4 better" if dw <= -4 else ("±4 no difference" if abs(dw) < 4 else ">+4 worse" if dw > 4 else "+4 worse(boundary)"))
    return a, ("|Δ|>2SE" if abs(dw) > 2 * se else "|Δ|≤2SE")
rows = []
base = "S3_main"
for k in sorted(W):
    if k in (base,) or k.startswith("S17_") or k.startswith("EX_F_MI"): continue
    rows.append(cmp(k, base))
T = pd.DataFrame(rows)
rerun = [r for r in ["EX_orig_rerun", "S26_main_bym2_typeI", "S17_common"] if r in W]
var_rerun = max(abs(W[r].local_waic.sum() - W[base].local_waic.sum()) for r in rerun) if rerun else np.nan
T["rerun_variation_abs"] = rhu(var_rerun, 2)
pit = pd.read_csv(os.path.join(R2, "S3_main_pit_summary.csv")).iloc[0]
T["cpo_fail_frac_main"] = pit.cpo_fail_frac
T[["judgement_pm4_6", "judgement_SE"]] = [judge(a, b) for a, b in zip(T.dWAIC, T.SE_dWAIC)]
T.to_csv(os.path.join(R2, "P1_waic_pointwise_SE.csv"), index=False)
cand = [c for c in ["S3_main", "S8_KH_II", "S8_KH_III", "S8_KH_IV"] if c in W]
best = min(cand, key=lambda c: W[c].local_waic.sum())
r17 = []
for k in sorted(x for x in W if x.startswith("S17_region_")):
    a = cmp(k, "S17_common"); b = cmp(k, best)
    r17.append(dict(covariate=k.replace("S17_region_", ""), WAIC_region=a["WAIC_model"], WAIC_common=a["WAIC_baseline"], dWAIC_vs_common=a["dWAIC"], SE_vs_common=a["SE_dWAIC"],
                    best_baseline=best, WAIC_best=b["WAIC_baseline"], dWAIC_vs_best=b["dWAIC"], SE_vs_best=b["SE_dWAIC"],
                    judgement_vs_common=judge(a["dWAIC"], a["SE_dWAIC"])[0], judgement_vs_best=judge(b["dWAIC"], b["SE_dWAIC"])[0], cpo_fail_frac_main=pit.cpo_fail_frac, rerun_variation_abs=rhu(var_rerun, 2)))
pd.DataFrame(r17).to_csv(os.path.join(R2, "P2_S17_vs_common_and_best.csv"), index=False)
sp = pd.read_csv(os.path.join(R2, "S3_spatial_OR.csv")); L = pd.read_csv(os.path.join(R2, "S11_spatialOR_long.csv"))
jr = []
for s, g in L.groupby("left_out"):
    keep = sp.sido != s; g = g.set_index("fbd_region").loc[sp.fbd_region[keep]]
    for cl in ["hotspot", "coldspot"]:
        a = set(sp.fbd_region[keep & (sp["class"] == cl)]); b = set(g.index[g["class"] == cl]); u = a | b
        jr.append(dict(left_out=s, cls=cl, n_common_regions=int(keep.sum()), n_main=len(a), n_leftout_fit=len(b), n_both=len(a & b), jaccard=rhu(len(a & b) / len(u), 3) if u else np.nan))
pd.DataFrame(jr).to_csv(os.path.join(R2, "P3_S11_hotspot_jaccard.csv"), index=False)
ms = pd.read_csv(os.path.join(R2, "S14_measured_spatialOR.csv")); P = pd.read_csv(PANEL, encoding="utf-8-sig")
FORC = ["urban_population_share", "clinics_per_100k", "sex_ratio", "fiscal_independence", "pop80_share"]
x = np.log(sp.set_index("fbd_region").OR); y = np.log(ms.set_index("fbd_region").OR.loc[x.index]); Z = P[P.year == 2018].set_index("fbd_region")[FORC].loc[x.index]
rk = lambda v: pd.Series(v).rank().values
def partial_spearman(a, b, Zm):
    Zr = np.column_stack([np.ones(len(a))] + [rk(Zm[:, j]) for j in range(Zm.shape[1])]); ra, rb = rk(a), rk(b)
    ea = ra - Zr @ np.linalg.lstsq(Zr, ra, rcond=None)[0]; eb = rb - Zr @ np.linalg.lstsq(Zr, rb, rcond=None)[0]; return np.corrcoef(ea, eb)[0, 1]
rs = pd.Series(x.values).corr(pd.Series(y.values), method="spearman"); pr = partial_spearman(x.values, y.values, Z.values)
rng = np.random.default_rng(226); B = 2000; bs = []; bp = []
for _ in range(B):
    i = rng.integers(0, len(x), len(x)); bs.append(pd.Series(x.values[i]).corr(pd.Series(y.values[i]), method="spearman")); bp.append(partial_spearman(x.values[i], y.values[i], Z.values[i]))
pd.DataFrame([dict(n=len(x), spearman=rhu(rs, 3), boot_lo=rhu(np.quantile(bs, .025), 3), boot_hi=rhu(np.quantile(bs, .975), 3), partial_spearman_forced2018=rhu(pr, 3),
                   partial_boot_lo=rhu(np.quantile(bp, .025), 3), partial_boot_hi=rhu(np.quantile(bp, .975), 3), n_boot=B, seed=226,
                   note="nominal: 지역 독립 가정 재표집(공간 상관 무시), 공간 OR은 사후평균이라 추정 불확실성 미반영")]).to_csv(os.path.join(R2, "P4_selfreport_vs_measured_spatialOR.csv"), index=False)
FO = pd.read_csv(FOUT, encoding="utf-8-sig"); FO = FO[FO.year != 2019]
p_ = FO.SR / 100; neff_raw = np.round(p_ * (1 - p_) / (FO.SR_SE.clip(lower=0.05) / 100) ** 2)
pd.DataFrame([dict(n_cells=len(FO), n_SR_SE_le0=int((FO.SR_SE <= 0).sum()), n_SR_SE_lt005_clipped=int((FO.SR_SE < 0.05).sum()), n_CR_SE_lt005_clipped=int((FO.CR_SE < 0.05).sum()),
                   n_neff_raw_lt10_floored=int((neff_raw < 10).sum()), neff_min_raw=float(neff_raw.min()), neff_median=float(FO.n_eff.median()),
                   n_meas_SE_lt005=int((FO.meas_SE.dropna() < 0.05).sum()))]).to_csv(os.path.join(R2, "P5_outcome_floor_counts.csv"), index=False)
sl = pd.concat([pd.read_csv(p) for p in glob.glob(os.path.join(R2, "S7_speclog_part_*.csv"))], ignore_index=True).sort_values("spec")
sl.to_csv(os.path.join(R2, "S7_speclog_all.csv"), index=False)
pd.DataFrame([dict(n_specs=len(sl), n_ok=int((sl.status == "ok").sum()), n_failed=int((sl.status != "ok").sum()), n_fallback=int((sl.fit_setting != "default").sum()),
                   sec_median=rhu(sl.seconds.median(), 1), sec_max=rhu(sl.seconds.max(), 1), total_min=rhu(sl.seconds.sum() / 60, 1), per_300_specs_min=rhu(sl.seconds.mean() * 300 / 60, 1), timeout_s=600,
                   seeds=f"20261019+spec ({sl.seed.min()}–{sl.seed.max()})")]).to_csv(os.path.join(R2, "P6_S7_speclog_summary.csv"), index=False)
fs = []
for p in sorted(glob.glob(os.path.join(R2, "*.csv")) + glob.glob(os.path.join(EXR, "fits", "*_OR.csv"))):
    try: d = pd.read_csv(p)
    except Exception: continue
    if "fit_setting" in d:
        for k, v in d.fit_setting.value_counts().items(): fs.append(dict(file=os.path.relpath(p, os.path.dirname(R2)), fit_setting=k, n_rows=int(v)))
pd.DataFrame(fs).to_csv(os.path.join(R2, "P7_fit_setting_counts.csv"), index=False)
print(T.to_string()); print(pd.DataFrame(r17).to_string()); print(pd.DataFrame(jr).groupby("cls").jaccard.describe()); print(open(os.path.join(R2, "P4_selfreport_vs_measured_spatialOR.csv")).read())
