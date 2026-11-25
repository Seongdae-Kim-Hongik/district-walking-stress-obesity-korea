#!/usr/bin/env python3
# Numeric additions without new fits: E-values, two-sided tail BH q, within-district OR per 1 SD, spatial OR counts, window P10-P90, yearly Moran's I. Environment: OBES_RES, OBES_EXTRA, OBES_PANEL, OBES_META.
import os, numpy as np, pandas as pd
from decimal import Decimal as _D, ROUND_HALF_UP as _RHU
import builtins as _bi
def round(x, nd=None):
    try:
        q = _D(repr(float(x))).quantize(_D(1).scaleb(-(nd or 0)), rounding=_RHU)
    except Exception:
        return _bi.round(x, nd) if nd is not None else _bi.round(x)
    return float(q) if nd else int(q)

OB = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
R3 = os.environ.get("OBES_RES", os.path.join(OB, "results_std3")); OUT = os.environ.get("OBES_EXTRA", os.path.join(OB, "results_extra"))
d = pd.read_csv(os.environ.get("OBES_PANEL", os.path.join(OB, "panel_OBES_std3.csv")), encoding="utf-8-sig")
M = pd.read_csv(os.environ.get("OBES_META", os.path.join(OB, "panel_META_OBES_std3.csv")), encoding="utf-8-sig").set_index("variable")
S3 = pd.read_csv(os.path.join(R3, "S3_main_OR.csv")); H = S3[S3.role == "hypothesis"].copy()
def ev(rr):
    rr = np.where(rr < 1, 1 / rr, rr); return rr + np.sqrt(rr * (rr - 1))
def ev_ci(lo, hi):
    out = []
    for a, b in zip(lo, hi):
        if a <= 1 <= b: out.append(1.0)
        elif a > 1: out.append(float(ev(np.array([a]))[0]))
        else: out.append(float(ev(np.array([b]))[0]))
    return np.array(out)
for tag, f in [("rare", lambda x: x), ("common_sqrt", np.sqrt)]:
    H[f"E_point_{tag}"] = ev(f(H.OR.values)).round(3); H[f"E_ci_{tag}"] = ev_ci(f(H.lo.values), f(H.hi.values)).round(3)
H["p_tail2"] = (2 * np.minimum(H.P_gt1, 1 - H.P_gt1)).clip(upper=1)
p = H.p_tail2.values; o = np.argsort(p); m = len(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
H["q_BH"] = np.minimum(q, 1).round(4); H["p_tail2"] = H.p_tail2.round(4)
H["note_p0"] = np.where(H.p_tail2 == 0, "P_gt1 반올림 0/1 → p<1e-4", "")
H[["cov", "label_en", "OR", "lo", "hi", "cred", "P_gt1", "p_tail2", "q_BH", "note_p0", "E_point_rare", "E_ci_rare", "E_point_common_sqrt", "E_ci_common_sqrt"]].to_csv(os.path.join(OUT, "G1_evalue_BHq.csv"), index=False)
S9 = pd.read_csv(os.path.join(R3, "S9_between_within.csv")); rows = []
for v in H["cov"]:
    x = d[v]; b = x.groupby(d.fbd_region).transform("mean"); w1 = x - b
    for comp, comp_sd in [("between", b.std(ddof=1)), ("within", w1.std(ddof=1))]:
        r = S9[(S9["cov"] == v) & (S9.component == comp)].iloc[0]
        rows.append(dict(cov=v, component=comp, OR_per_totalSD=r.OR, lo=r.lo, hi=r.hi, cred=r.cred, sd_component_in_totalSD=round(comp_sd, 4),
                         OR_per_componentSD=round(np.exp(np.log(r.OR) * comp_sd), 4), lo_c=round(np.exp(np.log(r.lo) * comp_sd), 4), hi_c=round(np.exp(np.log(r.hi) * comp_sd), 4)))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "G2_S9_per_componentSD.csv"), index=False)
SP = pd.read_csv(os.path.join(R3, "S3_spatial_OR.csv"))
g3 = pd.DataFrame([dict(rule="OR >= 1.10", n=int((SP.OR >= 1.10).sum())), dict(rule="OR >= 1.10 & P_gt1 > 0.95", n=int(((SP.OR >= 1.10) & (SP.P_gt1 > 0.95)).sum())),
                   dict(rule="OR <= 1/1.10", n=int((SP.OR <= 1 / 1.10).sum())), dict(rule="OR <= 1/1.10 & P_gt1 < 0.05", n=int(((SP.OR <= 1 / 1.10) & (SP.P_gt1 < 0.05)).sum())),
                   dict(rule="hotspot (P_gt1 > 0.95)", n=int((SP["class"] == "hotspot").sum())), dict(rule="coldspot (P_gt1 < 0.05)", n=int((SP["class"] == "coldspot").sum())),
                   dict(rule="OR range", n=f"{SP.OR.min():.4f}-{SP.OR.max():.4f}")])
g3.to_csv(os.path.join(OUT, "G3_spatial_OR_counts.csv"), index=False)
SP[SP.OR >= 1.10].sort_values("OR", ascending=False).to_csv(os.path.join(OUT, "G3_spatial_OR_ge110_list.csv"), index=False)
INV = {"log": np.exp, "log1p": np.expm1}
win = {"2011-15": (2011, 2015), "2016-19 (2016-18 observed)": (2016, 2019), "2020-24": (2020, 2024)}
rows = []
for v in list(H["cov"]) + list(S3[S3.role == "forced"]["cov"]):
    tr = M.loc[v, "transform"]; raw = d[v + "_raw"]; nat = INV[tr](raw) if tr in INV else raw
    for wn, (a, b) in win.items():
        k = d.year.between(a, b); s = nat[k]
        rows.append(dict(cov=v, label_en=M.loc[v, "label_en"], transform=tr, window=wn, n=int(k.sum()), P10=round(s.quantile(.1), 3), P50=round(s.quantile(.5), 3), P90=round(s.quantile(.9), 3),
                         P90_P10_ratio=round(s.quantile(.9) / s.quantile(.1), 3) if s.quantile(.1) > 0 else np.nan, SD_model_scale=round(raw[k].std(ddof=0), 4), SD_model_scale_pooled=round(M.loc[v, "sd"], 4)))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "G4_window_P10_P90.csv"), index=False)
MO = pd.read_csv(os.path.join(R3, "S6_moran_year.csv")); MO["E_I"] = round(-1 / 228, 4); MO["perm"] = 999
MO.to_csv(os.path.join(OUT, "G5_moran_by_year.csv"), index=False)
print(H[["cov", "OR", "p_tail2", "q_BH", "E_point_common_sqrt", "E_ci_common_sqrt"]].to_string()); print(g3.to_string()); print(pd.DataFrame(rows).query("cov in ['walking_practice_rate','tertiary_educated_share']").to_string())
