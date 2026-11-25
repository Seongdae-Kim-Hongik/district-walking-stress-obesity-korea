#!/usr/bin/env python3
# Fit inputs for the data-gate threshold alternatives: newly flagged main-model cells set to missing, interpolated within district on the raw scale and standardized with the panel mean and SD. Environment: E2_RUN (run root with 02_패널/ and 04_분석/OBES/ model inputs), E2_OUT (results folder, default 04_분석/OBES/results_extra2).
import os, json, numpy as np, pandas as pd, unicodedata
NFC = lambda s: unicodedata.normalize("NFC", str(s))
OBD = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.environ.get("E2_OUT", os.path.join(OBD, "results_extra2")); RUN = os.environ["E2_RUN"]; T = os.path.join(OUT, "tables")
P = pd.read_csv(os.path.join(RUN, "04_분석/OBES/panel_OBES_std3.csv"), encoding="utf-8-sig"); P["fbd_region"] = P.fbd_region.map(NFC); P["sigungu_cd"] = P.sigungu_cd.astype(str)
CFG = json.load(open(os.path.join(RUN, "04_분석/OBES/config_OBES_std3.json"))); MAINV = list(dict.fromkeys(CFG["sets_after_vif"][CFG["rules"]["main_set"]] + CFG["forced"]))
N = pd.read_csv(os.path.join(T, "gatethr_scan_new_cells_by_alternative.csv"), dtype={"sigungu_cd": str})
v25 = pd.read_csv(os.path.join(T, "gatethr_village_rules_area0.05_cnt0.25_spk2.5.csv"), dtype={"sigungu_cd": str}); v30 = pd.read_csv(os.path.join(T, "gatethr_village_rules_area0.05_cnt0.25_spk3.0.csv"), dtype={"sigungu_cd": str})
k6 = lambda d: set(zip(d[d.step.str.startswith("S6")].sigungu_cd, d[d.step.str.startswith("S6")].year))
vil_extra = k6(v25) - k6(v30); print("village rule 2.5-fold extra cells:", vil_extra)
X = P[["year", "fbd_region", "sigungu_cd"]].copy(); rep = []
def variant(tag, cells):
    Q = P.copy()
    for v in MAINV:
        cv = {(c, y) for vv, c, y in cells if vv == v}
        if not cv: continue
        k = [(c, y) in cv for c, y in zip(Q.sigungu_cd, Q.year)]
        r = Q[v + "_raw"].astype(float).copy(); r[k] = np.nan
        r = r.groupby(Q.fbd_region).transform(lambda s: s.sort_index().interpolate(limit_direction="both")) if False else r
        tmp = pd.DataFrame({"fbd_region": Q.fbd_region, "year": Q.year, "r": r}).sort_values(["fbd_region", "year"])
        tmp["r"] = tmp.groupby("fbd_region").r.transform(lambda s: s.interpolate(limit_direction="both")); r = tmp.r.reindex(Q.index)
        mu, sd = P[v + "_raw"].mean(), P[v + "_raw"].std(ddof=0); Q[v] = (r - mu) / sd
        rep.append(dict(variant=tag, variable=v, n_cells_set_missing=int(sum(k)), max_abs_z_change=round(float((Q[v] - P[v]).abs().max()), 4)))
    for v in MAINV: X[f"{tag}__{v}"] = Q[v].values
for fold, lo in [(2.5, 0.25), (2.0, 0.333)]:
    sub = N[(N.rule_bcd_fold == fold) & (N.rule_a_low.round(3) == lo)]
    variant(f"a205k{str(fold).replace('.', '')}", list(zip(sub.variable, sub.sigungu_cd, sub.year)))
variant("a205vil25", [("village_sports_facilities_per_100k", c, y) for c, y in vil_extra])
X.drop(columns="sigungu_cd").to_csv(os.path.join(OUT, "inputs_extra2_a205.csv"), index=False, encoding="utf-8-sig")
pd.DataFrame(rep).to_csv(os.path.join(T, "gatethr_fit_input_changes.csv"), index=False); print(pd.DataFrame(rep).to_string())
