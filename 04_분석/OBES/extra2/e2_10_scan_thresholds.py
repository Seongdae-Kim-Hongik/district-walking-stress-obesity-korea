#!/usr/bin/env python3
# Year-jump scan of the data gate (panel_gate/02_year_jump_scan.py) with alternative thresholds as arguments. Read-only. Environment: E2_PANEL_WORKDIR, ST_SQLITE, E2_TABLES.
import sqlite3, json, os, unicodedata
import numpy as np, pandas as pd
B = os.environ["E2_PANEL_WORKDIR"]
DB = os.environ["ST_SQLITE"]
NFC = lambda s: unicodedata.normalize("NFC", s) if isinstance(s, str) else s
c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
cfg = json.load(open(os.path.join(B, "input_first/config_OBES_std3.json"))); V = cfg["vars"]
cols = sorted({v["column"] for v in V.values() if not v["column"].startswith("ext__")} | {"kchs__DT_H_OBE_OBE__SR", "kchs__DT_H_OBE_OBE__SR_SE"})
R = pd.read_sql("select year,key,fbd_region,sigungu_cd,population," + ",".join(f'"{x}"' for x in cols) + " from region_year where year between 2011 and 2024", c)
R.columns = [NFC(x) for x in R.columns]; R["fbd_region"] = R.fbd_region.map(NFC)
X = pd.read_csv(os.path.join(B, "src_data/EXT_canonical_pool_OBES.csv"), encoding="utf-8-sig"); X["fbd_region"] = X.fbd_region.map(NFC)
R = R.merge(X.rename(columns={k: "x_" + k for k in X.columns if k not in ("year", "fbd_region")}), on=["year", "fbd_region"], how="left", validate="one_to_one")
H = pd.read_csv(os.path.join(B, "src_data/HTN_Estd2.csv"), encoding="utf-8-sig"); H["fbd_region"] = H.fbd_region.map(NFC)
R = R.merge(H[["year", "fbd_region", "pop80_share_raw", "pop65_share_raw"]], on=["year", "fbd_region"], how="left", validate="one_to_one")
bl = pd.read_sql("select year,sigungu_cd,basic_livelihood_rate_raw as bl_sqlite from disease_panel where tag='CANM'", c)
R = R.merge(bl, on=["year", "sigungu_cd"], how="left", validate="one_to_one")
pop = R.population.astype(float)
LEVEL = {}
for k, v in V.items():
    col, tr = v["column"], v["transform"]
    if col.startswith("ext__"):
        s = {"ext__pop80_share": R.pop80_share_raw, "ext__restaurant_per10k": R.x_restaurant_per10k, "ext__instant_food_per10k": R.x_instant_food_per10k,
             "ext__tertiary_educated_share": R.x_tertiary_educated_share, "ext__basic_livelihood_rate": R.bl_sqlite}[col]
    else:
        s = R[col].astype(float)
        if tr == "per100k": s = s / pop * 1e5
        elif tr == "pct_of_pop": s = s / pop * 100
    LEVEL[k] = s.astype(float)
LEVEL["basic_livelihood_rate(99c EXT 보간판)"] = R.x_basic_livelihood_rate
LEVEL["OUTCOME_SR"] = R["kchs__DT_H_OBE_OBE__SR"].astype(float); LEVEL["OUTCOME_SR_SE"] = R["kchs__DT_H_OBE_OBE__SR_SE"].astype(float)
LEVEL["population"] = pop; LEVEL["pop65_share(S12 조율 모형용)"] = R.pop65_share_raw

import itertools
OUTT = os.environ["E2_TABLES"]
def scan(alo, ahi, bk):
    flags = []
    for name, s in LEVEL.items():
        df = pd.DataFrame({"sigungu_cd": R.sigungu_cd, "fbd_region": R.fbd_region, "year": R.year, "v": s.values}).sort_values(["sigungu_cd", "year"])
        for cd, g in df.groupby("sigungu_cd"):
            g = g[g.v.notna()].reset_index(drop=True)
            for i in range(len(g)):
                v = g.v[i]; prv = g.v[i - 1] if i > 0 else np.nan; nxt = g.v[i + 1] if i < len(g) - 1 else np.nan
                if v == 0: continue
                rule = []
                if not np.isnan(prv) and not np.isnan(nxt):
                    nb = (prv + nxt) / 2
                    if nb > 0 and (v < alo * nb or v > ahi * nb): rule.append("a")
                    if prv > 0 and nxt > 0 and ((v >= bk * prv and v >= bk * nxt) or (v <= prv / bk and v <= nxt / bk)): rule.append("b")
                    if prv > 0 and (v >= bk * prv or v <= prv / bk) and nxt > 0 and (nxt >= bk * prv or nxt <= prv / bk) and "b" not in rule: rule.append("c")
                else:
                    nb1 = prv if not np.isnan(prv) else nxt
                    if nb1 and nb1 > 0 and (v >= bk * nb1 or v <= nb1 / bk): rule.append("d")
                if rule: flags.append((name, cd, int(g.year[i]), ";".join(rule)))
    return flags
base = scan(0.2, 5, 3)
T02 = pd.read_csv(os.path.join(B, "tables/T02_연도급변_전수스캔_플래그.csv"), encoding="utf-8-sig", dtype={"sigungu_cd": str})
bset = {(a, str(b), c) for a, b, c, _ in base}; tset = {(a, str(b), int(c)) for a, b, c in zip(T02.variable, T02.sigungu_cd, T02.year)}
print("baseline reproduces T02:", bset == tset, len(bset), len(tset))
rows = []; det = []
for alo, ahi, bk in [(0.2, 5, 3), (0.25, 4, 2.5), (1/3, 3, 2), (0.1, 10, 4), (0.1, 10, 5), (0.2, 5, 2), (0.2, 5, 4)]:
    f = scan(alo, ahi, bk); fs = {(a, str(b), c) for a, b, c, _ in f}
    new = fs - bset; gone = bset - fs
    model_vars = [k for k in LEVEL if k in V]
    rows.append(dict(rule_a_low=round(alo, 3), rule_a_high=ahi, rule_bcd_fold=bk, n_flag=len(fs), n_flag_model_vars=sum(1 for x in fs if x[0] in V), n_new_vs_baseline=len(new), n_new_model_vars=sum(1 for x in new if x[0] in V), n_dropped_vs_baseline=len(gone),
                     new_by_var=";".join(f"{k}:{v}" for k, v in pd.Series([x[0] for x in new]).value_counts().items()) if new else ""))
    for x in sorted(new): det.append(dict(rule_a_low=round(alo, 3), rule_bcd_fold=bk, variable=x[0], sigungu_cd=x[1], year=x[2]))
pd.DataFrame(rows).to_csv(os.path.join(OUTT, "gatethr_scan_threshold_alternatives.csv"), index=False, encoding="utf-8-sig")
pd.DataFrame(det).to_csv(os.path.join(OUTT, "gatethr_scan_new_cells_by_alternative.csv"), index=False, encoding="utf-8-sig")
print(pd.DataFrame(rows).drop(columns="new_by_var").to_string()); print(pd.DataFrame(rows).new_by_var.to_string())
