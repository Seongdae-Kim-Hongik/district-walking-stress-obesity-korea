#!/usr/bin/env python3
# Inputs of the round-2 additional analyses and checks without fits: road pavement, outside-district care, healthcare institutions and hospitals, urbanicity fifths and rw2 bins,
# 3-year moving averages, extreme values, calendar-year lags with 2019 source check, imputed-cell flags, instant-food counts against the official yearbook, outcome floor check,
# crude vs standardized prevalence by window, and a search for age-specific obesity columns in the SQLite.
# Environment: E2_RUN (run root with 02_패널/ and 04_분석/OBES/ model inputs), E2_OUT (results folder, default 04_분석/OBES/results_extra2). Also PANEL_WORKDIR (data-gate work folder), ST_SQLITE (district-year SQLite), LICENCE_PANEL_DIR (local-licensing yearly counts).
import os, sys, json, hashlib, unicodedata
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
OBD = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.environ.get("E2_OUT", os.path.join(OBD, "results_extra2")); O = OBD; PW = os.environ.get("PANEL_WORKDIR", os.path.join(O, "panel_work"))
RUN = os.environ["E2_RUN"]
T = os.path.join(OUT, "tables"); os.makedirs(T, exist_ok=True)
P = pd.read_csv(os.path.join(RUN, "04_분석/OBES/panel_OBES_std3.csv"), encoding="utf-8-sig"); P["fbd_region"] = P.fbd_region.map(NFC)
META = pd.read_csv(os.path.join(RUN, "04_분석/OBES/panel_META_OBES_std3.csv"), encoding="utf-8-sig").set_index("variable")
CFG = json.load(open(os.path.join(RUN, "04_분석/OBES/config_OBES_std3.json")))
SETS = CFG["sets_after_vif"]; HYP = SETS[CFG["rules"]["main_set"]]; FORCED = CFG["forced"]; MAINV = list(dict.fromkeys(HYP + FORCED))
M = pd.read_parquet(os.path.join(RUN, "02_패널/MASTER_시군구_연도_wide.parquet")); M.columns = [NFC(c) for c in M.columns]; M["fbd_region"] = M.fbd_region.map(NFC)
z = lambda s: (s - s.mean()) / s.std(ddof=0)
log = []
def L(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); log.append(s)
assert len(P) == 229 * 13 and not P.duplicated(["year", "fbd_region"]).any()
X = P[["year", "fbd_region", "sigungu_cd", "population", "urban_population_share_raw"]].copy()
road = [c for c in M.columns if c.startswith(NFC("eloc__DT_1YL20721__도로포장률"))]; assert len(road) == 1, road
need = ["year", "fbd_region", road[0], "nhis__관내전체__진료실인원수", "nhis__관외전체__진료실인원수", "chk__요양기관_계", "chk__요양기관_종합병원", "chk__요양기관_병원"]
Mx = M[need].rename(columns={road[0]: "road_raw0"})
X = X.merge(Mx, on=["year", "fbd_region"], how="left", validate="one_to_one"); assert len(X) == len(P)
X["road_pavement_rate_raw"] = X.road_raw0.astype(float)
X["outside_care_share_raw"] = 100 * X["nhis__관외전체__진료실인원수"] / (X["nhis__관내전체__진료실인원수"] + X["nhis__관외전체__진료실인원수"])
X["care_institutions_per_100k_raw"] = X["chk__요양기관_계"] / X.population * 1e5
X["hospitals_per_100k_raw"] = (X["chk__요양기관_종합병원"].astype(float) + X["chk__요양기관_병원"].astype(float)) / X.population * 1e5
NEWV = ["road_pavement_rate", "outside_care_share", "care_institutions_per_100k", "hospitals_per_100k"]
X = X.sort_values(["fbd_region", "year"]).reset_index(drop=True)
gate = []
for v in NEWV:
    c = v + "_raw"; s = X[c].astype(float)
    prv = X.groupby("fbd_region")[c].shift(1); nxt = X.groupby("fbd_region")[c].shift(-1); nbm = (prv + nxt) / 2
    ruleb = (prv > 0) & (nxt > 0) & (((s >= 3 * prv) & (s >= 3 * nxt)) | ((s <= prv / 3) & (s <= nxt / 3)))
    rulea = (nbm > 0) & ((s < 0.2 * nbm) | (s > 5 * nbm))
    zero_iso = (s == 0) & (prv > 0) & (nxt > 0)
    bad = (ruleb | rulea | zero_iso).fillna(False)
    nmiss0 = int(s.isna().sum()); X.loc[bad, c] = np.nan
    X[c] = X.groupby("fbd_region")[c].transform(lambda q: q.interpolate(limit_direction="both"))
    X[c] = X[c].fillna(X.groupby("year")[c].transform("mean"))
    assert X[c].notna().all(), v
    X[v] = z(X[c])
    gate.append(dict(variable=v, n_missing_raw=nmiss0, n_flag_rule_a=int(rulea.fillna(False).sum()), n_flag_rule_b=int(ruleb.fillna(False).sum()), n_isolated_zero=int(zero_iso.fillna(False).sum()),
                     n_set_missing=int(bad.sum()), n_imputed_total=nmiss0 + int(bad.sum()), mean=round(X[c].mean(), 3), sd=round(X[c].std(ddof=0), 3), min=round(X[c].min(), 3), max=round(X[c].max(), 3)))
G = pd.DataFrame(gate); G.to_csv(os.path.join(T, "access_new_variable_gate.csv"), index=False, encoding="utf-8-sig"); L("[access gate]\n" + G.to_string())
cor = {v: round(X[v + "_raw"].corr(X.urban_population_share_raw, method="spearman"), 3) for v in NEWV}; L("[access] Spearman with urban share", cor)
ur = X.urban_population_share_raw
X["urban_q5"] = pd.qcut(ur.rank(method="first"), 5, labels=False) + 1
X["urban_rw2_bin"] = pd.qcut(ur.rank(method="first"), 20, labels=False) + 1
P2 = P.sort_values(["fbd_region", "year"]).reset_index(drop=True); assert (P2[["year", "fbd_region"]].values == X[["year", "fbd_region"]].values).all()
KCHS = [v for v in MAINV if META.loc[v, "source_column"].startswith("kchs__")]; L("[88g] KCHS covariates", KCHS)
for v in KCHS:
    r = P2.groupby("fbd_region")[v + "_raw"].transform(lambda q: q.rolling(3, center=True, min_periods=2).mean())
    X["ma3_" + v] = z(r)
ext = []
for v in HYP:
    r = P2[v + "_raw"].astype(float).copy()
    q1 = P2.groupby("year")[v + "_raw"].transform(lambda q: q.quantile(.25)); q3 = P2.groupby("year")[v + "_raw"].transform(lambda q: q.quantile(.75)); iqr = q3 - q1
    med = P2.groupby("year")[v + "_raw"].transform("median")
    k = (r < med - 5 * iqr) | (r > med + 5 * iqr); ext.append(dict(variable=v, n_extreme_5IQR=int(k.sum())))
    r[k] = np.nan; r = r.groupby(P2.fbd_region).transform(lambda q: q.interpolate(limit_direction="both"))
    X["xt_" + v] = z(r)
E = pd.DataFrame(ext); E.to_csv(os.path.join(T, "extreme_extreme_5IQR_counts.csv"), index=False); L("[extreme]", E.to_dict("records"))
EXT = pd.read_csv(os.path.join(PW, "src_data/EXT_canonical_pool_OBES.csv"), encoding="utf-8-sig"); EXT["fbd_region"] = EXT.fbd_region.map(NFC)
HT = pd.read_csv(os.path.join(PW, "src_data/HTN_Estd2.csv"), encoding="utf-8-sig"); HT["fbd_region"] = HT.fbd_region.map(NFC)
VIL = pd.read_csv(os.path.join(PW, "tables/T07b_마을체육시설_계_정정_전칸.csv"), dtype={"sigungu_cd": str}).merge(P.drop_duplicates("sigungu_cd")[["sigungu_cd", "fbd_region"]].assign(sigungu_cd=lambda q: q.sigungu_cd.astype(str)), on="sigungu_cd", how="left")
FCv = pd.read_csv(os.path.join(PW, "out/fill_cells_OBES_std3.csv"), encoding="utf-8-sig"); FCv["fbd_region"] = FCv.key.map(P.drop_duplicates("key").set_index("key").fbd_region).map(NFC)
def raw_from_source(v, yrs):
    sc, tr = META.loc[v, "source_column"], META.loc[v, "transform"]
    if sc.startswith("ext__"):
        k = sc.replace("ext__", "")
        if k == "pop80_share": d = HT[HT.year.isin(yrs)][["year", "fbd_region", "pop80_share_raw"]].rename(columns={"pop80_share_raw": "x"})
        else: d = EXT[EXT.year.isin(yrs)][["year", "fbd_region", k]].rename(columns={k: "x"})
    elif sc == "vil__계":
        d = VIL[VIL.year.isin(yrs)].merge(M[["year", "fbd_region", "population"]], on=["year", "fbd_region"], how="left"); d["x"] = d.vil_n_corrected.astype(float) / d.population * 1e5; d = d[["year", "fbd_region", "x"]]
    else:
        d = M[M.year.isin(yrs)][["year", "fbd_region", sc, "population"]].rename(columns={sc: "x"}); d["x"] = d.x.astype(float)
        if tr == "per100k": d["x"] = d.x / d.population * 1e5
        elif tr == "pct_of_pop": d["x"] = d.x / d.population * 100
        d = d[["year", "fbd_region", "x"]]
    d = d.copy(); d["x"] = d.x.astype(float)
    if tr == "log": d["x"] = np.log(d.x)
    elif tr == "log1p": d["x"] = np.log1p(d.x)
    return d
val, fill19 = [], {}
for v in MAINV:
    a = raw_from_source(v, [2018, 2020]).merge(P[["year", "fbd_region", v + "_raw"]], on=["year", "fbd_region"], how="inner")
    imp = set(map(tuple, FCv[FCv.variable == v][["year", "fbd_region"]].values)); a = a[[(y, r) not in imp for y, r in zip(a.year, a.fbd_region)]]
    dif = (a.x - a[v + "_raw"]).abs(); okc = int((dif < 1e-6).sum()); n_ = len(a)
    s19 = raw_from_source(v, [2019]).set_index("fbd_region").x
    ok = okc == n_ and n_ > 0
    pm = P[P.year.isin([2018, 2020])].groupby("fbd_region")[v + "_raw"].mean()
    s19 = s19.reindex(pm.index); good = s19.notna() & np.isfinite(s19)
    fill19[v] = s19.where(good, pm) if ok else pm
    val.append(dict(variable=v, transform=META.loc[v, "transform"], source=META.loc[v, "source_column"], cells_2018_2020=n_, exact_match=okc, max_abs_diff=round(float(dif.max()), 6) if n_ else np.nan,
                    n_2019_source_missing=int(s19.isna().sum() + (~np.isfinite(s19)).sum()) if len(s19) else 229, used_for_2019=("source(검증 일치)" + (f", 결측 {int((~good).sum())}칸은 2018·2020 평균" if (~good).sum() else "")) if ok else "2018·2020 패널 원값 평균"))
V = pd.DataFrame(val); V.to_csv(os.path.join(T, "callag_2019_fill_validation.csv"), index=False, encoding="utf-8-sig"); L("[callag 2019 validation]\n" + V.to_string())
for v in HYP:
    piv = P.pivot(index="fbd_region", columns="year", values=v + "_raw"); piv[2019] = fill19[v].reindex(piv.index)
    mu, sd = P[v + "_raw"].mean(), P[v + "_raw"].std(ddof=0)
    for lag in (1, 2):
        lagA, lagB = [], []
        for _, r in X[["year", "fbd_region"]].iterrows():
            yl = r.year - lag
            b = piv.at[r.fbd_region, yl] if yl in piv.columns else np.nan
            a = b if yl != 2019 else np.nan
            lagA.append(a); lagB.append(b)
        X[f"clA{lag}_{v}"] = (np.array(lagA, float) - mu) / sd; X[f"clB{lag}_{v}"] = (np.array(lagB, float) - mu) / sd
FC = pd.read_csv(os.path.join(PW, "out/fill_cells_OBES_std3.csv"), encoding="utf-8-sig")
key2reg = P.drop_duplicates("key").set_index("key").fbd_region
FC["fbd_region"] = FC.key.map(key2reg).map(NFC); assert FC.fbd_region.notna().all()
FCm = FC[FC.variable.isin(MAINV)]
IMPV = sorted(FCm.variable.unique()); L("[imputed] fill-cell variables in main model", IMPV, "cells", len(FCm))
X["imp_any_main"] = X.set_index(["year", "fbd_region"]).index.isin(pd.MultiIndex.from_frame(FCm[["year", "fbd_region"]])).astype(int)
L("[imputed] district-years with any imputed main-model cell", int(X.imp_any_main.sum()))
json.dump({"imputed_main_vars_fill_cells": IMPV, "pool_interpolated_vars": ["tertiary_educated_share", "basic_livelihood_rate"]}, open(os.path.join(T, "imputed_imputed_vars.json"), "w"), ensure_ascii=False, indent=1)
keep = ["year", "fbd_region"] + NEWV + [v + "_raw" for v in NEWV] + ["urban_q5", "urban_rw2_bin", "imp_any_main"] + [c for c in X.columns if c.startswith(("ma3_", "xt_", "clA", "clB"))]
X[keep].sort_values(["year", "fbd_region"]).to_csv(os.path.join(OUT, "inputs_extra2.csv"), index=False, encoding="utf-8-sig")
L("[save] inputs_extra2.csv", X[keep].shape)
LIC = os.environ["LICENCE_PANEL_DIR"]
OF = pd.read_csv(os.path.join(LIC, "national_yearend_vs_official.csv"), encoding="utf-8-sig")
e = EXT.merge(M[["year", "fbd_region", "population"]], on=["year", "fbd_region"], how="left", validate="one_to_one")
e["n_inst"] = e.instant_food_per10k * e.population / 1e4
nat = e.groupby("year").n_inst.sum().rename("panel_count_reconstructed").reset_index()
lp = pd.read_csv(os.path.join(LIC, "instant_food_yearend_panel_229x2011_2024.csv"), encoding="utf-8-sig")
nat = nat.merge(lp.groupby("year")[["n_open_all", "n_open_excl_short30"]].sum().reset_index(), on="year", how="left").merge(OF[["year", "official_mfds_yearbook"]], on="year", how="left")
nat["ratio_panel_to_official"] = (nat.panel_count_reconstructed / nat.official_mfds_yearbook).round(3)
nat["ratio_panel_to_licence_excl30"] = (nat.panel_count_reconstructed / nat.n_open_excl_short30).round(4)
nat["ratio_panel_to_licence_all"] = (nat.panel_count_reconstructed / nat.n_open_all).round(4)
nat["official_source"] = "MFDS 식품의약품통계연보 2026(제28호) 표 1-1-1-1 즉석판매제조·가공업 열(연말 개소)"
nat.to_csv(os.path.join(T, "instant_instant_food_vs_official.csv"), index=False, encoding="utf-8-sig"); L("[instant]\n" + nat.round(3).to_string())
tb, mb = "kchs__DT_H_OBE_OBE__", "kchs__DT_117075_MEA_OBE__"
Mo = M[(M.year >= 2011) & (M.year <= 2024)]
d = pd.DataFrame({"year": Mo.year.values, "fbd_region": Mo.fbd_region.values})
for it in ["SR", "SR_SE", "N", "CR", "CR_SE"]: d[it] = Mo[tb + it].astype(float).values
d["meas_SR"] = Mo[mb + "SR"].astype(float).values; d["meas_SE"] = Mo[mb + "SR_SE"].astype(float).values
p = d.SR / 100; d["y_logit"] = np.log(p / (1 - p)); d["v_logit"] = (d.SR_SE / 100) ** 2 / (p * (1 - p)) ** 2
pc = d.CR / 100; d["v_logit_cr"] = (d.CR_SE / 100) ** 2 / (pc * (1 - pc)) ** 2
d["n_eff"] = np.round(p * (1 - p) / (d.SR_SE / 100) ** 2); d["y_eff"] = np.round(p * d.n_eff)
pm = d.meas_SR / 100; d["v_logit_meas"] = (d.meas_SE / 100) ** 2 / (pm * (1 - pm)) ** 2
FO = pd.read_csv(os.path.join(RUN, "02_패널/OBES_fh_outcomes.csv"), encoding="utf-8-sig"); FO["fbd_region"] = FO.fbd_region.map(NFC)
c = FO.merge(d, on=["year", "fbd_region"], suffixes=("", "_nofloor"), validate="one_to_one"); obs = c.year != 2019; obsm = c.meas_SR.notna()
rows = []
for col, msk in [("v_logit", obs), ("v_logit_cr", obs), ("n_eff", obs), ("y_eff", obs), ("y_logit", obs), ("v_logit_meas", obsm)]:
    a, b = c.loc[msk, col], c.loc[msk, col + "_nofloor"]; rows.append(dict(column=col, n_cells=int(msk.sum()), n_cells_different=int(((a - b).abs() > 1e-12).sum()), max_abs_diff=float((a - b).abs().max())))
rows.append(dict(column="SR_SE<=0.05 (SE clip binds)", n_cells=int(obs.sum()), n_cells_different=int((c.loc[obs, "SR_SE"] <= 0.05).sum()), max_abs_diff=np.nan))
rows.append(dict(column="CR_SE<=0.05", n_cells=int(obs.sum()), n_cells_different=int((c.loc[obs, "CR_SE"] <= 0.05).sum()), max_abs_diff=np.nan))
rows.append(dict(column="meas_SE<=0.05", n_cells=int(obsm.sum()), n_cells_different=int((c.loc[obsm, "meas_SE"] <= 0.05).sum()), max_abs_diff=np.nan))
rows.append(dict(column="n_eff_unfloored<10 (n_eff floor binds)", n_cells=int(obs.sum()), n_cells_different=int((c.loc[obs, "n_eff_nofloor"] < 10).sum()), max_abs_diff=np.nan))
nofloor = pd.DataFrame(rows); nofloor.to_csv(os.path.join(T, "nofloor_nofloor_vs_canonical_outcomes.csv"), index=False); L("[nofloor]\n" + nofloor.to_string())
L("[nofloor] min SR_SE", c.loc[obs, "SR_SE"].min(), "min CR_SE", c.loc[obs, "CR_SE"].min(), "min meas_SE", c.loc[obsm, "meas_SE"].min(), "min n_eff(no floor)", c.loc[obs, "n_eff_nofloor"].min())
FO2 = FO[FO.year != 2019].merge(P[["year", "fbd_region", "window", "pop65_share_raw_fbd5", "pop80_share_raw"]], on=["year", "fbd_region"], validate="one_to_one")
def q(s): return f"{s.median():.1f} ({s.quantile(.25):.1f}–{s.quantile(.75):.1f})"
rows = []
for w, g in list(FO2.groupby("window")) + [("All years", FO2)]:
    rows.append({"window": w, "district_years": len(g), "crude_prevalence_median_IQR": q(g.CR), "standardized_prevalence_median_IQR": q(g.SR),
                 "crude_popweighted": round(np.average(g.CR, weights=g.population), 1), "standardized_popweighted": round(np.average(g.SR, weights=g.population), 1),
                 "SR_minus_CR_median_IQR": q(g.SR - g.CR), "spearman_CR_SR": round(g.CR.corr(g.SR, method="spearman"), 3),
                 "share65_median_IQR": q(g.pop65_share_raw_fbd5), "share80_median_IQR": q(g.pop80_share_raw)})
T1 = pd.DataFrame(rows); T1.to_csv(os.path.join(T, "crude_crude_vs_standardized_by_window.csv"), index=False, encoding="utf-8-sig"); L("[crude]\n" + T1.to_string())
DB = os.environ["ST_SQLITE"]
import sqlite3
cn = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
cols = [r[1] for r in cn.execute("pragma table_info(region_year)")]
obe = [x for x in cols if "OBE" in x.upper()]
tabs = [r[0] for r in cn.execute("select name from sqlite_master where type='table'")]
age_hits = []
for tname in tabs:
    for r in cn.execute(f"pragma table_info('{tname}')"):
        nm = str(r[1])
        if ("OBE" in nm.upper() or "비만" in nm) and any(k in nm for k in ["AGE", "age", "연령", "_19", "_30", "_40", "_50", "_60", "_65", "_70"]): age_hits.append(f"{tname}.{nm}")
json.dump({"sqlite": DB, "region_year_obesity_columns": obe, "age_specific_obesity_columns_found": age_hits, "n_tables_scanned": len(tabs)}, open(os.path.join(T, "crude_age_specific_source_search.json"), "w"), ensure_ascii=False, indent=1)
L("[crude] SQLite obesity columns", obe, "age-specific hits", age_hits)
open(os.path.join(OUT, "logs", "e2_00_prep.log"), "w").write("\n".join(log) + "\n")
