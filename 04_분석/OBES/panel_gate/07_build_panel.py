# Builds the corrected analytic panel with the rules of 80_build_disease_panel_v.py (std3); first rebuilds the first-release panel without corrections and checks it cell by cell, then applies the village-facility corrections and the official basic livelihood security rates.
# Environment: ST_SQLITE (district-year SQLite, read-only), PANEL_WORKDIR (work folder with input_v1/, src_data/, tables/, out/, gate/; default current folder).
import sqlite3, json, os, hashlib, unicodedata, sys
import numpy as np, pandas as pd
B = os.environ.get("PANEL_WORKDIR", os.getcwd()); OUT = os.path.join(B, "out"); os.makedirs(OUT, exist_ok=True)
DB = os.environ["ST_SQLITE"]
NFC = lambda s: unicodedata.normalize("NFC", s) if isinstance(s, str) else s
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
h_db0 = sha(DB); c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
cfg = json.load(open(os.path.join(B, "input_v1/config_OBES_std3.json"))); V = cfg["vars"]; RU = cfg["rules"]
FORCED, CAND, NEG, SENSV = cfg["forced"], cfg["candidates"], cfg["negative_controls"], cfg.get("sensitivity_only_vars", [])
ALLV = FORCED + CAND + SENSV + NEG
P1 = pd.read_csv(os.path.join(B, "input_v1/panel_OBES_std3.csv"), encoding="utf-8-sig")
for k in ["key", "fbd_region", "sido_std", "sigungu"]: P1[k] = P1[k].map(NFC)
srccols = sorted({V[k]["column"] for k in ALLV if not V[k]["column"].startswith("ext__")})
R = pd.read_sql("select year,key,sigungu_cd,boundary_flag," + ",".join(f'"{x}"' for x in srccols) + " from region_year", c)
R.columns = [NFC(x) for x in R.columns]; R["key"] = R.key.map(NFC)
kc = R[["key", "sigungu_cd"]].drop_duplicates(); assert kc.key.is_unique and kc.sigungu_cd.is_unique and len(kc) == 229
D = P1[["year", "key", "sido_std", "sigungu", "fbd_region", "population", "area_m2", "offset", "cases", "E", "SIR", "window", "post"]].merge(kc, on="key", how="left", validate="many_to_one")
D = D.merge(R.drop(columns=["key"]), on=["year", "sigungu_cd"], how="left", validate="one_to_one"); assert len(D) == 2977 and D.sigungu_cd.notna().all()
D = D.sort_values(["key", "year"]).reset_index(drop=True); assert (D[["year", "key"]].values == P1.sort_values(["key", "year"])[["year", "key"]].values).all()
X = pd.read_csv(os.path.join(B, "src_data/EXT_canonical_pool_OBES_261111.csv"), encoding="utf-8-sig"); X["fbd_region"] = X.fbd_region.map(NFC)
D = D.merge(X.rename(columns={k: "x_" + k for k in X.columns if k not in ("year", "fbd_region")}), on=["year", "fbd_region"], how="left", validate="one_to_one")
dp = lambda tag, col: pd.read_sql(f"select year, sigungu_cd, {col} as v from disease_panel where tag=?", c, params=(tag,))
D = D.merge(dp("HTN", "pop80_share_raw").rename(columns={"v": "x_pop80"}), on=["year", "sigungu_cd"], how="left", validate="one_to_one")
D = D.merge(dp("CANM", "basic_livelihood_rate_raw").rename(columns={"v": "x_bl_sqlite"}), on=["year", "sigungu_cd"], how="left", validate="one_to_one")
VIL = pd.read_csv(os.path.join(B, "tables/T07b_마을체육시설_계_v2_전칸.csv"), dtype={"sigungu_cd": str})
D = D.merge(VIL, on=["year", "sigungu_cd"], how="left", validate="one_to_one"); assert (D.vil_n_sqlite.fillna(-1) == D["vil__계"].fillna(-1)).all()
def transform(raw, tr, pop):
    if tr == "none": return raw
    if tr == "per100k": return raw / pop * 1e5
    if tr == "pct_of_pop": return raw / pop * 100
    if tr == "log": return np.log(raw.where(raw > 0))
    if tr == "log1p": return np.log1p(raw.clip(lower=0))
    raise ValueError(tr)
EXTMAP = {"ext__pop80_share": "x_pop80", "ext__restaurant_per10k": "x_restaurant_per10k", "ext__instant_food_per10k": "x_instant_food_per10k",
          "ext__tertiary_educated_share": "x_tertiary_educated_share", "ext__basic_livelihood_rate": "x_basic_livelihood_rate"}
def build(corrected):
    out = D[["year", "key", "sido_std", "sigungu", "fbd_region", "population", "area_m2", "offset", "cases", "E", "SIR", "window", "post"]].copy(); meta = []; fills = []
    for k in ALLV:
        col, tr = V[k]["column"], V[k]["transform"]
        if col.startswith("ext__"):
            src = EXTMAP[col]
            if corrected and col == "ext__basic_livelihood_rate": src = "x_bl_sqlite"
            raw = D[src].astype(float)
        else:
            raw = D[col].astype(float)
            if corrected and col == "vil__계": raw = D["vil_n_v2"].astype(float)
        x = transform(raw, tr, D.population.astype(float))
        t = pd.DataFrame({"key": D.key, "year": D.year, "x": x.values})
        t["x1"] = t.groupby("key")["x"].transform(lambda s: s.interpolate(method="linear", limit_area="inside")); n1 = int(t.x.isna().sum() - t.x1.isna().sum())
        t["x2"] = t.x1.fillna(t.groupby("key")["x1"].transform("mean")); n2 = int(t.x1.isna().sum() - t.x2.isna().sum())
        t["x3"] = t.x2.fillna(t.groupby("year")["x2"].transform("mean")); n3 = int(t.x2.isna().sum() - t.x3.isna().sum()); assert t.x3.notna().all()
        out[k + "_raw"] = t.x3.values; mu, sd = float(t.x3.mean()), float(t.x3.std(ddof=0)); out[k] = (t.x3.values - mu) / sd
        for stage, m in [("interp", t.x.isna() & t.x1.notna()), ("sgg_mean", t.x1.isna() & t.x2.notna()), ("year_mean", t.x2.isna() & t.x3.notna())]:
            for i in np.where(m)[0]: fills.append(dict(variable=k, year=int(t.year.iat[i]), key=t.key.iat[i], stage=stage, filled_raw=float(t.x3.iat[i])))
        yrs = t[t.x.notna()].year
        meta.append(dict(variable=k, source_column=col, transform=tr, label_en=V[k]["label_en"], label_ko=V[k]["label_ko"], role=V[k]["role"], expected_direction=V[k]["expected_direction"],
                         n_obs=len(t), n_missing_raw=int(raw.isna().sum()), n_missing_after_transform=int(x.isna().sum()), coverage=round(float(x.notna().mean()), 3),
                         obs_year_min=int(yrs.min()), obs_year_max=int(yrs.max()), n_sgg_observed=int(t[t.x.notna()].key.nunique()), n_imputed_interp=n1, n_imputed_sgg_mean=n2, n_imputed_year_mean=n3,
                         mean=mu, sd=sd, min=float(t.x3.min()), median=float(t.x3.median()), max=float(t.x3.max()), sgg_all_missing=";".join(sorted(t.loc[t.x2.isna(), "key"].unique()))))
    M = pd.DataFrame(meta); M["low_coverage"] = M.coverage < 0.7
    return out, M, pd.DataFrame(fills)
o0, m0, _ = build(False); P1s = P1.sort_values(["key", "year"]).reset_index(drop=True)
mx = max(float(np.nanmax(np.abs(o0[cname].values - P1s[cname].values))) for cname in o0.columns if o0[cname].dtype != object and cname in P1s)
assert mx < 1e-9, mx; print(f"[07] 재현 검증: 정정 없이 다시 만든 패널 = v1 전 {o0.shape[1]}열·{len(o0)}행, 최대 절대차 {mx:.2e}")
M1 = pd.read_csv(os.path.join(B, "input_v1/panel_META_OBES_std3.csv"), encoding="utf-8-sig")
for col in ["n_missing_raw", "n_imputed_interp", "n_imputed_sgg_mean", "n_imputed_year_mean"]: assert (m0[col].values == M1[col].values).all(), col
o2, m2, F2 = build(True)
def vif_table(df, cols):
    np.seterr(all="ignore"); Xm = np.column_stack([np.ones(len(df)), df[cols].values.astype(np.float64)]); v = {}
    for j, cc in enumerate(cols):
        y = Xm[:, j + 1]; Xo = np.delete(Xm, j + 1, axis=1); b, *_ = np.linalg.lstsq(Xo, y, rcond=None); r2 = 1 - ((y - Xo @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum(); v[cc] = 1 / max(1 - r2, 1e-12)
    return pd.Series(v)
keep = FORCED + CAND; vt0 = vif_table(o2, FORCED + CAND + SENSV); vtf = pd.concat([vif_table(o2, keep), vif_table(o2, keep + SENSV)[SENSV]]) if SENSV else vif_table(o2, keep)
m2["vif_initial"] = m2.variable.map(vt0).round(2); m2["vif_final"] = m2.variable.map(vtf).round(2); m2["vif_dropped"] = False; m2["vif_at_drop"] = np.nan
m2["vif_over_threshold_retained"] = (RU["vif_action"] == "report") & (m2.vif_final > RU["vif_threshold"])
assert m2.vif_final.max() <= 10, "VIF>10"
bf = D[["year", "key", "sigungu_cd", "boundary_flag"]]
o2 = o2.merge(bf, on=["year", "key"], how="left", validate="one_to_one")
BND = {"37310": "군위군: 2023-07-01 경북→대구 편입(전 연도 '대구시|군위군' 키로 통일, 2023 KOSIS 코드 개편 x3xx→x5xx는 옛 코드로 통일)",
       "33010": "청주시: 2014-07-01 청원군 통합(통합 전 연도는 청주+청원 합산, 비율 지표는 청주 값 유지·플래그)",
       "23030": "인천 남구(현 미추홀구): 2018-07-01 명칭 변경(패널 키는 옛 이름 인천시남구, 코드 23030↔23090 별칭)", "29010": "세종시: 2012-07-01 출범(연기군 승계)",
       "38110": "창원시: 2010-07-01 마산·진해 통합(2010 전 연도 비율 결측), 2023 인구 1,009,038 교정", "31280": "여주시: 2013 군→시 승격", "34080": "당진시: 2012 군→시 승격"}
assert set(BND) <= set(o2.sigungu_cd), set(BND) - set(o2.sigungu_cd)
o2["boundary_note"] = o2.sigungu_cd.map(BND).fillna("")
o2["boundary_change_flag"] = (o2.boundary_note != "").astype(int)
H = pd.read_csv(os.path.join(B, "src_data/HTN_Estd2.csv"), encoding="utf-8-sig"); H["fbd_region"] = H.fbd_region.map(NFC)
o2 = o2.merge(H[["year", "fbd_region", "pop65_share_raw"]].rename(columns={"pop65_share_raw": "pop65_share_raw_fbd5"}), on=["year", "fbd_region"], how="left", validate="one_to_one"); assert o2.pop65_share_raw_fbd5.notna().all()
o2 = o2.sort_values(["key", "year"]).reset_index(drop=True)
fp = os.path.join(OUT, "panel_OBES_std3.csv"); o2.to_csv(fp, index=False, encoding="utf-8-sig")
m2.to_csv(os.path.join(OUT, "panel_META_OBES_std3.csv"), index=False, encoding="utf-8-sig")
F2.to_csv(os.path.join(OUT, "fill_cells_OBES_std3.csv"), index=False, encoding="utf-8-sig")
INFO = json.load(open(os.path.join(B, "input_v1/panel_INFO_OBES_std3.json")))
INFO.update(n_imputed_total=int(m2[["n_imputed_interp", "n_imputed_sgg_mean", "n_imputed_year_mean"]].sum().sum()), vif_over_threshold=[v for v in keep if vtf[v] > RU["vif_threshold"]],
            panel_version="v2_261125", panel_v2_note="마을체육시설 계 정정(57칸)·기초생활수급률 공식 관측 복원값(418칸) 반영, 그 밖 v1과 동일",
            source_sqlite_sha256=h_db0)
json.dump(INFO, open(os.path.join(OUT, "panel_INFO_OBES_std3.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
P1s2 = P1s.merge(o2[["year", "key"]], on=["year", "key"])
diff = []
for k in ALLV:
    a, b = P1s[k + "_raw"].values, o2[k + "_raw"].values; za, zb = P1s[k].values, o2[k].values
    nd = int((np.abs(a - b) > 1e-9).sum()); nz = int((np.abs(za - zb) > 1e-9).sum())
    diff.append(dict(variable=k, n_raw_cells_changed=nd, n_z_cells_changed=nz, max_abs_raw_change=float(np.abs(a - b).max()), v1_mean=float(a.mean()), v2_mean=float(b.mean()), v1_sd=float(a.std()), v2_sd=float(b.std()),
                     v1_max=float(a.max()), v2_max=float(b.max()), vif_v1=float(M1.set_index("variable").vif_final.get(k, np.nan)), vif_v2=float(m2.set_index("variable").vif_final.get(k, np.nan))))
DF = pd.DataFrame(diff); DF.to_csv(os.path.join(B, "tables/T08_v1대v2_변수별변화.csv"), index=False, encoding="utf-8-sig")
cells = []
for k in ["village_sports_facilities_per_100k", "basic_livelihood_rate"]:
    m = np.abs(P1s[k + "_raw"].values - o2[k + "_raw"].values) > 1e-9
    for i in np.where(m)[0]: cells.append(dict(variable=k, year=int(o2.year.iat[i]), sigungu_cd=o2.sigungu_cd.iat[i], fbd_region=o2.fbd_region.iat[i], v1_raw=P1s[k + "_raw"].iat[i], v2_raw=o2[k + "_raw"].iat[i], ratio_v2_v1=o2[k + "_raw"].iat[i] / P1s[k + "_raw"].iat[i] if P1s[k + "_raw"].iat[i] else np.nan))
pd.DataFrame(cells).to_csv(os.path.join(B, "tables/T08b_v1대v2_변경칸.csv"), index=False, encoding="utf-8-sig")
assert sha(DB) == h_db0
print(DF[DF.n_raw_cells_changed > 0].to_string()); print("VIF 최대 v2", m2.vif_final.max(), "· 파일", fp, sha(fp))
