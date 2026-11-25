# Zero and year audit, structural zero vs missing classification, fill table, livelihood-rate change log, year-anomaly gate for the model variables, rescan after correction.
# Environment: ST_SQLITE (district-year SQLite, read-only), PANEL_WORKDIR (work folder with input_first/, src_data/, tables/, out/, gate/; default current folder).
import os, json, unicodedata, sqlite3
import numpy as np, pandas as pd
B = os.environ.get("PANEL_WORKDIR", os.getcwd())
NFC = lambda s: unicodedata.normalize("NFC", s) if isinstance(s, str) else s
DB = os.environ["ST_SQLITE"]
c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
cfg = json.load(open(os.path.join(B, "input_first/config_OBES_std3.json"))); V = cfg["vars"]
ALLV = cfg["forced"] + cfg["candidates"] + cfg.get("sensitivity_only_vars", []) + cfg["negative_controls"]
cols = sorted({V[k]["column"] for k in ALLV if not V[k]["column"].startswith("ext__")} | {"kchs__DT_H_OBE_OBE__SR"})
R = pd.read_sql("select year,key,fbd_region,sigungu_cd,population," + ",".join(f'"{x}"' for x in cols) + " from region_year where year between 2011 and 2024", c)
R.columns = [NFC(x) for x in R.columns]; R["fbd_region"] = R.fbd_region.map(NFC); R["key"] = R.key.map(NFC)
X = pd.read_csv(os.path.join(B, "src_data/EXT_canonical_pool_OBES.csv"), encoding="utf-8-sig"); X["fbd_region"] = X.fbd_region.map(NFC)
R = R.merge(X.rename(columns={k: "x_" + k for k in X.columns if k not in ("year", "fbd_region")}), on=["year", "fbd_region"], how="left")
for tag, col, nm in [("HTN", "pop80_share_raw", "x_pop80"), ("CANM", "basic_livelihood_rate_raw", "x_bl")]:
    R = R.merge(pd.read_sql(f"select year,sigungu_cd,{col} as {nm} from disease_panel where tag='{tag}'", c), on=["year", "sigungu_cd"], how="left")
VIL = pd.read_csv(os.path.join(B, "tables/T07b_마을체육시설_계_정정_전칸.csv"), dtype={"sigungu_cd": str}); R = R.merge(VIL, on=["year", "sigungu_cd"], how="left")
pop = R.population.astype(float)
EXT = {"ext__pop80_share": "x_pop80", "ext__restaurant_per10k": "x_restaurant_per10k", "ext__instant_food_per10k": "x_instant_food_per10k", "ext__tertiary_educated_share": "x_tertiary_educated_share", "ext__basic_livelihood_rate": "x_bl"}
LV = {}
for k in ALLV:
    col, tr = V[k]["column"], V[k]["transform"]
    s = R[EXT[col]].astype(float) if col.startswith("ext__") else (R["vil_n_corrected"].astype(float) if col == "vil__계" else R[col].astype(float))
    if tr == "per100k": s = s / pop * 1e5
    elif tr == "pct_of_pop": s = s / pop * 100
    LV[k] = s
LV["OUTCOME_SR"] = R["kchs__DT_H_OBE_OBE__SR"].astype(float)
rows = []
for k, s in LV.items():
    for y in range(2011, 2025):
        m = R.year == y; v = s[m]
        rows.append(dict(variable=k, year=y, n=int(m.sum()), n_obs=int(v.notna().sum()), n_missing=int(v.isna().sum()), n_zero=int((v == 0).sum()), all_missing=bool(v.isna().all()), all_zero=bool((v.fillna(0) == 0).all()),
                         mean=float(v.mean()) if v.notna().any() else np.nan, median=float(v.median()) if v.notna().any() else np.nan, min=float(v.min()) if v.notna().any() else np.nan, max=float(v.max()) if v.notna().any() else np.nan))
A = pd.DataFrame(rows); A.to_csv(os.path.join(B, "tables/T13_영점_연도_감사표.csv"), index=False, encoding="utf-8-sig")
g = []
for k in LV:
    a = A[A.variable == k].set_index("year")
    for y in range(2012, 2025):
        p, q = a.at[y - 1, "median"], a.at[y, "median"]
        if pd.notna(p) and pd.notna(q) and p != 0 and (q / p > 2 or q / p < 0.5): g.append(dict(variable=k, year=y, type="중앙값 급변", prev=p, value=q, ratio=round(q / p, 3)))
        if a.at[y - 1, "n_obs"] >= 100 and a.at[y, "n_obs"] < a.at[y - 1, "n_obs"] * 0.5: g.append(dict(variable=k, year=y, type="시군구 수 급감", prev=a.at[y - 1, "n_obs"], value=a.at[y, "n_obs"], ratio=None))
    for y in range(2011, 2025):
        if 0 < a.at[y, "n_obs"] < 100: g.append(dict(variable=k, year=y, type="시군구 수 100 미만", prev=None, value=a.at[y, "n_obs"], ratio=None))
GY = pd.DataFrame(g, columns=["variable", "year", "type", "prev", "value", "ratio"]); GY.to_csv(os.path.join(B, "gate/GATE_year_anomaly_OBES_std3.csv"), index=False, encoding="utf-8-sig")
M70 = pd.read_csv(os.path.join(B, "src_data/MASTER_연도이상_게이트.csv"), encoding="utf-8-sig"); M70["column"] = M70.column.map(NFC)
M70f = M70[M70.column.isin(cols) & M70.year.between(2011, 2024)]; M70f.to_csv(os.path.join(B, "gate/MASTER_연도이상_게이트_OBES변수_발췌.csv"), index=False, encoding="utf-8-sig")
Z = []
for k, s in LV.items():
    df = pd.DataFrame({"key": R.key, "sigungu_cd": R.sigungu_cd, "fbd_region": R.fbd_region, "year": R.year, "v": s}).sort_values(["key", "year"])
    for key_, gg in df.groupby("key"):
        z = gg[gg.v == 0]
        if z.empty: continue
        allzero = (gg.v.dropna() == 0).all()
        for r in z.itertuples():
            prv = gg[(gg.year < r.year) & gg.v.notna()].v; nxt = gg[(gg.year > r.year) & gg.v.notna()].v
            p = prv.iat[-1] if len(prv) else np.nan; q = nxt.iat[0] if len(nxt) else np.nan
            if allzero or ((p == 0 or np.isnan(p)) and (q == 0 or np.isnan(q))): typ = "구조적 0(그 지역 전 기간 또는 연속 0)"
            elif p > 0 and q > 0: typ = "0 입력 결측 의심(앞뒤 양수 사이 단독 0)"
            else: typ = "경계(앞 또는 뒤만 0)"
            Z.append(dict(variable=k, key=key_, sigungu_cd=r.sigungu_cd, fbd_region=r.fbd_region, year=r.year, prev=p, next=q, zero_type=typ))
ZD = pd.DataFrame(Z)
MECH = {"urban_population_share": "도시지역(용도지역) 미지정 군 — 도시지역 인구 0이 실제 0(국토계획법 용도지역 없음)", "clinics_per_100k": "의원 0곳 군(울릉·옹진 등 도서·산간) — 실제 0",
        "eloc__dt_1yl20421e__1___1": "녹지지역 미지정(도시지역 녹지 0) — 실제 0(log1p 변환으로 유지)", "village_sports_facilities_per_100k": "울릉군 2018·2019–24 마을체육시설 0 — 원표 0(앞뒤도 0, 실제 0 판정); 울진 2023 0은 경북 시트 행 밀림(울릉 값) → 정정에서 16으로 재배치",
        "multicultural_marriage_share": "소규모 군의 그해 다문화 혼인 0건 — 실제 0(사건 수 작음, 음성대조 변수)", "fires_per_10000": "화재 0건"}
if len(ZD): ZD["mechanism"] = ZD.variable.map(MECH).fillna("확인 필요")
ZD.to_csv(os.path.join(B, "tables/T14_구조적0_입력결측_분류표.csv"), index=False, encoding="utf-8-sig")
MISS = []
for k, s in LV.items():
    for r in R[s.isna()].itertuples():
        y = r.year
        if k == "OUTCOME_SR" and y == 2019: why = "원표 없음(2019 KCHS 계측 해, 자가보고 표 미공표) — 결과 행 제외"
        elif k in ("public_sports_facilities_per_100k", "village_sports_facilities_per_100k") and y == 2015: why = "문체부 2015년말 기준 표 미공표(전국 결측) — 시군구 안 선형 보간"
        elif k == "moderate_pa_rate" and y == 2018: why = "KCHS 2018 중강도 신체활동 표 없음 — 시군구 안 선형 보간"
        elif k == "village_sports_facilities_per_100k": why = "정정 결측(원표 결함: 이중 집계·행 어긋남 → 결측 처리) — 시군구 안 선형 보간"
        elif r.sigungu_cd == "29010" and y == 2011: why = "세종 출범 전(2012-07) — 시군구 평균 채움"
        elif k == "rear_seatbelt_wearing_rate" and y < 2014: why = "KCHS 뒷좌석 안전벨트 문항 2014년 시작 — 시군구 평균 채움"
        elif k == "rear_seatbelt_wearing_rate": why = "KCHS 문항 격년·결측 — 보간/시군구 평균"
        elif k == "fiscal_independence" and r.key in ("제주도|서귀포시", "제주도|제주시"): why = "제주 행정시는 재정자주도 미산출(전 연도) — 연도 평균 채움"
        else: why = "원표 결측(지역 단위) — 80 규칙 채움"
        MISS.append(dict(variable=k, key=r.key, year=y, mechanism=why))
MD = pd.DataFrame(MISS); MD.to_csv(os.path.join(B, "tables/T15_결측칸_기전분류.csv"), index=False, encoding="utf-8-sig")
F = pd.read_csv(os.path.join(B, "out/fill_cells_OBES_std3.csv"), encoding="utf-8-sig")
FT = F.pivot_table(index="variable", columns="stage", values="year", aggfunc="count", fill_value=0).reindex(columns=["interp", "sgg_mean", "year_mean"], fill_value=0)
FT["total"] = FT.sum(axis=1); FT["ge100"] = FT.total >= 100; FT.to_csv(os.path.join(B, "tables/T16_채움표_변수x단계.csv"), encoding="utf-8-sig")
P_first = pd.read_csv(os.path.join(B, "input_first/panel_OBES_std3.csv"), encoding="utf-8-sig"); P_corr = pd.read_csv(os.path.join(B, "out/panel_OBES_std3.csv"), encoding="utf-8-sig", dtype={"sigungu_cd": str})
P_first["key"] = P_first.key.map(NFC); P_corr["key"] = P_corr.key.map(NFC)
L = P_corr[["year", "key", "sigungu_cd", "fbd_region", "basic_livelihood_rate_raw"]].merge(P_first[["year", "key", "basic_livelihood_rate_raw"]], on=["year", "key"], suffixes=("_corrected", "_first"))
L["ratio_corrected_over_first"] = L.basic_livelihood_rate_raw_corrected / L.basic_livelihood_rate_raw_first
H33 = pd.read_sql("select code as sigungu_cd, year, persons, population as basis_pop, percent from official_recipient_33_observed_basis", c)
for tb, lab in [("official_recipient_33_observed_basis", "33키 공식 관측 복원(2021-23 일반구 시·창원)"), ("official_recipient_historical3_basis", "과거 3키(2017 서귀포·2020 창원·2023 군위)"), ("official_recipient_2018_basis", "2018 공식 원문 교정"), ("official_recipient_2024_2025_basis", "2024 공식 반영")]:
    t = pd.read_sql(f"select code as sigungu_cd, year from {tb}", c); t["basis"] = lab
    L = L.merge(t, on=["sigungu_cd", "year"], how="left", suffixes=("", "_" + tb[-6:]))
bcols = [x for x in L.columns if x.startswith("basis")]; L["basis_all"] = L[bcols].apply(lambda r: "; ".join(sorted({v for v in r if isinstance(v, str)})), axis=1); L = L.drop(columns=bcols)
b21 = pd.read_sql("select code as sigungu_cd, year, original_raw, corrected_raw from bokjiro_2021_analysis_corrections where tag='CANM'", c)
L = L.merge(b21.assign(basis21="2021 복지로 3지역(양산·함양·거창)")[["sigungu_cd", "year", "basis21"]], on=["sigungu_cd", "year"], how="left")
hist = pd.read_sql("select code as sigungu_cd, year, old_raw as collapsed_or_prior_raw from official_recipient_33_observed_history where tag='CANM'", c)
L = L.merge(hist, on=["sigungu_cd", "year"], how="left")
L["changed"] = (L.basic_livelihood_rate_raw_corrected - L.basic_livelihood_rate_raw_first).abs() > 1e-9
L.to_csv(os.path.join(B, "tables/T17_수급률_정정전후비_로그.csv"), index=False, encoding="utf-8-sig")
s = LV["village_sports_facilities_per_100k"]; rs = []
df = pd.DataFrame({"key": R.key, "fbd_region": R.fbd_region, "year": R.year, "v": s}).sort_values(["key", "year"])
for key_, gg in df.groupby("key"):
    gg = gg[gg.v.notna()].reset_index(drop=True)
    for i in range(1, len(gg) - 1):
        p, v, q = gg.v[i - 1], gg.v[i], gg.v[i + 1]
        if p > 0 and q > 0 and ((v >= 3 * p and v >= 3 * q) or (v <= p / 3 and v <= q / 3)): rs.append(dict(fbd_region=gg.fbd_region[i], year=int(gg.year[i]), prev=p, value=v, next=q))
RS = pd.DataFrame(rs); RS.to_csv(os.path.join(B, "tables/T18_마을체육시설_정정후_잔여단년일탈.csv"), index=False, encoding="utf-8-sig")
U = pd.read_csv(os.path.join(B, "src_data/MASTER_병합_unmatched.csv"), encoding="utf-8-sig")
U.to_csv(os.path.join(B, "tables/T19_MASTER_병합_unmatched_전체.csv"), index=False, encoding="utf-8-sig")
out = dict(zero_year_rows=len(A), all_missing_var_years=A[A.all_missing][["variable", "year"]].values.tolist(), year_gate_flags=len(GY), year_gate=GY.values.tolist(),
           master70_flags_29vars=len(M70f), zero_cells=ZD.groupby(["variable", "zero_type"]).size().reset_index().values.tolist() if len(ZD) else [],
           missing_mech=MD.groupby(["variable", "mechanism"]).size().reset_index().values.tolist(), fill_table=FT.reset_index().values.tolist(),
           livelihood_changed=int(L.changed.sum()), livelihood_ratio_range=[float(L.loc[L.changed, "ratio_corrected_over_first"].min()), float(L.loc[L.changed, "ratio_corrected_over_first"].max())],
           livelihood_changed_by_year=L[L.changed].groupby("year").size().to_dict(), village_residual_spikes=RS.values.tolist(),
           unmatched_rows=len(U), unmatched_by_source=U.source.value_counts().to_dict(), unmatched_vil=U[U.source == "vil"].values.tolist(), unmatched_fac=U[U.source == "fac"].values.tolist())
json.dump(out, open(os.path.join(B, "tables/T20_감사요약_09.json"), "w"), ensure_ascii=False, indent=1, default=str)
print(json.dumps({k: v for k, v in out.items() if k not in ("missing_mech",)}, ensure_ascii=False, default=str)[:6000])
