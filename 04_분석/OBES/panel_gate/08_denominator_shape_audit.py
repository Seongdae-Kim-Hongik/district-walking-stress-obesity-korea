# Denominator audit (population against KOSIS resident registration tables), boundary-file code, name and area check, 65+ end-fill counts. Also needs SIGUNGU_KEY_DIR and ST_SHP.
# Environment: ST_SQLITE (district-year SQLite, read-only), PANEL_WORKDIR (work folder with input_first/, src_data/, tables/, out/, gate/; default current folder).
import os, sys, json, glob, unicodedata, sqlite3
import numpy as np, pandas as pd
B = os.environ.get("PANEL_WORKDIR", os.getcwd()); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import sigungu_key as sk
NFC = lambda s: unicodedata.normalize("NFC", s) if isinstance(s, str) else s
DB = os.environ["ST_SQLITE"]
c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
P = pd.read_csv(os.path.join(B, "out/panel_OBES_std3.csv"), encoding="utf-8-sig", dtype={"sigungu_cd": str})
RY = pd.read_sql('select year,sigungu_cd,fbd_region,population,"eloc__DT_1YL20631__고령인구비율＜br＞(A÷B×100)" as old65 from region_year where year between 2011 and 2024', c)
RY.columns = [NFC(x) for x in RY.columns]
rep = {}
rows = []
for y in range(2011, 2025):
    f = pd.read_csv(os.path.join(B, f"src_pop/{y}.csv"), encoding="cp949"); f = f[f.iloc[:, 1] == "계"]
    nat_fbd = float(f[f.iloc[:, 0] == "전국"]["전체"].iloc[0])
    s229_ry = float(RY[RY.year == y].population.sum())
    s229_p = float(P[P.year == y].population.sum()) if y != 2019 else np.nan
    rows.append(dict(year=y, fbd_national=nat_fbd, sum229_sqlite=s229_ry, sum229_panel=s229_p, ratio_sqlite_to_fbd=s229_ry / nat_fbd, ratio_panel_to_fbd=(s229_p / nat_fbd) if y != 2019 else np.nan))
D1 = pd.DataFrame(rows)
def kosis(y):
    fp = os.path.join(B, f"kosis/raw/DT_1B04005N_T2_{y}12.json") if y > 2011 else os.path.join(B, "kosis/raw/DT_1B040A3_T20_201112.json")
    return pd.DataFrame(json.load(open(fp))).assign(year=y, src=os.path.basename(fp))
K = pd.concat([kosis(y) for y in range(2011, 2025)], ignore_index=True); K["DT"] = pd.to_numeric(K.DT, errors="coerce")
D1["kosis_national"] = D1.year.map(K[K.C1 == "00"].set_index("year").DT)
D1["ratio_sqlite_to_kosis"] = D1.sum229_sqlite / D1.kosis_national
D1.to_csv(os.path.join(B, "tables/T09_분모_총인구비_229합대전국.csv"), index=False, encoding="utf-8-sig")
assert (D1.ratio_sqlite_to_fbd.round(4) == 1.0).all(), D1[["year", "ratio_sqlite_to_fbd"]]
assert (D1.dropna(subset=["ratio_panel_to_fbd"]).ratio_panel_to_fbd.round(4) == 1.0).all()
assert ((D1.ratio_sqlite_to_kosis - 1).abs() < 0.0005).all(), D1[["year", "ratio_sqlite_to_kosis"]]
print("[08] 총인구 비(패널/FBD) 1.0000 assert 통과, 229 합/KOSIS 전국:", D1.ratio_sqlite_to_kosis.round(6).tolist())
K5 = K[K.C1.str.len() == 5].copy()
K5["sido_cd"] = K5.C1.str[:2]; sidonm = K[K.C1.str.len() == 2].drop_duplicates("C1").set_index("C1").C1_NM.to_dict()
K5["sido_nm"] = K5.sido_cd.map(sidonm); K5["name"] = (K5.sido_nm.fillna("") + " " + K5.C1_NM.str.replace(r"^.*-", "", regex=True)).map(NFC)
kind = pd.Series(sk.name_kind(K5.name.tolist()), index=K5.index); K5["kind"] = kind
K5 = K5[~K5.kind.isin(["일반구(하위구역)", "시도", "전국"])].copy()
bad = K5[~K5.kind.isin(["시군구", "통합전"])]
K5 = K5[K5.kind.isin(["시군구", "통합전"])].copy()
K5["sigungu_cd"] = sk.to_code(K5.name.tolist(), year=K5.year.tolist(), subdistrict="error")
KG = K5.groupby(["year", "sigungu_cd"], as_index=False).DT.sum()
CMP = RY.merge(KG, on=["year", "sigungu_cd"], how="left"); CMP["ratio_sqlite_to_kosis"] = CMP.population / CMP.DT
CMP.to_csv(os.path.join(B, "tables/T10_인구_KOSIS주민등록원표_시군구대조.csv"), index=False, encoding="utf-8-sig")
off = CMP[(CMP.ratio_sqlite_to_kosis - 1).abs() > 0.001]
print("[08] KOSIS 시군구 대조: 칸", len(CMP), "결합 결측", int(CMP.DT.isna().sum()), "|비-1|>0.1%", len(off), "· 미분류 KOSIS 이름", sorted(bad.name.unique())[:10])
print(off.sort_values("ratio_sqlite_to_kosis").to_string() if len(off) else "  (차이 없음)")
rep["kosis_unclassified_names"] = sorted(bad.name.unique())
import geopandas as gpd
G = gpd.read_file(os.path.join(B, "src_shp/final.shp")); G["region"] = G.region.map(NFC)
G["area_geom_m2"] = G.geometry.area
G["sigungu_cd_by_name"] = sk.to_code(G.region.tolist(), subdistrict="error")
T = sk.load_tables(); key = T["key"].set_index("sigungu_cd")
pa = P.drop_duplicates("sigungu_cd")[["sigungu_cd", "fbd_region", "area_m2"]]
S = G[["region", "sigungu_cd_by_name", "area_geom_m2"]].merge(pa, left_on="region", right_on="fbd_region", how="outer", indicator=True)
S["code_match"] = S.sigungu_cd_by_name == S.sigungu_cd
S["key_area_km2"] = S.sigungu_cd.map(key.area_km2.astype(float))
S["panel_over_geom"] = S.area_m2 / S.area_geom_m2; S["geom_over_key_km2"] = S.area_geom_m2 / 1e6 / S.key_area_km2
S["key_name"] = S.sigungu_cd.map(key.canonical_region)
S.to_csv(os.path.join(B, "tables/T11_셰이프_코드이름면적_대조.csv"), index=False, encoding="utf-8-sig")
print("[08] 셰이프:", len(G), "도형 · 결합", S._merge.value_counts().to_dict(), "· 코드 일치", int(S.code_match.sum()), "· 패널면적/기하 범위", S.panel_over_geom.min(), S.panel_over_geom.max(),
      "· 기하/key km² 범위", round(S.geom_over_key_km2.min(), 5), round(S.geom_over_key_km2.max(), 5), "· CRS", G.crs)
assert len(G) == 229 and (S._merge == "both").all() and S.code_match.all()
d = P[["year", "key", "fbd_region", "sigungu_cd"]].merge(RY[["year", "sigungu_cd", "old65"]], on=["year", "sigungu_cd"], how="left").sort_values(["key", "year"])
n_na = int(d.old65.isna().sum()); inside = ends = 0; ends_cells = []
for k, g in d.groupby("key"):
    s = g.old65.values; ok = ~np.isnan(s)
    if ok.sum() == 0: continue
    first, last = np.argmax(ok), len(s) - 1 - np.argmax(ok[::-1])
    for i in np.where(~ok)[0]:
        if first < i < last: inside += 1
        else: ends += 1; ends_cells.append(dict(key=k, year=int(g.year.iat[i]), fill_value=float(s[first] if i < first else s[last])))
pd.DataFrame(ends_cells).to_csv(os.path.join(B, "tables/T12_65세이상_rule2_끝단채움칸.csv"), index=False, encoding="utf-8-sig")
rep["old65"] = dict(n_rows=len(d), n_missing=n_na, interp_inside=inside, rule2_end_fill=ends, regions_end_fill=sorted({e["key"] for e in ends_cells}))
print("[08] 65+ 원열 결측", n_na, "· 안쪽 선형", inside, "· rule=2 끝단 상수 채움", ends, rep["old65"]["regions_end_fill"])
json.dump(rep, open(os.path.join(B, "tables/T12b_감사요약.json"), "w"), ensure_ascii=False, indent=1, default=str)
