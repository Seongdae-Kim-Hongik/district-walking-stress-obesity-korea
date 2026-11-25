# Cell-by-cell comparison of the first-release panel (input_v1/) with the district-year SQLite for every variable (read-only, no fits).
# Environment: ST_SQLITE (district-year SQLite, read-only), PANEL_WORKDIR (work folder with input_v1/, src_data/, tables/, out/, gate/; default current folder).
import sqlite3, json, os, sys, hashlib, unicodedata
import numpy as np, pandas as pd
B = os.environ.get("PANEL_WORKDIR", os.getcwd())
DB = os.environ["ST_SQLITE"]
NFC = lambda s: unicodedata.normalize("NFC", s) if isinstance(s, str) else s
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
h0 = sha(DB)
c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
cfg = json.load(open(os.path.join(B, "input_v1/config_OBES_std3.json")))
V = cfg["vars"]
P = pd.read_csv(os.path.join(B, "input_v1/panel_OBES_std3.csv"), encoding="utf-8-sig")
for col in ["key", "fbd_region", "sido_std", "sigungu"]: P[col] = P[col].map(NFC)
src_cols = sorted({v["column"] for v in V.values() if not v["column"].startswith("ext__")} | {"kchs__DT_H_OBE_OBE__SR", "kchs__DT_H_OBE_OBE__SR_SE", "kchs__DT_H_OBE_OBE__N", "kchs__DT_H_OBE_OBE__CR", "kchs__DT_H_OBE_OBE__CR_SE"})
q = ",".join(['year', 'key', 'fbd_region', 'sigungu_cd', 'population', 'boundary_flag'] + [f'"{s}"' for s in src_cols])
R = pd.read_sql(f"select {q} from region_year", c)
R.columns = [NFC(x) for x in R.columns]
for col in ["key", "fbd_region"]: R[col] = R[col].map(NFC)
kc = R[["key", "sigungu_cd"]].drop_duplicates()
assert kc.key.is_unique and kc.sigungu_cd.is_unique and len(kc) == 229, "key-code 1:1 아님"
P = P.merge(kc, on="key", how="left", validate="many_to_one"); assert P.sigungu_cd.notna().all()
M = P.merge(R, on=["year", "sigungu_cd"], how="left", suffixes=("", "_db"), validate="one_to_one"); assert len(M) == len(P)
assert (M.fbd_region == M.fbd_region_db).all() and (M.key == M.key_db).all()
rows = []
def cmp(name, a, b, note):
    a = pd.to_numeric(a, errors="coerce").astype(float); b = pd.to_numeric(b, errors="coerce").astype(float)
    both = a.notna() & b.notna(); d = (a - b).abs()
    rel = (d / b.abs().where(b.abs() > 0)).where(both)
    ndiff = int(((d > 1e-6 * np.maximum(1, b.abs())) & both).sum())
    rows.append(dict(variable=name, n_cells=len(a), n_both=int(both.sum()), n_v1_only=int((a.notna() & b.isna()).sum()), n_db_only=int((a.isna() & b.notna()).sum()),
                     n_diff=ndiff, max_abs_diff=float(d[both].max()) if both.any() else np.nan, max_rel_diff=float(rel.max()) if both.any() else np.nan, note=note))
    return d.where(both)
DIFF = []
def keep_diff(name, a, b, d):
    m = d > 1e-6 * np.maximum(1, pd.to_numeric(b, errors="coerce").abs())
    for i in np.where(m.fillna(False))[0]:
        DIFF.append(dict(variable=name, year=int(M.year.iat[i]), sigungu_cd=M.sigungu_cd.iat[i], fbd_region=M.fbd_region.iat[i], v1=float(a.iat[i]), sqlite=float(b.iat[i])))
for nm, a, b, note in [("population", M.population, M.population_db, "region_year.population(창원 2023 교정 반영 여부)"),
                       ("cases=kchs__DT_H_OBE_OBE__SR", M.cases, M["kchs__DT_H_OBE_OBE__SR"], "결과 SR")]:
    d = cmp(nm, a, b, note); keep_diff(nm, a, b, d)
O = pd.read_csv(os.path.join(B, "src_data/OBES_fh_outcomes.csv"), encoding="utf-8-sig"); O["fbd_region"] = O.fbd_region.map(NFC)
O = O.merge(R[["year", "fbd_region", "sigungu_cd"]], on=["year", "fbd_region"], how="left", validate="one_to_one")
MO = O.merge(R, on=["year", "sigungu_cd"], how="left", suffixes=("", "_db"))
for a, b in [("SR", "kchs__DT_H_OBE_OBE__SR"), ("SR_SE", "kchs__DT_H_OBE_OBE__SR_SE"), ("N", "kchs__DT_H_OBE_OBE__N"), ("CR", "kchs__DT_H_OBE_OBE__CR"), ("CR_SE", "kchs__DT_H_OBE_OBE__CR_SE"), ("population", "population_db")]:
    cmp(f"OBES_fh_outcomes.{a}", MO[a], MO[b], "결과 파일(99O) 3,206행 대 region_year")
def transform(raw, tr, pop):
    if tr == "none": return raw
    if tr == "per100k": return raw / pop * 1e5
    if tr == "pct_of_pop": return raw / pop * 100
    if tr == "log": return np.log(raw.where(raw > 0))
    if tr == "log1p": return np.log1p(raw.clip(lower=0))
    raise ValueError(tr)
for k, v in V.items():
    if v["column"].startswith("ext__"): continue
    x = transform(M[v["column"]].astype(float), v["transform"], M.population.astype(float))
    d = cmp(f"{k}_raw", M[k + "_raw"], x, f"{v['column']} ({v['transform']}); SQLite 결측 칸 = v1 대치 칸")
    keep_diff(f"{k}_raw", M[k + "_raw"], x, d)
X = pd.read_csv(os.path.join(B, "src_data/EXT_canonical_pool_OBES_261111.csv"), encoding="utf-8-sig"); X["fbd_region"] = X.fbd_region.map(NFC); X = X.rename(columns={c_: "x_" + c_ for c_ in X.columns if c_ not in ("year", "fbd_region")})
H = pd.read_csv(os.path.join(B, "src_data/HTN_Estd2.csv"), encoding="utf-8-sig"); H["fbd_region"] = H.fbd_region.map(NFC)
MX = M.merge(X, on=["year", "fbd_region"], how="left", validate="one_to_one").merge(H[["year", "fbd_region", "pop80_share_raw", "pop65_share_raw"]].rename(columns={"pop80_share_raw": "h80", "pop65_share_raw": "h65"}), on=["year", "fbd_region"], how="left", validate="one_to_one")
def dp(tag, col):
    t = pd.read_sql(f"select year, sigungu_cd, {col} as v from disease_panel where tag=?", c, params=(tag,))
    return MX[["year", "sigungu_cd"]].merge(t, on=["year", "sigungu_cd"], how="left")["v"]
pairs = [("pop80_share_raw", MX.h80, "HTN_Estd2.csv"), ("pop80_share_raw", dp("HTN", "pop80_share_raw"), "SQLite disease_panel HTN"),
         ("ext__restaurant_per10k_raw", np.log(MX.x_restaurant_per10k), "EXT 풀 log"), ("ext__restaurant_per10k_raw", np.log(dp("CRC", "restaurants_per10k_raw")), "SQLite CRC restaurants_per10k(log)"),
         ("ext__instant_food_per10k_raw", np.log(MX.x_instant_food_per10k), "EXT 풀 log(SQLite 미수록)"),
         ("tertiary_educated_share_raw", MX.x_tertiary_educated_share, "EXT 풀"), ("tertiary_educated_share_raw", dp("CANM", "tertiary_educated_share_raw"), "SQLite CANM"),
         ("basic_livelihood_rate_raw", MX.x_basic_livelihood_rate, "EXT 풀(99c 33칸 보간)"), ("basic_livelihood_rate_raw", dp("CANM", "basic_livelihood_rate_raw"), "SQLite CANM(공식 복지로 복원 현재값)")]
for k, b, note in pairs:
    d = cmp(f"{k} vs {note}", MX[k], b, note); keep_diff(f"{k} vs {note}", MX[k], pd.Series(np.asarray(b, dtype=float)), d)
S = pd.DataFrame(rows); D = pd.DataFrame(DIFF)
S.to_csv(os.path.join(B, "tables/T01_SQLite대조_변수요약.csv"), index=False, encoding="utf-8-sig")
D.to_csv(os.path.join(B, "tables/T01b_SQLite대조_차이칸.csv"), index=False, encoding="utf-8-sig")
r19 = R[R.year == 2019]["kchs__DT_H_OBE_OBE__SR"]
print("SQLite 2019 OBE SR 비결측", int(r19.notna().sum()))
print(S.to_string())
print("차이 칸", len(D)); print(D.groupby("variable").size().to_string() if len(D) else "")
assert sha(DB) == h0, "SQLite 변경됨"
print("SQLite sha256 전후 동일", h0)
