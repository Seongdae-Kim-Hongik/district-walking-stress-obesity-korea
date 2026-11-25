# Full scan for year-to-year jumps and zero cells on the raw scale, 2011-2024: (a) <0.2x or >5x the neighbour-year mean, (b) single-year 3x spikes, (c) 3x steps that persist, (d) end years.
# Environment: ST_SQLITE (district-year SQLite, read-only), PANEL_WORKDIR (work folder with input_first/, src_data/, tables/, out/, gate/; default current folder).
import sqlite3, json, os, unicodedata
import numpy as np, pandas as pd
B = os.environ.get("PANEL_WORKDIR", os.getcwd())
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
GENGU = ["고양", "성남", "수원", "안산", "안양", "용인", "창원", "포항", "전주", "천안", "청주"]
R["gengu_city"] = R.fbd_region.apply(lambda r: any(g in r for g in GENGU) and r.endswith("시"))
flags, zeros, summ = [], [], []
for name, s in LEVEL.items():
    df = pd.DataFrame({"sigungu_cd": R.sigungu_cd, "fbd_region": R.fbd_region, "year": R.year, "v": s.values, "gengu": R.gengu_city}).sort_values(["sigungu_cd", "year"])
    nflag = 0
    for cd, g in df.groupby("sigungu_cd"):
        g = g[g.v.notna()].reset_index(drop=True)
        for i in range(len(g)):
            v = g.v[i]; prv = g.v[i - 1] if i > 0 else np.nan; nxt = g.v[i + 1] if i < len(g) - 1 else np.nan
            yp = g.year[i - 1] if i > 0 else None; yn = g.year[i + 1] if i < len(g) - 1 else None
            if v == 0:
                zeros.append(dict(variable=name, sigungu_cd=cd, fbd_region=g.fbd_region[i], year=int(g.year[i]), prev=prv, next=nxt,
                                  zero_type=("구조적 0 후보(앞뒤도 0 또는 끝단 0 연속)" if ((prv == 0 or np.isnan(prv)) and (nxt == 0 or np.isnan(nxt))) else "0 입력 결측 후보(앞뒤 양수 사이 단독 0)" if (prv > 0 and nxt > 0) else "경계(한쪽만 0)")))
                continue
            rule = []
            if not np.isnan(prv) and not np.isnan(nxt):
                nb = (prv + nxt) / 2
                if nb > 0 and (v < 0.2 * nb or v > 5 * nb): rule.append("a_이웃평균_0.2배미만_5배초과")
                if prv > 0 and nxt > 0 and ((v >= 3 * prv and v >= 3 * nxt) or (v <= prv / 3 and v <= nxt / 3)): rule.append("b_단년급변_3배")
                if prv > 0 and (v >= 3 * prv or v <= prv / 3) and nxt > 0 and (nxt >= 3 * prv or nxt <= prv / 3) and "b_단년급변_3배" not in rule: rule.append("c_계단_3배_유지")
            else:
                nb1 = prv if not np.isnan(prv) else nxt
                if nb1 and nb1 > 0 and (v >= 3 * nb1 or v <= nb1 / 3): rule.append("d_끝단_한쪽이웃_3배")
            if rule:
                nflag += 1
                flags.append(dict(variable=name, sigungu_cd=cd, fbd_region=g.fbd_region[i], year=int(g.year[i]), value=v, prev_year=yp, prev=prv, next_year=yn, next=nxt,
                                  ratio_prev=v / prv if prv else np.nan, ratio_next=v / nxt if nxt else np.nan, rules=";".join(rule), gengu_city=bool(g.gengu[i])))
    summ.append(dict(variable=name, n_obs=int(df.v.notna().sum()), n_missing=int(df.v.isna().sum()), n_zero=int((df.v == 0).sum()), n_flag=nflag))
F = pd.DataFrame(flags); Z = pd.DataFrame(zeros); S = pd.DataFrame(summ)
S = S.merge(F.groupby("variable").rules.apply(lambda x: ";".join(sorted(set(";".join(x).split(";"))))).rename("rules_hit"), on="variable", how="left") if len(F) else S
F.to_csv(os.path.join(B, "tables/T02_연도급변_전수스캔_플래그.csv"), index=False, encoding="utf-8-sig")
Z.to_csv(os.path.join(B, "tables/T02b_영점칸_분류.csv"), index=False, encoding="utf-8-sig")
S.to_csv(os.path.join(B, "tables/T02c_연도급변_변수요약.csv"), index=False, encoding="utf-8-sig")
G = R[R.gengu_city][["sigungu_cd", "fbd_region", "year"]].copy()
for name, s in LEVEL.items(): G[name] = s[R.gengu_city].values
G.sort_values(["fbd_region", "year"]).to_csv(os.path.join(B, "tables/T03_일반구시11곳_전변수_연도시계열.csv"), index=False, encoding="utf-8-sig")
print(S.to_string()); print("플래그", len(F), "영점", len(Z))
print(F.sort_values(["variable", "year"]).to_string(max_colwidth=40) if len(F) else "")
print(Z.groupby(["variable", "zero_type"]).size().to_string() if len(Z) else "")
print("일반구 시 곳 수", G.sigungu_cd.nunique(), sorted(G.fbd_region.unique()))
