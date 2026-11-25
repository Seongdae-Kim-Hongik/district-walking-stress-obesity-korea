# Village sports facilities: flagged cells checked for row shifts or duplicated rows in the source tables by fingerprint (count, area).
# Environment: ST_SQLITE (district-year SQLite, read-only), PANEL_WORKDIR (work folder with input_v1/, src_data/, tables/, out/, gate/; default current folder).
import sqlite3, os, unicodedata
import pandas as pd, numpy as np
B = os.environ.get("PANEL_WORKDIR", os.getcwd())
DB = os.environ["ST_SQLITE"]
NFC = lambda s: unicodedata.normalize("NFC", s) if isinstance(s, str) else s
c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
R = pd.read_sql('select year,sido_std,fbd_region,sigungu_cd,population,"vil__계" as n,"vil__시설조성면적_m2" as area,"vil__체력단련시설_점" as fit from region_year where year between 2010 and 2025', c)
R["fbd_region"] = R.fbd_region.map(NFC); R["sido_std"] = R.sido_std.map(NFC)
F = pd.read_csv(os.path.join(B, "tables/T02_연도급변_전수스캔_플래그.csv"), dtype={"sigungu_cd": str})
F = F[F.variable == "village_sports_facilities_per_100k"]
out = []
for _, f in F.iterrows():
    me = R[(R.sigungu_cd == f.sigungu_cd) & (R.year == f.year)].iloc[0]
    same = R[(R.sido_std == me.sido_std) & (R.year.isin([f.year - 1, f.year, f.year + 1])) & (R.sigungu_cd != me.sigungu_cd)]
    m = same[(same.n == me.n) & (same.area == me.area)]
    own = R[(R.sigungu_cd == me.sigungu_cd) & (R.year.isin([f.year - 1, f.year + 1]))]
    out.append(dict(fbd_region=me.fbd_region, sigungu_cd=me.sigungu_cd, year=int(f.year), n=me.n, area=me.area, fit=me.fit, rules=f.rules,
                    own_prev=tuple(own[own.year == f.year - 1][["n", "area"]].values.ravel()), own_next=tuple(own[own.year == f.year + 1][["n", "area"]].values.ravel()),
                    match_other=";".join(f"{r.fbd_region}@{r.year}(n={r.n:g},area={r.area:g})" for r in m.itertuples())))
O = pd.DataFrame(out).sort_values(["year", "fbd_region"])
O.to_csv(os.path.join(B, "tables/T04_마을체육시설_행밀림_지문대조.csv"), index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250)
print(O[O.match_other != ""].to_string(max_colwidth=70)); print("지문 일치", int((O.match_other != "").sum()), "/", len(O))
