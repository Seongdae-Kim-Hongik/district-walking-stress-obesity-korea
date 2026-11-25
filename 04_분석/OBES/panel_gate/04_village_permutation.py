# Village sports facilities: detection of label-value misalignment (row shifts) by matching fingerprints with neighbouring years of other districts in the same province.
# Environment: ST_SQLITE (district-year SQLite, read-only), PANEL_WORKDIR (work folder with input_v1/, src_data/, tables/, out/, gate/; default current folder).
import sqlite3, os, unicodedata, json
import pandas as pd, numpy as np
B = os.environ.get("PANEL_WORKDIR", os.getcwd())
DB = os.environ["ST_SQLITE"]
NFC = lambda s: unicodedata.normalize("NFC", s) if isinstance(s, str) else s
c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
R = pd.read_sql('select year,sido_std,fbd_region,sigungu_cd,population,"vil__계" as n,"vil__시설조성면적_m2" as area from region_year where year between 2009 and 2025', c)
R["fbd_region"] = R.fbd_region.map(NFC); R["sido_std"] = R.sido_std.map(NFC)
R = R[R.n.notna() & R.area.notna()]
fp = {(r.sigungu_cd, r.year): (r.n, r.area) for r in R.itertuples()}
reg = R.drop_duplicates("sigungu_cd").set_index("sigungu_cd")[["fbd_region", "sido_std"]]
rows = []
for r in R.itertuples():
    if not (2011 <= r.year <= 2024): continue
    me = (r.n, r.area); own = [fp.get((r.sigungu_cd, r.year + d)) for d in (-1, 1)]
    if me in own: continue
    if me == (0.0, 0.0): continue
    cand = R[(R.sido_std == r.sido_std) & (R.sigungu_cd != r.sigungu_cd) & (R.year.isin([r.year - 1, r.year + 1])) & (R.n == r.n) & (R.area == r.area)]
    if len(cand) == 0: continue
    qs = sorted(set(cand.sigungu_cd))
    rows.append(dict(year=r.year, sido=r.sido_std, cell_region=r.fbd_region, cell_cd=r.sigungu_cd, n=r.n, area=r.area,
                     own_prev=own[0], own_next=own[1], data_belongs_to=";".join(reg.loc[q, "fbd_region"] for q in qs), data_belongs_to_cd=";".join(qs)))
P = pd.DataFrame(rows).sort_values(["year", "sido", "cell_region"])
P.to_csv(os.path.join(B, "tables/T05_마을체육시설_행밀림_탐지.csv"), index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250)
print(P.to_string(max_colwidth=40)); print(P.groupby(["year", "sido"]).size())
