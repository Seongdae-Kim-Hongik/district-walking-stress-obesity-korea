# Village sports facilities: detection of swapped count cells between two districts of the same province.
# Environment: ST_SQLITE (district-year SQLite, read-only), PANEL_WORKDIR (work folder with input_first/, src_data/, tables/, out/, gate/; default current folder).
import sqlite3, os, unicodedata, pandas as pd
B = os.environ.get("PANEL_WORKDIR", os.getcwd())
DB = os.environ["ST_SQLITE"]
NFC = lambda s: unicodedata.normalize("NFC", s) if isinstance(s, str) else s
c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
R = pd.read_sql('select year,sido_std,fbd_region,sigungu_cd,"vil__계" as n from region_year where year between 2010 and 2024', c)
R["fbd_region"] = R.fbd_region.map(NFC); R["sido_std"] = R.sido_std.map(NFC)
W = R.pivot(index="sigungu_cd", columns="year", values="n"); info = R.drop_duplicates("sigungu_cd").set_index("sigungu_cd")
out = []
for t in range(2011, 2025):
    if t not in W or (t - 1) not in W: continue
    prev = t - 1 if W[t - 1].notna().sum() else t - 2
    for sido, g in info.groupby("sido_std"):
        ids = list(g.index)
        for i, r in enumerate(ids):
            for q in ids[i + 1:]:
                a, b, a0, b0 = W.at[r, t], W.at[q, t], W.at[r, prev], W.at[q, prev]
                if pd.isna(a) or pd.isna(b) or pd.isna(a0) or pd.isna(b0): continue
                if a != a0 and a == b0 and b == a0 and min(a, b) >= 10:
                    out.append(dict(year=t, prev_year=prev, sido=sido, r=info.at[r, "fbd_region"], q=info.at[q, "fbd_region"], r_prev=a0, r_t=a, q_prev=b0, q_t=b,
                                    r_next=W.at[r, t + 1] if t + 1 in W else None, q_next=W.at[q, t + 1] if t + 1 in W else None))
O = pd.DataFrame(out); O.to_csv(os.path.join(B, "tables/T06_마을체육시설_계열교차_탐지.csv"), index=False, encoding="utf-8-sig"); print(O.to_string())
