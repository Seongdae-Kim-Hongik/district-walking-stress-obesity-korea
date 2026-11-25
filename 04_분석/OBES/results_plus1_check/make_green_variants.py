# Rebuilds the green-zone covariate under alternative zero handling (log(x+1), log(x+c) with c = half the smallest positive value, zero as missing) with the fill and standardization rules of 80_build_disease_panel_v.py. Usage: python3 make_green_variants.py <project root>
import sys, unicodedata, numpy as np, pandas as pd, json
X = sys.argv[1]; NFC = lambda s: unicodedata.normalize("NFC", s)
P = pd.read_csv(f"{X}/04_분석/OBES/panel_OBES_std3.csv", encoding="utf-8-sig"); P["fbd_region"] = P.fbd_region.map(NFC)
M = pd.read_csv(f"{X}/02_패널/green_raw_from_master.csv", encoding="utf-8-sig"); M.columns = [NFC(c) for c in M.columns]; M["fbd_region"] = M.fbd_region.map(NFC)
src = NFC("eloc__DT_1YL20421E__1인당도시지역면적현황__1인당면적녹지지역")
D = P[["year", "key", "fbd_region"]].merge(M[["year", "fbd_region", src]], on=["year", "fbd_region"], how="left", validate="one_to_one")
D = D.sort_values(["key", "year"]).reset_index(drop=True)
raw = D[src].astype(float)
pos_min = float(raw[raw > 0].min()); c = pos_min / 2
def fill_std(x):
    t = pd.DataFrame({"key": D.key, "year": D.year, "x": x.values})
    t["x1"] = t.groupby("key")["x"].transform(lambda s: s.interpolate(method="linear", limit_area="inside"))
    t["x2"] = t.x1.fillna(t.groupby("key")["x1"].transform("mean"))
    t["x3"] = t.x2.fillna(t.groupby("year")["x2"].transform("mean"))
    n = dict(missing=int(t.x.isna().sum()), interp=int(t.x.isna().sum() - t.x1.isna().sum()), sgg=int(t.x1.isna().sum() - t.x2.isna().sum()), year=int(t.x2.isna().sum() - t.x3.isna().sum()))
    z = (t.x3 - t.x3.mean()) / t.x3.std(ddof=0)
    return z.values, t.x3.values, n
zp1, xp1, np1 = fill_std(np.log1p(raw.clip(lower=0)))
zc, xc, nc = fill_std(np.log(raw.clip(lower=0) + c))
zna, xna, nna = fill_std(np.log(raw.where(raw > 0)))
D["green_p1"], D["green_halfmin"], D["green_zeroNA"] = zp1, zc, zna
D["green_zero_flag"] = (raw == 0).astype(int)
chk = D.merge(P[["year", "fbd_region", "eloc__dt_1yl20421e__1___1"]], on=["year", "fbd_region"])
dev = float((chk.green_p1 - chk["eloc__dt_1yl20421e__1___1"]).abs().max())
info = dict(n_cells=len(D), n_zero=int((raw == 0).sum()), n_missing_raw=int(raw.isna().sum()), positive_min=pos_min, c=c,
            zero_cells=D.loc[raw == 0, ["year", "fbd_region"]].astype(str).agg(" ".join, axis=1).tolist(),
            missing_cells=D.loc[raw.isna(), ["year", "fbd_region"]].astype(str).agg(" ".join, axis=1).tolist(),
            fill_plus1=np1, fill_halfmin=nc, fill_zeroNA=nna, max_abs_dev_rebuilt_vs_panel=dev,
            corr_p1_halfmin=float(np.corrcoef(zp1, zc)[0, 1]), corr_p1_zeroNA=float(np.corrcoef(zp1, zna)[0, 1]),
            sd_log_scale=dict(plus1=float(np.std(xp1)), halfmin=float(np.std(xc)), zeroNA=float(np.std(xna))))
assert dev < 1e-9, dev
D[["year", "fbd_region", "green_p1", "green_halfmin", "green_zeroNA", "green_zero_flag"]].to_csv(f"{X}/04_분석/OBES/results_plus1_check/green_variants.csv", index=False, encoding="utf-8")
json.dump(info, open(f"{X}/04_분석/OBES/results_plus1_check/green_variants_info.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(info, ensure_ascii=False, indent=1))
