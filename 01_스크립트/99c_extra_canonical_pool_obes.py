#!/bin/bash
# District covariates from the canonical covariate table (benefit receipt with 2021-2023 correction, education, restaurant and instant-food outlet densities), joined by district code.
import os, unicodedata
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, ".."))
ZD = os.environ.get("ST_CANON", os.path.join(os.path.expandvars("$ST_DATA_ROOT"), "기타/R_Studio/FBD_DATA_ZIP/통합자료"))
PAN = os.path.join(ROOT, "02_패널")
C = pd.read_csv(os.path.join(ZD, "01_공변량_229x2011_2024.csv"), encoding="utf-8-sig", dtype={"sigungu_cd": str})
K = pd.read_csv(os.path.join(ZD, "00_시군구_229_대응표.csv"), encoding="utf-8-sig", dtype={"sigungu_cd": str}); assert K.sigungu_cd.is_unique and len(K) == 229
K.loc[K.sigungu_cd == "37310", "region"] = "대구시군위군"
C = C.drop(columns=["region"]).merge(K[["sigungu_cd", "region"]], on="sigungu_cd", how="left", validate="many_to_one"); assert C.region.notna().all()
out = pd.DataFrame({"year": C.year, "fbd_region": C.region.map(NFC),
                    "basic_livelihood_rate": C.basic_livelihood_recipient_rate.astype(float),
                    "tertiary_educated_share": C.tertiary_educated_share.astype(float),
                    "restaurant_per10k": C.lic_general_restaurant_per10k.astype(float),
                    "instant_food_per10k": C.lic_instant_food_per10k.astype(float)})
out = out[(out.year >= 2011) & (out.year <= 2024)].sort_values(["fbd_region", "year"]).copy()
ref = out[out.year.isin([2020, 2024])].groupby("fbd_region").basic_livelihood_rate.mean()
yy = out.year.isin([2021, 2022, 2023])
bad = yy & ((out.basic_livelihood_rate < 0.2 * out.fbd_region.map(ref)) | (out.basic_livelihood_rate > 5 * out.fbd_region.map(ref)) | (out.fbd_region == NFC("경상남도창원시")))
print("[99c-obes] 수급률 결함 칸", int(bad.sum()), sorted(out.loc[bad, "fbd_region"].unique()))
assert bad.sum() == 33 and out.loc[bad, "fbd_region"].nunique() == 11
out.loc[bad, "basic_livelihood_rate"] = np.nan
VARS = ["basic_livelihood_rate", "tertiary_educated_share", "restaurant_per10k", "instant_food_per10k"]
for c in VARS:
    out[c] = out.groupby("fbd_region")[c].transform(lambda s: s.interpolate(method="linear", limit_area="inside"))
    out[c] = out[c].fillna(out.groupby("fbd_region")[c].transform("mean"))
    nb = out.groupby("fbd_region")[c].transform(lambda s: (s.shift(1) + s.shift(-1)) / 2)
    lo = (out[c] < 0.2 * nb).sum(); hi = (out[c] > 5 * nb).sum()
    print(f"[99c-obes] 게이트 {c}: 이웃 평균 20% 미만 {int(lo)}칸, 5배 초과 {int(hi)}칸, 0값 {int((out[c] == 0).sum())}칸")
    med = out.groupby("year")[c].median(); r = (med / med.shift(1)).dropna()
    print(f"[99c-obes]   연도 중앙값 전년비 범위 {r.min():.3f}-{r.max():.3f}")
_r = out.loc[yy, "basic_livelihood_rate"] / out.loc[yy, "fbd_region"].map(ref); assert _r.between(0.2, 5).all(), "수급률 2021-23 정정 후에도 양쪽 기준 위반"
cw = out[(out.fbd_region == NFC("경상남도창원시")) & yy].basic_livelihood_rate; assert cw.between(2.3, 2.4).all(), cw.tolist(); print("[99c-obes] 창원 2021-23 정정값", cw.round(3).tolist())
assert out.groupby("fbd_region").year.nunique().eq(14).all() and out.fbd_region.nunique() == 229 and not out.duplicated(["year", "fbd_region"]).any()
assert out.notna().all().all(), out.isna().sum()
fp = os.path.join(PAN, "EXT_canonical_pool_OBES.csv"); out.to_csv(fp, index=False, encoding="utf-8-sig")
print(f"[99c-obes] 저장 {fp} {out.shape}"); print(out.describe().T[["min", "50%", "max"]].round(2))
