#!/usr/bin/env python3
# Inputs for the additional analyses A-F: imputed-cell flags (from the panel fill ledger), facility densities per km2, boundary-change lists and multiple-imputation noise SDs.
# Usage: python3 82_extra_inputs.py <run root (parent of 01_스크립트)> <fill_cells_OBES_std3.csv> <output folder>
import os, sys, json, unicodedata
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
ROOT, FILL, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
d = pd.read_csv(os.path.join(ROOT, "04_분석", "OBES", "panel_OBES_std3.csv"), encoding="utf-8-sig"); d["fbd_region"] = d.fbd_region.map(NFC); d["key"] = d.key.map(NFC)
assert len(d) == 2977
F = pd.read_csv(FILL, encoding="utf-8-sig"); F["key"] = F.key.map(NFC)
def flag(var, stages=None):
    g = F[(F.variable == var) & (F.stage.isin(stages) if stages else True)]
    s = set(zip(g.year, g.key)); return np.array([(y, k) in s for y, k in zip(d.year, d.key)]).astype(int)
f = pd.DataFrame({"year": d.year, "fbd_region": d.fbd_region})
f["imp_mvpa"] = flag("moderate_pa_rate")
f["imp_fac"] = flag("public_sports_facilities_per_100k")
f["imp_vil"] = flag("village_sports_facilities_per_100k")
f["imp_tert"] = (~d.year.isin([2015, 2020])).astype(int)
f["imp_bl"] = 0
print("[prep] imputed cells", {k: int(f[k].sum()) for k in ["imp_mvpa", "imp_fac", "imp_vil", "imp_tert", "imp_bl"]})
assert f.imp_mvpa.sum() == 229 and f.imp_fac.sum() == 239
km2 = d.area_m2 / 1e6
for v, nm in [("public_sports_facilities_per_100k", "pub"), ("village_sports_facilities_per_100k", "vil")]:
    dens = d[v + "_raw"] * d.population / 1e5 / km2
    c = dens[dens > 0].min() / 2
    x = np.log(dens + c); f[f"{nm}_km2_log"] = x; f[f"{nm}_km2_z"] = (x - x.mean()) / x.std(ddof=0); f[f"{nm}_km2_raw"] = dens
    print(f"[prep] {nm} per km2: median {dens.median():.3f}, offset c={c:.4g}, zeros {(dens==0).sum()}, corr(z per100k, z per km2)={np.corrcoef(d[v], f[f'{nm}_km2_z'])[0,1]:.3f}")
TERR = {"세종시세종시": "2012-07 세종특별자치시 출범(연기군 폐지, 공주시·청원군 일부 편입); 2011 행은 연기군",
        "충청북도청주시": "2014-07 청주시·청원군 통합",
        "충청남도공주시": "2012-07 일부 읍면(장기면 등) 세종시로 편입"}
NAME = {"대구시군위군": "2023-07 경상북도→대구광역시 편입(경계 불변, 시도 변경)",
        "인천시남구": "2018-07 남구→미추홀구 명칭 변경(경계 불변)",
        "충청남도당진시": "2012-01 당진군→당진시 승격(경계 불변)",
        "경기도여주시": "2013-09 여주군→여주시 승격(경계 불변)",
        "경기도부천시": "2016-07 일반구(원미·소사·오정) 폐지(시 경계 불변)"}
f["bnd_terr"] = f.fbd_region.isin(TERR).astype(int); f["bnd_any"] = f.fbd_region.isin(set(TERR) | set(NAME)).astype(int)
assert f.groupby("fbd_region").bnd_terr.first().sum() == 3 and f.groupby("fbd_region").bnd_any.first().sum() == 8
bc = d.groupby("fbd_region").boundary_change_flag.max(); print("[prep] panel boundary_change_flag regions:", sorted(bc[bc == 1].index.tolist()))
VARS = {"imp_mvpa": "moderate_pa_rate", "imp_fac": "public_sports_facilities_per_100k", "imp_vil": "village_sports_facilities_per_100k",
        "imp_tert": "tertiary_educated_share", "imp_bl": "basic_livelihood_rate"}
for flg, v in VARS.items():
    obs = d[v].where(f[flg] == 0)
    sd = obs.groupby(d.fbd_region).transform(lambda s: s.std(ddof=1))
    sd = sd.fillna(obs.groupby(d.fbd_region).transform("std").median())
    f[f"sd_{v}"] = sd
    print(f"[prep] MI noise SD {v}: median {sd.median():.3f} (z), cells imputed {int(f[flg].sum())}")
f.to_csv(os.path.join(OUT, "extra_inputs.csv"), index=False, encoding="utf-8-sig")
json.dump({"boundary_territorial": TERR, "boundary_name_status": NAME, "bl_regions": [], "bl_note": "수급률 대체 칸 없음(공식 관측 복원값)",
           "panel_boundary_change_flag_regions": sorted(bc[bc == 1].index.tolist())},
          open(os.path.join(OUT, "D_boundary_list.json"), "w"), ensure_ascii=False, indent=1)
print("[prep] saved", f.shape)
