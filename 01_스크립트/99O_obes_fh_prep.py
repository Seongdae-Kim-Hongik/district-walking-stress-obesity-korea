#!/bin/bash
# Survey outcomes (logit, delta-method sampling variance), effective sample sizes, and measured obesity.
import os, unicodedata
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, ".."))
M = pd.read_parquet(os.path.join(ROOT, "02_패널", "MASTER_시군구_연도_wide.parquet")); M.columns = [NFC(c) for c in M.columns]
M = M[(M.year >= 2011) & (M.year <= 2024)].copy(); M["fbd_region"] = M.fbd_region.map(NFC)
tb, mb = "kchs__DT_H_OBE_OBE__", "kchs__DT_117075_MEA_OBE__"
d = pd.DataFrame({"year": M.year.values, "fbd_region": M.fbd_region.values, "population": M.population.astype(float).values})
for it in ["SR", "SR_SE", "N", "CR", "CR_SE"]: d[it] = M[tb + it].astype(float).values
d["meas_SR"] = M[mb + "SR"].astype(float).values; d["meas_SE"] = M[mb + "SR_SE"].astype(float).values
assert len(d) == 229 * 14 and d.loc[d.year == 2019, "SR"].isna().all()
obs = d.year != 2019
assert d.loc[obs, ["SR", "SR_SE", "N", "CR", "CR_SE"]].notna().all().all() and (d.loc[obs, "SR"] > 0).all() and (d.loc[obs, "SR"] < 100).all()
n0 = int((d.SR_SE <= 0).sum()); d["SE_used"] = d.SR_SE.clip(lower=0.05); d["CRSE_used"] = d.CR_SE.clip(lower=0.05)
p = d.SR / 100; d["y_logit"] = np.log(p / (1 - p)); d["v_logit"] = (d.SE_used / 100) ** 2 / (p * (1 - p)) ** 2
pc = d.CR / 100; d["y_logit_cr"] = np.log(pc / (1 - pc)); d["v_logit_cr"] = (d.CRSE_used / 100) ** 2 / (pc * (1 - pc)) ** 2
d["n_eff"] = np.maximum(np.round(p * (1 - p) / (d.SE_used / 100) ** 2), 10); d["y_eff"] = np.round(p * d.n_eff)
pm = d.meas_SR / 100; d["y_logit_meas"] = np.log(pm / (1 - pm)); d["v_logit_meas"] = (d.meas_SE.clip(lower=0.05) / 100) ** 2 / (pm * (1 - pm)) ** 2
d["report_gap"] = d.y_logit_meas - d.y_logit
print(f"[99O] rows {len(d)} · 관측 {int(obs.sum())} · SE≤0 칸 {n0} · SR {d.SR.min():.1f}-{d.SR.max():.1f} · n_eff 중앙 {d.n_eff.median():.0f}")
print(f"[99O] 측정 2018 {int(d.loc[d.year==2018,'meas_SR'].notna().sum())} · 2019 {int(d.loc[d.year==2019,'meas_SR'].notna().sum())} 시군구, 2018 보고 격차 평균 {d.loc[d.year==2018,'report_gap'].mean():.3f}")
fp = os.path.join(ROOT, "02_패널", "OBES_fh_outcomes.csv"); d.to_csv(fp, index=False, encoding="utf-8-sig"); print("[99O] 저장", fp, d.shape)
