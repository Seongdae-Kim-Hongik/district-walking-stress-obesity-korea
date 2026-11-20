# Reads result tables and computes every reported number.
import os, json, unicodedata
import pandas as pd, numpy as np
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); AN = os.path.normpath(os.path.join(HERE, "..", "..", "..")); ROOT = os.path.dirname(AN)
R = os.path.join(AN, "OBES", "results_std3")
rd = lambda f: pd.read_csv(os.path.join(R, f))
N = {}
FO = pd.read_csv(os.path.join(ROOT, "02_패널", "OBES_fh_outcomes.csv"), encoding="utf-8-sig"); FO = FO[FO.year != 2019]
assert len(FO) == 2977 and FO.SR.notna().all()
N["n_cells"] = len(FO); N["n_resp_total"] = int(FO.N.sum()); N["n_resp_median"] = int(FO.N.median()); N["n_resp_min"] = int(FO.N.min())
pw = pd.Series({y: np.average(g.SR, weights=g.population) for y, g in FO.groupby("year")})
N["sr_2011"] = round(pw[2011], 1); N["sr_2024"] = round(pw[2024], 1); N["sr_2017"] = round(pw[2017], 1); N["sr_2018"] = round(pw[2018], 1)
p = FO.groupby("fbd_region").SR.mean(); N["p10"] = round(p.quantile(.1), 1); N["p90"] = round(p.quantile(.9), 1); N["p9010"] = round(p.quantile(.9) / p.quantile(.1), 2)
N["se_median"] = round(FO.SR_SE.median(), 1)
s0 = rd("S0_screen.csv").iloc[0]; N["moran"] = round(s0.moran_I_pooled_logitSR, 2); N["moran_p"] = s0.moran_p; N["moran_ymin"] = round(s0.moran_I_year_min, 2); N["moran_ymax"] = round(s0.moran_I_year_max, 2)
N["persist"] = round(s0.persistence_rank_rho, 2); N["phi_forced"] = round(s0.phi_forced_only, 2)
N["kendallW"] = round(rd("S6_kendallW.csv").kendall_W[0], 2)
f3 = rd("S3_main_fit.csv").iloc[0]; N["phi_main"] = round(f3.phi, 2); N["n_hot"] = int(f3.n_hotspot); N["n_cold"] = int(f3.n_coldspot); N["resid_moran"] = round(f3.resid_moran_I, 2); N["resid_moran_p"] = round(f3.resid_moran_p, 2)
N["largest_hot"] = int(f3.largest_hot); N["largest_cold"] = int(f3.largest_cold)
sp = rd("S3_spatial_OR.csv"); hot = sp[sp["class"] == "hotspot"].sido.value_counts(); cold = sp[sp["class"] == "coldspot"].sido.value_counts()
N["hot_sido"] = hot.to_dict(); N["cold_sido"] = cold.to_dict()
N["or_min"] = sp.sort_values("OR").head(3)[["fbd_region", "OR"]].values.tolist(); N["or_max"] = sp.sort_values("OR").tail(3)[["fbd_region", "OR"]].values.tolist()
vf = rd("S3_vs_forced_only.csv").iloc[0]; N["sp_rho_vs_forced"] = round(vf.spearman_spatial_OR, 2); N["sd_main"] = round(vf.sd_logOR_main, 3); N["sd_forced"] = round(vf.sd_logOR_forced, 3)
sf = rd("S0_forced_only_spatialOR.csv"); N["seoul_cold_forced"] = int(((sf["class"] == "coldspot") & (sf.sido == "서울시")).sum()); N["seoul_cold_main"] = int(((sp["class"] == "coldspot") & (sp.sido == "서울시")).sum())
vp = rd("S3_varpart.csv").set_index("component").share; N["vp"] = {k: round(100 * v, 0) for k, v in vp.items()}
nt = vp.drop("temporal"); nt = nt / nt.sum(); N["vp_notime"] = {k: round(100 * v, 0) for k, v in nt.items()}
m = rd("S3_main_OR.csv").set_index("cov"); N["main"] = {k: dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for k, r in m.iterrows()}
mv = rd("S7_multiverse_summary.csv").set_index("cov"); N["mv"] = {k: dict(cred=r.cred_frac, sign=r.sign_agree, robust=bool(r.robust)) for k, r in mv.iterrows()}
kh = rd("S8_KH.csv"); N["kh"] = {f"{r.model}|{r['cov']}": dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for _, r in kh.iterrows()}
ROB = [k for k in mv.index if mv.loc[k, "robust"] and kh[(kh["cov"] == k) & kh.model.isin(["KH_II", "KH_IV"])].cred.all()]
assert set(ROB) == {"walking_practice_rate", "ext__restaurant_per10k", "high_risk_drinking_rate_drinkers", "current_smoking_rate", "perceived_stress_rate", "tertiary_educated_share"}, ROB
N["robust"] = ROB
bw = rd("S9_between_within.csv"); N["bw"] = {f"{r['cov']}|{r.component}": dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for _, r in bw.iterrows()}
nc = rd("S3_negcontrol.csv"); N["neg"] = {r["cov"]: dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for _, r in nc.iterrows()}; assert not nc.cred.any()
s2 = rd("S2_setcompare.csv").set_index("set").WAIC; N["waic"] = s2.round(1).to_dict()
lo1 = rd("S11_leave_one_sido.csv"); N["loo"] = lo1.groupby("cov").cred.sum().astype(int).to_dict(); N["loo_n"] = int(lo1.left_out.nunique())
s14 = rd("S14_measured_summary.csv").iloc[0]; N["meas"] = {k: (float(v) if isinstance(v, (int, float, np.floating, np.integer)) else v) for k, v in s14.items()}
mo = rd("S14_measured_OR.csv").set_index("cov"); N["meas_OR"] = {k: dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for k, r in mo.iterrows()}
rp = rd("S44_reverse_posthoc.csv"); N["rev"] = {f"{r['cov']}|{r.test}": dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for _, r in rp.iterrows()}
ye = rd("S3_year_effect.csv").set_index("year"); N["rw1_2017"] = round(ye.rw1[2017], 3); N["rw1_2018"] = round(ye.rw1[2018], 3)
w = rd("S4_window_OR.csv"); N["win"] = {f"{r.window}|{r['cov']}": dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for _, r in w.iterrows()}
l = rd("S10_lag.csv"); N["lag"] = {f"{r.lag}|{r['cov']}": dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for _, r in l.iterrows()}
li = rd("S12_likelihood.csv"); N["lik"] = {f"{r.likelihood}|{r['cov']}": dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for _, r in li.iterrows()}
rs = rd("S16_RSR.csv"); N["rsr"] = {r["cov"]: dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for _, r in rs.iterrows()}
rg = rd("S17_region_slopes.csv"); N["reg"] = {f"{r.covariate}|{r.region}": dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred), dW=r.dWAIC_vs_common) for _, r in rg.iterrows()}
inf = rd("S18_influence.csv"); N["infl_n"] = int(inf.variant.iloc[0].split("_")[1])
cv = rd("S5_covid.csv"); N["covid"] = {r["cov"]: dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for _, r in cv.iterrows()}
s19 = rd("S19_step2018.csv"); N["s19"] = {f"{r.variant}|{r['cov']}": dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for _, r in s19.iterrows()}
mvr = rd("S7_multiverse_results.csv"); N["mv_n"] = mvr.groupby("cov").spec.nunique().astype(int).to_dict()
so = rd("S14_selfreport2018_OR.csv"); N["self18_OR"] = {r["cov"]: dict(OR=r.OR, lo=r.lo, hi=r.hi, cred=bool(r.cred)) for _, r in so.iterrows()}
META = pd.read_csv(os.path.join(AN, "OBES", "panel_META_OBES_std3.csv"), encoding="utf-8-sig").set_index("variable")
N["filled"] = {v: int(META.loc[v, "n_imputed_interp"] + META.loc[v, "n_imputed_sgg_mean"] + META.loc[v, "n_imputed_year_mean"]) for v in META.index}
N["sr_mean_all"] = round(float(FO.SR.mean()), 1)
json.dump(N, open(os.path.join(HERE, "numbers_obes.json"), "w"), ensure_ascii=False, indent=1, default=float)
print("robust", ROB); print({k: N[k] for k in ["n_resp_total", "n_resp_median", "sr_2011", "sr_2024", "p10", "p90", "moran", "persist", "kendallW", "phi_main", "n_hot", "n_cold", "vp", "vp_notime", "infl_n"]})
