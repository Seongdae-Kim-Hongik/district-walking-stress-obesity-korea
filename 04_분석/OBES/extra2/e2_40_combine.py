#!/usr/bin/env python3
# Combination of all sensitivity specifications: long table, retention of the six robust associations, supplementary table drafts, summary figures. Environment: E2_RUN (run root with 02_패널/ and 04_분석/OBES/ model inputs), E2_OUT (results folder, default 04_분석/OBES/results_extra2). Also E2_MAINRES (default results_std3) and E2_EXTRA (default results_extra).
import os, glob, json, math, unicodedata
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
NFC = lambda s: unicodedata.normalize("NFC", str(s))
OBD = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.environ.get("E2_OUT", os.path.join(OBD, "results_extra2")); O = OBD; MAINRES = os.environ.get("E2_MAINRES", os.path.join(O, "results_std3")); EXV = os.environ.get("E2_EXTRA", os.path.join(O, "results_extra"))
T = os.path.join(OUT, "tables"); FIG = os.path.join(OUT, "figures"); RUN = os.environ["E2_RUN"]
mpl.rcParams.update({"font.family": ["Arial", "Helvetica Neue"], "font.size": 8, "pdf.fonttype": 42, "axes.linewidth": 0.6})
CFG = json.load(open(os.path.join(RUN, "04_분석/OBES/config_OBES_std3.json"))); HYP = CFG["sets_after_vif"][CFG["rules"]["main_set"]]; FORCED = CFG["forced"]; MAINV = HYP + [f for f in FORCED if f not in HYP]
META = pd.read_csv(os.path.join(RUN, "04_분석/OBES/panel_META_OBES_std3.csv"), encoding="utf-8-sig").set_index("variable"); LAB = META.label_en.to_dict()
LAB.update({"eloc__dt_1yl20421e__1___1": "Green zone area per capita (log)", "ext__restaurant_per10k": "General restaurants per 10,000 (log)", "ext__instant_food_per10k": "Instant-food outlets per 10,000 (log)"})
SH = {"walking_practice_rate": "Walking", "tertiary_educated_share": "Tertiary education", "perceived_stress_rate": "Perceived stress", "high_risk_drinking_rate_drinkers": "High-risk drinking", "ext__restaurant_per10k": "Restaurant density", "current_smoking_rate": "Current smoking"}
ROB6 = list(SH)
M3 = pd.read_csv(os.path.join(MAINRES, "S3_main_OR.csv")).set_index("cov"); sgn = np.sign(np.log(M3.OR))
f3 = lambda r: f"{r.OR:.3f} ({r.lo:.3f}–{r.hi:.3f}){'*' if bool(r.cred) else ''}"
L = []
def add(df, group, name_col=None, name=None, prefix=""):
    if df is None or not len(df) or "cov" not in df: return
    df = df.copy(); df["analysis"] = (prefix + df[name_col].astype(str)) if name_col else name; df["group"] = group
    L.append(df[["group", "analysis", "cov", "OR", "lo", "hi", "cred"]])
rd = lambda f: pd.read_csv(os.path.join(MAINRES, f)) if os.path.exists(os.path.join(MAINRES, f)) else None
add(M3.reset_index(), "Principal", name="S3 main model")
add(rd("S8_KH.csv"), "Space-time interaction", "model"); add(rd("S13_graph_prior.csv"), "Graph and priors", "variant"); add(rd("S25_bym2_prior3.csv"), "Graph and priors", "phi_prior", prefix="phi prior ")
add(rd("S20_knn_graph.csv"), "Graph and priors", "variant"); add(rd("S16_RSR.csv"), "Spatial structure", "variant"); add(rd("S24_spatialplus.csv"), "Spatial structure", "variant"); add(rd("S26_spec_OR.csv").query("variant in ['iid_area','icar','no_spatial','parsimonious_S5core']") if rd("S26_spec_OR.csv") is not None else None, "Spatial structure", "variant")
add(rd("S12_likelihood.csv"), "Likelihood and time", "likelihood"); add(rd("S29_standard_extras.csv"), "Likelihood and time", "variant")
s10 = rd("S10_lag.csv"); add(s10, "Timing and measurement", "lag", prefix="survey-wave lag ") ; add(rd("S22_period.csv"), "Timing and measurement", "variant"); add(rd("S19_step2018.csv"), "Timing and measurement (post hoc)", "variant")
add(rd("S18_influence.csv"), "Influence and extremes", "variant"); add(rd("S21_leave_one_year.csv"), "Leave one year out", "left_out_year", prefix="drop ")
add(rd("S11_leave_one_sido.csv"), "Leave one province out", "left_out", prefix="drop ")
add(rd("S27_kchs_covariates.csv"), "Confounding and covariates", "variant"); add(rd("S28_forced_drop.csv"), "Confounding and covariates", "variant"); add(rd("S31_transform_single.csv"), "Confounding and covariates", "variant")
for f in sorted(glob.glob(os.path.join(EXV, "fits", "fit_*_OR.csv"))):
    d = pd.read_csv(f); nm = os.path.basename(f)[4:-7]
    if nm.startswith("F_MI") or nm == "orig_rerun" or "cov" not in d: continue
    if "component" in d: continue
    add(d, "Extra (A–E)", name="extra " + nm)
mi = pd.read_csv(os.path.join(EXV, "F_MI_rubin.csv")); add(mi, "Timing and measurement", name="imputation noise, Rubin (20)")
GRP = {"access": "Confounding and access (new)", "care": "Confounding and access (new)", "imputed": "Imputed covariates (new)", "instant": "Timing and measurement (new)", "s2": "Confounding and covariates (new)", "residspec": "Spec. for residual Moran (refit)",
       "callag": "Timing and measurement (new)", "G88": "Likelihood, smoothing, measurement (new)", "gatethr": "Data gate thresholds (new)", "W85": "Calendar windows (refit with maps)"}
for f in sorted(glob.glob(os.path.join(OUT, "fits", "fit_*_OR.csv"))):
    d = pd.read_csv(f); nm = os.path.basename(f)[4:-7]
    if "cov" not in d or nm.startswith(("residspec_S9", "nocov")): continue
    add(d, GRP.get(nm.split("_")[0], "new"), name=nm)
LONG = pd.concat(L, ignore_index=True); LONG["cred"] = LONG.cred.astype(str).str.upper().eq("TRUE")
LONG.to_csv(os.path.join(OUT, "E2_sensitivity_long_all.csv"), index=False)
R6 = LONG[LONG["cov"].isin(ROB6)].copy(); R6["kept"] = R6.cred & (np.sign(np.log(R6.OR)) == R6["cov"].map(sgn))
ret = R6.groupby("cov").agg(n_specs=("analysis", "nunique"), kept=("kept", "sum"), OR_min=("OR", "min"), OR_max=("OR", "max")).reindex(ROB6); ret["kept_frac"] = (ret.kept / ret.n_specs).round(3)
ret["main"] = [f3(M3.loc[v].rename("x").to_frame().T.assign(cov=v).iloc[0]) for v in ROB6]; ret.to_csv(os.path.join(OUT, "E2_retention_robust6.csv"))
lost = R6[~R6.kept].sort_values(["cov", "group"]); lost.to_csv(os.path.join(OUT, "E2_robust6_not_kept_list.csv"), index=False)
NEWR = R6[R6.group.str.contains(r"\(new\)|refit")]; byan = NEWR.groupby("analysis").apply(lambda d: ", ".join(SH[c] for c in ROB6 if c in set(d["cov"]) and not d[d["cov"] == c].kept.iloc[0])).rename("lost")
def wide(names, label=None, covs=None):
    rows = []
    for nm in names:
        d = LONG[LONG.analysis == nm].set_index("cov")
        if not len(d): continue
        r = {"analysis": label.get(nm, nm) if label else nm}
        for v in (covs or ROB6): r[SH.get(v, LAB.get(v, v))] = f3(d.loc[v]) if v in d.index else "—"
        rows.append(r)
    return pd.DataFrame(rows)
A = ["S3 main model", "access_plus_road", "access_urban_q5", "access_urban_rw2", "access_urban_q5_plus_road", "access_plus_outside_care", "care_plus_care_institutions", "care_plus_hospitals", "care_plus_outside_and_hospitals"]
addv = {"access_plus_road": "road_pavement_rate", "access_urban_q5_plus_road": "road_pavement_rate", "access_plus_outside_care": "outside_care_share", "care_plus_care_institutions": "care_institutions_per_100k", "care_plus_hospitals": "hospitals_per_100k", "care_plus_outside_and_hospitals": "outside_care_share"}
tb = wide(A); tb["added covariate OR"] = [f3(LONG[(LONG.analysis == a) & (LONG["cov"] == addv[a])].iloc[0]) if a in addv else "—" for a in A]
tb.loc[tb.analysis == "care_plus_outside_and_hospitals", "added covariate OR"] = "outside " + f3(LONG[(LONG.analysis == "care_plus_outside_and_hospitals") & (LONG["cov"] == "outside_care_share")].iloc[0]) + "; hospitals " + f3(LONG[(LONG.analysis == "care_plus_outside_and_hospitals") & (LONG["cov"] == "hospitals_per_100k")].iloc[0])
tb.to_csv(os.path.join(T, "SI_access_care_access_urbanicity.csv"), index=False)
wide(["S3 main model", "imputed_drop_all9_imputed_vars", "imputed_drop_all11_incl_census_interp", "imputed_cc_no_imputed_cell"]).to_csv(os.path.join(T, "SI_imputed_imputed_all_at_once.csv"), index=False)
wide(["S3 main model", "survey-wave lag 1", "survey-wave lag 2", "callag_cal_lag1_A", "callag_nolag_matched_lag1_A", "callag_cal_lag1_B", "callag_cal_lag2_A", "callag_nolag_matched_lag2_A", "callag_cal_lag2_B"]).to_csv(os.path.join(T, "SI_callag_calendar_lag.csv"), index=False)
wide(["S3 main model", "gatethr_a205k25", "gatethr_a205k20"]).to_csv(os.path.join(T, "SI_gatethr_gate_threshold_fits.csv"), index=False)
wide(["S3 main model", "G88_poisson_lognormal", "G88_smooth_sexratio_pop80", "G88_kchs_cov_ma3", "G88_extreme_5IQR_missing", "s2_no_pop80_forced", "residspec_S13_queen_s577", "instant_years_2016_2024"]).to_csv(os.path.join(T, "SI_extreme_stdsets_standard_suite.csv"), index=False)
wide(["S3 main model"] + [f"W85_{w}" for w in ["W1_2011-15", "W2_2016-18", "W3_2020-24"]], covs=MAINV).to_csv(os.path.join(T, "SI_wingrid_window_ORs_all18.csv"), index=False)
S8 = pd.read_csv(os.path.join(MAINRES, "S8_KH.csv")); S9 = pd.read_csv(os.path.join(MAINRES, "S9_between_within.csv")); rows = []
for v in HYP + [f for f in FORCED]:
    r = {"Characteristic": LAB.get(v, v), "Main (type I)": f3(M3.loc[v].rename("x").to_frame().T.iloc[0])}
    for mdl, nm in [("KH_II", "Type II"), ("KH_III", "Type III"), ("KH_IV", "Type IV")]:
        q = S8[(S8.model == mdl) & (S8["cov"] == v)]; r[nm] = f3(q.iloc[0]) if len(q) else "—"
    for cmp_, nm in [("between", "Between-district"), ("within", "Within-district")]:
        q = S9[(S9.component == cmp_) & (S9["cov"] == v)]; r[nm] = f3(q.iloc[0]) if len(q) else "(forced)"
    rows.append(r)
pd.DataFrame(rows).to_csv(os.path.join(T, "SI_intervals_Table2_with_CrI.csv"), index=False)
p0 = 0.290; o0 = p0 / (1 - p0); pp = lambda o: 100 * (o0 * o / (1 + o0 * o) - p0)
PPD = pd.DataFrame([dict(cov=v, label=SH[v], OR=M3.loc[v, "OR"], pp_diff=round(pp(M3.loc[v, "OR"]), 2), pp_lo=round(pp(M3.loc[v, "lo"]), 2), pp_hi=round(pp(M3.loc[v, "hi"]), 2)) for v in ROB6]); PPD.to_csv(os.path.join(T, "intervals_pp_difference_per_SD.csv"), index=False)
S5 = pd.read_csv(os.path.join(MAINRES, "S5_covid.csv")); S5["label"] = S5["cov"].map(LAB); S5["COVID interaction OR (95% CrI)"] = S5.apply(f3, axis=1)
S5[["label", "COVID interaction OR (95% CrI)", "P_gt1"]].to_csv(os.path.join(T, "SI_covidint_slopes_covid_interaction_13.csv"), index=False)
S4i = pd.read_csv(os.path.join(MAINRES, "S4_interaction.csv")); S4i["label"] = S4i["cov"].map(LAB); S4i["v"] = S4i.apply(f3, axis=1)
S4i.pivot(index="label", columns="window", values="v").reindex([LAB[v] for v in HYP]).to_csv(os.path.join(T, "SI_covidint_window_interaction_26.csv"))
S17 = pd.concat([pd.read_csv(os.path.join(MAINRES, "S17_region_slopes.csv")), pd.read_csv(os.path.join(OUT, "slopes_S17_region_slopes_remaining3.csv"))], ignore_index=True)
S17["v"] = S17.apply(f3, axis=1); S17["label"] = S17.covariate.map(LAB)
s17w = S17.pivot(index="label", columns="region", values="v").reindex([LAB[v] for v in HYP]); s17w["dWAIC vs common"] = S17.groupby("label").dWAIC_vs_common.first().reindex(s17w.index)
s17w.to_csv(os.path.join(T, "SI_slopes_region_slopes_all13.csv"))
S26 = pd.read_csv(os.path.join(MAINRES, "S26_resid_moran_by_spec.csv"))[["variant", "WAIC", "phi", "resid_moran_I", "resid_moran_p"]]
NS = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(OUT, "fits", "fit_residspec_*_stat.csv")))], ignore_index=True).rename(columns={"fit": "variant"})[["variant", "WAIC", "phi", "resid_moran_I", "resid_moran_p"]]
pd.concat([S26, NS], ignore_index=True).to_csv(os.path.join(T, "SI_residspec_resid_moran_by_spec.csv"), index=False)
nh = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(OUT, "nullhot_null_hotspots_J*.csv")))], ignore_index=True); nh.to_csv(os.path.join(OUT, "nullhot_null_hotspots_all.csv"), index=False)
F3f = pd.read_csv(os.path.join(MAINRES, "S3_main_fit.csv")).iloc[0]
pd.DataFrame([dict(n_sims=len(nh), null_hot_mean=nh.n_hot.mean(), null_hot_max=nh.n_hot.max(), null_cold_mean=nh.n_cold.mean(), null_cold_max=nh.n_cold.max(), null_phi_mean=round(nh.phi.mean(), 3), null_sd_logOR_mean=round(nh.sd_logOR.mean(), 4), null_sd_logOR_max=round(nh.sd_logOR.max(), 4),
               obs_hot=int(F3f.n_hotspot), obs_cold=int(F3f.n_coldspot), obs_phi=F3f.phi, p_hot=(1 + (nh.n_hot >= F3f.n_hotspot).sum()) / (len(nh) + 1), p_cold=(1 + (nh.n_cold >= F3f.n_coldspot).sum()) / (len(nh) + 1))]).to_csv(os.path.join(T, "nullhot_null_hotspot_summary.csv"), index=False)
NC = pd.read_csv(os.path.join(MAINRES, "S3_negcontrol.csv")); NC["label"] = NC.label_en; MM = M3.reset_index(); MM["label"] = MM["cov"].map(LAB)
fig, ax = plt.subplots(figsize=(5.2, 6.0)); allr = pd.concat([MM.assign(kind=np.where(MM.role == "forced", "Control", "Characteristic")), NC.assign(kind="Negative control")], ignore_index=True)
allr = allr.iloc[::-1].reset_index(drop=True)
for i, r in allr.iterrows():
    c = {"Characteristic": "#1f4e79", "Control": "#7f7f7f", "Negative control": "#b2182b"}[r.kind]; ax.plot([r.lo, r.hi], [i, i], color=c, lw=1.1); ax.plot(r.OR, i, "o", color=c, mfc=c if r.cred else "white", ms=4)
ax.axvline(1, color="black", lw=0.6); ax.set_yticks(range(len(allr))); ax.set_yticklabels(allr.label, fontsize=6.8); ax.set_xlabel("Odds ratio per SD (95% CrI); filled = interval excludes 1")
ax.legend(handles=[mpl.lines.Line2D([0], [0], color=c, marker="o", lw=1, label=k) for k, c in [("Characteristic", "#1f4e79"), ("Control", "#7f7f7f"), ("Negative control (added one at a time)", "#b2182b")]], fontsize=6.5, frameon=False, loc="lower right")
plt.tight_layout(); fig.savefig(os.path.join(FIG, "FigS_negctrl_negative_controls.png"), dpi=300, facecolor="white"); plt.close(fig)
EV = pd.read_csv(os.path.join(EXV, "G1_evalue_BHq.csv")); EV = EV.sort_values("E_point_common_sqrt")
fig, ax = plt.subplots(figsize=(5.2, 4.2)); y = np.arange(len(EV))
ax.barh(y, EV.E_point_common_sqrt - 1, left=1, color=np.where(EV.cred, "#1f4e79", "#bdbdbd"), height=0.55, label="E-value, point estimate"); ax.plot(EV.E_ci_common_sqrt, y, "k|", ms=9, label="E-value, interval limit closer to 1")
ax.set_yticks(y); ax.set_yticklabels(EV.label_en, fontsize=6.8); ax.set_xlabel("E-value (risk-ratio scale, common-outcome square-root conversion of the OR)"); ax.legend(fontsize=6.5, frameon=False, loc="lower right")
plt.tight_layout(); fig.savefig(os.path.join(FIG, "FigS_evalue_evalues.png"), dpi=300, facecolor="white"); plt.close(fig)
mi2 = mi.merge(M3.reset_index()[["cov", "OR", "lo", "hi"]], on="cov", suffixes=("_mi", "_main")); mi2["label"] = mi2["cov"].map(LAB); mi2 = mi2.iloc[::-1].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(5.2, 5.4))
for i, r in mi2.iterrows():
    ax.plot([r.lo_main, r.hi_main], [i + .15, i + .15], color="#7f7f7f", lw=1.1); ax.plot(r.OR_main, i + .15, "o", color="#7f7f7f", ms=3.5)
    ax.plot([r.lo_mi, r.hi_mi], [i - .15, i - .15], color="#1f4e79", lw=1.1); ax.plot(r.OR_mi, i - .15, "s", color="#1f4e79", ms=3.5)
ax.axvline(1, color="black", lw=0.6); ax.set_yticks(range(len(mi2))); ax.set_yticklabels(mi2.label, fontsize=6.8); ax.set_xlabel("Odds ratio per SD (95% CrI)")
ax.legend(handles=[mpl.lines.Line2D([0], [0], color="#7f7f7f", marker="o", label="Main model (single imputation)"), mpl.lines.Line2D([0], [0], color="#1f4e79", marker="s", label="20 noise-added imputations, Rubin's rules")], fontsize=6.5, frameon=False, loc="lower right")
plt.tight_layout(); fig.savefig(os.path.join(FIG, "FigS_impunc_imputation_uncertainty.png"), dpi=300, facecolor="white"); plt.close(fig)
AT = R6[~R6.group.isin(["Leave one year out", "Leave one province out"])].copy(); AT["key"] = AT.group + " | " + AT.analysis
order = list(dict.fromkeys(AT.sort_values(["group"], key=lambda s: s.map({g: i for i, g in enumerate(dict.fromkeys(LONG.group))})).key))
Mx = AT.pivot_table(index="key", columns="cov", values="OR", aggfunc="first").reindex(order)[ROB6]; Kp = AT.pivot_table(index="key", columns="cov", values="kept", aggfunc="first").reindex(order)[ROB6]; Cr = AT.pivot_table(index="key", columns="cov", values="cred", aggfunc="first").reindex(order)[ROB6]
rel = np.log(Mx.astype(float)) / np.log(M3.loc[ROB6, "OR"].values)
fig, ax = plt.subplots(figsize=(6.6, 0.16 * len(order) + 1.2)); im = ax.imshow(rel.values.astype(float), cmap="RdYlGn", vmin=-0.5, vmax=2.5, aspect="auto")
for i in range(len(order)):
    for j in range(6):
        if pd.isna(Mx.iat[i, j]): continue
        k, c = Kp.iat[i, j], Cr.iat[i, j]; ax.text(j, i, "●" if k else ("×" if c else "○"), ha="center", va="center", fontsize=5.5)
ax.set_xticks(range(6)); ax.set_xticklabels([SH[v] for v in ROB6], rotation=30, ha="right", fontsize=7); ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=4.6)
fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01, label="log OR relative to main model (1 = same size)"); ax.set_title("● credible, same direction   ○ interval includes 1   × credible, opposite direction", fontsize=6.5)
plt.tight_layout(); fig.savefig(os.path.join(FIG, "FigS_atlas_sensitivity_atlas.png"), dpi=300, facecolor="white"); plt.close(fig)
lines = ["# OBES sensitivity summary", "", "Main conclusion = six robust associations (walking -, tertiary education -, stress +, high-risk drinking +, restaurant density +, smoking +; 95% CrI excludes 1 in the main model). An analysis keeps the conclusion when all six stay credible in the same direction.", "",
         f"All sensitivity specifications (including leave-one-province and leave-one-year): {LONG.analysis.nunique()}; kept (credible, same direction) / specifications containing the covariate:", ""]
lines += ["| Covariate | Main model | Kept / specifications | OR range |", "|---|---|---|---|"] + [f"| {SH[v]} | {ret.loc[v, 'main']} | {int(ret.loc[v, 'kept'])} / {int(ret.loc[v, 'n_specs'])} | {ret.loc[v, 'OR_min']:.3f}–{ret.loc[v, 'OR_max']:.3f} |" for v in ROB6]
lines += ["", "## By round-2 analysis", "", "| Analysis | Conclusion kept? | Lost among the six |", "|---|---|---|"]
for an in sorted(NEWR.analysis.unique()):
    lo_ = byan.get(an, ""); lines.append(f"| {an} | {'yes' if not lo_ else 'no (partly)'} | {lo_ or '—'} |")
open(os.path.join(OUT, "00_sensitivity_summary.md"), "w").write("\n".join(lines) + "\n")
print(ret.to_string()); print(byan.to_string())
