#!/usr/bin/env python3
# Supplementary table drafts (markdown) for the additional analyses. Environment: OBES_RES, OBES_EXTRA.
import os, numpy as np, pandas as pd
OB = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")); OUT = os.environ.get("OBES_EXTRA", os.path.join(OB, "results_extra"))
S3 = pd.read_csv(os.path.join(os.environ.get("OBES_RES", os.path.join(OB, "results_std3")), "S3_main_OR.csv")); HYP = list(S3[S3.role == "hypothesis"]["cov"]); LAB = dict(zip(S3["cov"], S3.label_en))
SHORT = {"walking_practice_rate": "Walking practice", "moderate_pa_rate": "Moderate-or-higher PA", "public_sports_facilities_per_100k": "Public sports facilities", "village_sports_facilities_per_100k": "Village sports facilities",
         "eloc__dt_1yl20421e__1___1": "Green zone area per capita", "cars_per_capita": "Cars per capita", "ext__restaurant_per10k": "General restaurants", "ext__instant_food_per10k": "Instant-food outlets",
         "high_risk_drinking_rate_drinkers": "High-risk drinking", "current_smoking_rate": "Current smoking", "perceived_stress_rate": "Perceived stress", "tertiary_educated_share": "Tertiary education", "basic_livelihood_rate": "Livelihood-security receipt",
         "pub_km2_z": "Public sports facilities per km² (log)", "vil_km2_z": "Village sports facilities per km² (log)"}
def md(df): return "| " + " | ".join(df.columns) + " |\n|" + "|".join(["---"] * len(df.columns)) + "|\n" + "\n".join("| " + " | ".join(str(x) for x in r) + " |" for r in df.values) + "\n"
C = pd.read_csv(os.path.join(OUT, "ABCDE_OR_compare.csv")); C["Characteristic"] = C["cov"].map(SHORT)
ST = pd.read_csv(os.path.join(OUT, "ABCDE_fitstats.csv")).set_index("fit")
def wa(f): s = ST.loc[f]; return f"n = {int(s.n_obs)}; WAIC {s.WAIC:.1f} (Δ {s.dWAIC_vs_main:+.1f})"
o = []
o.append("### Table S7. Alternative spatial structures (A)\n\nOR per 1 SD (95% CrI); * = 95% CrI excludes 1.\n\n")
t = C[C["cov"].isin(HYP)][["Characteristic", "main", "A1_icar", "A2_nospatial"]].copy(); t.columns = ["Characteristic", "Main (BYM2)", "ICAR (Besag) only", "No spatial effect"]
o.append(md(t)); o.append(f"\nFit: main {wa('orig_rerun')}; ICAR {wa('A1_icar')}; no spatial {wa('A2_nospatial')}.\n\n")
B = pd.read_csv(os.path.join(OUT, "B_twoway_vs_S9.csv")); B["Characteristic"] = B["cov"].map(SHORT)
o.append("### Table S8. Between- and within-district components: one-way (S9) versus two-way centring (B)\n\nOR per 1 SD of the pooled covariate; last column rescales the two-way within OR to 1 SD of the within component.\n\n")
o.append(md(B[["Characteristic", "S9_between", "S9_within_oneway", "B_between", "B_within_twoway", "within_twoway_SD", "B_within_per_withinSD"]].rename(columns={"S9_between": "Between (one-way)", "S9_within_oneway": "Within (one-way)", "B_between": "Between (two-way)", "B_within_twoway": "Within (two-way)", "within_twoway_SD": "SD of two-way within (pooled-SD units)", "B_within_per_withinSD": "Within OR per within-SD"})))
o.append("\n### Table S9. Sports-facility density per km² instead of per capita (C)\n\n")
t = C[["Characteristic", "main", "C1_km2", "C2_both"]].copy(); t.columns = ["Characteristic", "Main (per 100k)", "Per km² only", "Per 100k + per km²"]
o.append(md(t)); o.append(f"\nFit: per km² {wa('C1_km2')}; both {wa('C2_both')}.\n\n")
o.append("### Table S10. Exclusion of districts with boundary or status changes, 2011–2024 (D)\n\n")
t = C[C["cov"].isin(HYP)][["Characteristic", "main", "D1_bnd_terr", "D2_bnd_any"]].copy(); t.columns = ["Characteristic", "Main", "Excl. 3 territorial changes", "Excl. all 8 changes"]
o.append(md(t)); o.append(f"\nFit: {wa('D1_bnd_terr')}; {wa('D2_bnd_any')} (WAIC not comparable across different n).\n\n")
o.append("### Table S11. Sensitivity to imputed covariate cells (E)\n\n'drop' = variable removed; 'cc' = rows with an imputed value of that variable removed.\n\n")
t = C[C["cov"].isin(HYP)][["Characteristic", "main", "E1_drop_mvpa", "E1cc_mvpa", "E2_drop_fac", "E2cc_fac", "E3_drop_tert", "E3cc_tert", "E4_drop_bl"]].copy()
t.columns = ["Characteristic", "Main", "MVPA drop", "MVPA cc (−2018)", "Facilities drop", "Facilities cc (−2015 etc.)", "Tertiary drop", "Tertiary cc (2015, 2020 only)", "Livelihood drop"]
o.append(md(t)); o.append("\nn: " + "; ".join(f"{f} {int(ST.loc[f].n_obs)}" for f in ["E1cc_mvpa", "E2cc_fac", "E3cc_tert"]) + ".\n\n")
R = pd.read_csv(os.path.join(OUT, "F_MI_rubin.csv")); R["Characteristic"] = R["cov"].map(SHORT).fillna(R["cov"].map(LAB))
R["Pooled OR (95% CrI)"] = [f"{a:.3f} ({b:.3f}–{c:.3f}){'*' if d else ''}" for a, b, c, d in zip(R.OR, R.lo, R.hi, R.cred)]
o.append("### Table S12. Propagation of imputation uncertainty, M = 20 (F)\n\n")
o.append(md(R[["Characteristic", "role", "main", "Pooled OR (95% CrI)", "FMI", "n_cred_of_M"]].rename(columns={"main": "Main", "n_cred_of_M": "Credible in m of 20"})))
G = pd.read_csv(os.path.join(OUT, "G1_evalue_BHq.csv")); G["Characteristic"] = G["cov"].map(SHORT)
G["OR (95% CrI)"] = [f"{a:.3f} ({b:.3f}–{c:.3f})" for a, b, c in zip(G.OR, G.lo, G.hi)]
G["Two-sided tail"] = [("<0.0001" if p == 0 else f"{p:.4f}") for p in G.p_tail2]; G["BH q"] = [("<0.0001" if q == 0 else f"{q:.4f}") for q in G.q_BH]
G["E-value, OR≈RR"] = [f"{a:.2f} ({b:.2f})" for a, b in zip(G.E_point_rare, G.E_ci_rare)]; G["E-value, √OR (common outcome)"] = [f"{a:.2f} ({b:.2f})" for a, b in zip(G.E_point_common_sqrt, G.E_ci_common_sqrt)]
o.append("\n### Table S13. Posterior tail probabilities, Benjamini–Hochberg q and E-values for the 13 characteristics (G)\n\nE-value for the point estimate (for the CrI limit closest to 1), per 1 SD.\n\n")
o.append(md(G[["Characteristic", "OR (95% CrI)", "Two-sided tail", "BH q", "E-value, OR≈RR", "E-value, √OR (common outcome)"]]))
W = pd.read_csv(os.path.join(OUT, "G4_window_P10_P90.csv")); W["Characteristic"] = W["cov"].map(SHORT).fillna(W["label_en"])
W["v"] = [f"{a:g} / {b:g} / {c:g}" for a, b, c in zip(W.P10.round(2), W.P50.round(2), W.P90.round(2))]
P = W.pivot_table(index="Characteristic", columns="window", values="v", aggfunc="first", sort=False).reset_index()
o.append("\n### Table S14. Covariate distribution by period, P10 / P50 / P90 in natural units (G)\n\n"); o.append(md(P))
MO = pd.read_csv(os.path.join(OUT, "G5_moran_by_year.csv"))
o.append("\n### Table S15. Global Moran's I of logit standardized obesity prevalence by year (G)\n\n999 permutations, row-standardized 593-edge graph, E[I] = −0.0044.\n\n")
o.append(md(MO[["year", "moran_I", "p", "SR_popweighted"]].rename(columns={"year": "Year", "moran_I": "Moran's I", "p": "Permutation p", "SR_popweighted": "Population-weighted SR (%)"})))
open(os.path.join(OUT, "_si_tables.md"), "w").write("".join(o)); print("ok")
