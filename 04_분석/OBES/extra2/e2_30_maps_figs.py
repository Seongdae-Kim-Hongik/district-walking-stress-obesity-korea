#!/usr/bin/env python3
# Maps and figures of the round-2 analyses (window grid, window spatial OR and exceedance maps, yearly observed and smoothed maps with emerging hot spot analysis,
# two-way centred bins, observed, fitted and residual maps, hot and cold district list). All maps are joined by district code with a one-to-one and set-equality check.
# Environment: E2_RUN (run root with 02_패널/ and 04_분석/OBES/ model inputs), E2_OUT (results folder, default 04_분석/OBES/results_extra2). Also E2_SHP (boundary file), E2_MAINRES (main result folder, default 04_분석/OBES/results_std3).
import os, json, math, unicodedata
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from scipy.stats import kendalltau
NFC = lambda s: unicodedata.normalize("NFC", str(s))
OBD = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.environ.get("E2_OUT", os.path.join(OBD, "results_extra2")); O = OBD
RUN = os.environ["E2_RUN"]; SHP = os.environ.get("E2_SHP", os.environ.get("ST_SHP", "data/boundary/final.shp")); MAINRES = os.environ.get("E2_MAINRES", os.path.join(O, "results_std3")); FIG = os.path.join(OUT, "figures"); T = os.path.join(OUT, "tables")
mpl.rcParams.update({"font.family": ["Arial", "Helvetica Neue", "Apple SD Gothic Neo"], "font.size": 8, "pdf.fonttype": 42, "axes.linewidth": 0.6})
P = pd.read_csv(os.path.join(RUN, "04_분석/OBES/panel_OBES_std3.csv"), encoding="utf-8-sig"); P["fbd_region"] = P.fbd_region.map(NFC); P["sigungu_cd"] = P.sigungu_cd.astype(int)
FO = pd.read_csv(os.path.join(RUN, "02_패널/OBES_fh_outcomes.csv"), encoding="utf-8-sig"); FO["fbd_region"] = FO.fbd_region.map(NFC)
META = pd.read_csv(os.path.join(RUN, "04_분석/OBES/panel_META_OBES_std3.csv"), encoding="utf-8-sig").set_index("variable")
LAB = META.label_en.to_dict(); LAB.update({"eloc__dt_1yl20421e__1___1": "Green zone area per capita (log)", "ext__restaurant_per10k": "General restaurants per 10,000 (log)", "ext__instant_food_per10k": "Instant-food outlets per 10,000 (log)"})
CFG = json.load(open(os.path.join(RUN, "04_분석/OBES/config_OBES_std3.json"))); HYP = CFG["sets_after_vif"][CFG["rules"]["main_set"]]; FORCED = CFG["forced"]
G = gpd.read_file(SHP); G["region"] = G.region.map(NFC); G["sigungu_cd"] = G.sigungu_cd.astype(int)
key = P.drop_duplicates("fbd_region")[["fbd_region", "sigungu_cd"]]
chk = []
def cjoin(df, name, how_key="fbd_region"):
    """분석 표(df: fbd_region 또는 sigungu_cd)를 셰이프에 코드로 결합. 1:1·집합 일치 assert, 기록."""
    d = df.copy()
    if "sigungu_cd" not in d: d = d.merge(key, on="fbd_region", how="left", validate="many_to_one")
    d["sigungu_cd"] = d.sigungu_cd.astype(int)
    if "fbd_region" not in d: d = d.merge(key, on="sigungu_cd", how="left", validate="many_to_one")
    a, b = set(G.sigungu_cd), set(d.sigungu_cd)
    ok_set = a == b; ok_11 = (not d.sigungu_cd.duplicated().any()) and (not G.sigungu_cd.duplicated().any())
    nm = G.merge(d[["sigungu_cd", "fbd_region"]], on="sigungu_cd", how="left"); name_ok = bool((nm.region == nm.fbd_region).all())
    chk.append(dict(map=name, n_shape=len(a), n_analysis=len(b), set_equal=ok_set, one_to_one=ok_11, only_in_shape=";".join(map(str, sorted(a - b))), only_in_analysis=";".join(map(str, sorted(b - a))), code_join_region_name_agrees=name_ok))
    assert ok_set and ok_11, name
    return G.merge(d.drop(columns=["fbd_region"]), on="sigungu_cd", how="left", validate="one_to_one")
SIDO = G.assign(sido=G.sigungu_cd // 1000).dissolve("sido")
def base(ax):
    SIDO.boundary.plot(ax=ax, color="#555555", linewidth=0.25); ax.set_axis_off(); ax.set_aspect("equal")
def divmap(ax, g, col, title, vmax=None, cls=None):
    x = np.abs(np.log(g[col])).max() if vmax is None else vmax; norm = TwoSlopeNorm(vmin=math.exp(-x), vcenter=1, vmax=math.exp(x))
    g.plot(column=col, cmap="RdBu_r", norm=norm, ax=ax, edgecolor="#a0a0a0", linewidth=0.05)
    if cls is not None:
        g[g[cls] == "hotspot"].boundary.plot(ax=ax, color="black", linewidth=0.45); g[g[cls] == "coldspot"].boundary.plot(ax=ax, color="black", linewidth=0.45, linestyle=(0, (2, 1.2)))
    base(ax); ax.set_title(title, fontsize=8, loc="left", fontweight="bold"); return norm
SP = pd.read_csv(os.path.join(MAINRES, "S3_spatial_OR.csv")); SP["fbd_region"] = SP.fbd_region.map(NFC)
EN = {"서울시": "Seoul", "부산시": "Busan", "대구시": "Daegu", "인천시": "Incheon", "광주시": "Gwangju", "대전시": "Daejeon", "울산시": "Ulsan", "세종시": "Sejong", "경기도": "Gyeonggi", "강원도": "Gangwon", "충청북도": "North Chungcheong", "충청남도": "South Chungcheong",
      "전라북도": "North Jeolla", "전라남도": "South Jeolla", "경상북도": "North Gyeongsang", "경상남도": "South Gyeongsang", "제주도": "Jeju"}
SP["province"] = SP.sido.map(NFC).map(EN); hc = SP[SP["class"] != "neutral"].merge(key, on="fbd_region").sort_values(["class", "province", "OR"], ascending=[False, True, False])
hc[["class", "province", "sido", "fbd_region", "sigungu_cd", "OR", "lo", "hi", "P_gt1", "SR_mean"]].to_csv(os.path.join(T, "hotcold_hot_cold_97_list.csv"), index=False, encoding="utf-8-sig")
cnt = hc.groupby(["class", "province"]).size().unstack(0).fillna(0).astype(int); cnt.to_csv(os.path.join(T, "hotcold_hot_cold_by_province.csv"), encoding="utf-8-sig")
print("[hotcold]", hc["class"].value_counts().to_dict())
sr = FO.groupby("fbd_region").SR.mean().rename("v").reset_index()
gname = G.merge(sr, left_on="region", right_on="fbd_region", how="left", validate="one_to_one"); gcode = cjoin(sr, "Figure1_A_selfreport_mean")
spn = G.merge(SP[["fbd_region", "OR"]], left_on="region", right_on="fbd_region", how="left", validate="one_to_one"); spc = cjoin(SP[["fbd_region", "OR", "class"]], "Figure1_B_spatialOR")
MS = pd.read_csv(os.path.join(MAINRES, "S14_measured_spatialOR.csv")); MS["fbd_region"] = MS.fbd_region.map(NFC); msc = cjoin(MS[["fbd_region", "OR", "class"]], "Figure1_D_measured_spatialOR")
e092 = pd.DataFrame([dict(panel="A self-report mean SR", max_abs_diff_name_vs_code=float((gname.sort_values("sigungu_cd").v.values - gcode.sort_values("sigungu_cd").v.values).__abs__().max())),
                     dict(panel="B spatial OR", max_abs_diff_name_vs_code=float(np.abs(spn.sort_values("sigungu_cd").OR.values - spc.sort_values("sigungu_cd").OR.values).max()))])
e092.to_csv(os.path.join(T, "namejoin_figure1_name_vs_code_join.csv"), index=False); print("[namejoin]", e092.to_dict("records"))
ef = os.path.join(OUT, "obsfit_obs_fit_resid_cells.csv")
if os.path.exists(ef):
    C = pd.read_csv(ef); C["fbd_region"] = C.fbd_region.map(NFC)
    A = C.groupby("sigungu_cd").agg(obs=("SR", "mean"), fit=("fit_SR", "mean"), res=("resid_logit", "mean")).reset_index(); A.to_csv(os.path.join(T, "obsfit_district_obs_fit_resid.csv"), index=False)
    g = cjoin(A, "obsfit_obs_fit_resid")
    fig, axs = plt.subplots(1, 3, figsize=(7.1, 3.3)); vmin, vmax = np.quantile(np.r_[g.obs, g.fit], [.02, .98])
    for ax, c, t in [(axs[0], "obs", "A  Observed, standardized %"), (axs[1], "fit", "B  Fitted (main model, without cell term) %")]:
        g.plot(column=c, cmap="YlOrRd", vmin=vmin, vmax=vmax, ax=ax, edgecolor="#a0a0a0", linewidth=0.05); base(ax); ax.set_title(t, fontsize=7.5, loc="left", fontweight="bold")
    sm = plt.cm.ScalarMappable(cmap="YlOrRd", norm=plt.Normalize(vmin, vmax)); fig.colorbar(sm, ax=axs[:2], orientation="horizontal", fraction=0.04, pad=0.02, label="Obesity prevalence, 2011–2024 mean (%)")
    x = np.abs(g.res).max(); g.plot(column="res", cmap="PuOr_r", vmin=-x, vmax=x, ax=axs[2], edgecolor="#a0a0a0", linewidth=0.05); base(axs[2]); axs[2].set_title("C  Residual (logit), mean", fontsize=7.5, loc="left", fontweight="bold")
    fig.colorbar(plt.cm.ScalarMappable(cmap="PuOr_r", norm=plt.Normalize(-x, x)), ax=axs[2], orientation="horizontal", fraction=0.04, pad=0.02, label="Observed − fitted, logit")
    fig.savefig(os.path.join(FIG, "FigS_obsfit_observed_fitted_residual.png"), dpi=300, facecolor="white"); plt.close(fig)
S4 = pd.read_csv(os.path.join(MAINRES, "S4_window_OR.csv")); IX = pd.read_csv(os.path.join(MAINRES, "S4_interaction.csv")); ixc = {(r.cov, r.window): bool(r.cred) for r in IX.itertuples()}
ixc = {(c, {"W2_2016-18": "W2_2016-18", "W3_2020-24": "W3_2020-24"}.get(w, w)): v for (c, w), v in ixc.items()}
WIN = list(dict.fromkeys(S4.window)); order = HYP + FORCED; nr, nc = len(order), len(WIN)
npw = P.groupby("window").size().to_dict()
fig, ax = plt.subplots(figsize=(4.0 + 1.75 * nc, 0.33 * nr + 1.1)); RED, BLUE, RL, BL = "#b2182b", "#2166ac", "#f4d0cb", "#d3dff0"; rows = []
for i, v in enumerate(order):
    for j, w in enumerate(WIN):
        r = S4[(S4.window == w) & (S4["cov"] == v)].iloc[0]; cred = bool(r.cred); pos = r.OR > 1
        ax.add_patch(Rectangle((j, nr - 1 - i), 1, 1, facecolor=(RL if pos else BL) if cred else "white", edgecolor=(RED if pos else BLUE) if cred else "#cccccc", linewidth=1.2 if cred else 0.4))
        ax.text(j + .5, nr - .5 - i, f"{r.OR:.3f} ({r.lo:.3f}–{r.hi:.3f})", ha="center", va="center", fontsize=7, color=((RED if pos else BLUE) if cred else "#555555"), fontweight="bold" if cred else "normal")
        dag = ixc.get((v, w), False)
        if dag: ax.text(j + .97, nr - 1 - i + .93, "†", ha="right", va="top", fontsize=7.5)
        rows.append(dict(cov=v, label=LAB.get(v, v), role="hypothesis" if v in HYP else "forced", window=w, OR=r.OR, lo=r.lo, hi=r.hi, cred=cred, interaction_vs_W1_cred=dag, WAIC=r.WAIC))
ax.axhline(nr - len(HYP), color="#333333", linewidth=0.8); ax.set_xlim(0, nc); ax.set_ylim(0, nr); ax.set_xticks(np.arange(nc) + .5)
ax.set_xticklabels([f"{w.split('_')[1].replace('-', '–20') if False else w.split('_')[1]}\n(n = {npw[w]:,} district-years)" for w in WIN]); ax.xaxis.tick_top()
ax.set_yticks(np.arange(nr) + .5); ax.set_yticklabels([LAB.get(v, v) for v in order][::-1], fontsize=7)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0); plt.tight_layout(); fig.savefig(os.path.join(FIG, "FigS_wingrid_window_grid.png"), dpi=300, facecolor="white"); plt.close(fig)
pd.DataFrame(rows).to_csv(os.path.join(T, "wingrid_window_grid_table.csv"), index=False)
wf = os.path.join(OUT, "W85_window_spatialOR.csv")
if os.path.exists(wf):
    W = pd.read_csv(wf); W["fbd_region"] = W.fbd_region.map(NFC); x = np.abs(np.log(W.OR)).max()
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 5.6))
    for j, w in enumerate(WIN):
        g = cjoin(W[W.window == w][["fbd_region", "OR", "P_gt1", "class"]], f"W85_{w}")
        divmap(axs[0, j], g, "OR", f"{'abc'[j]}  {w.split('_')[1]}", vmax=x, cls="class")
        g.plot(column="P_gt1", cmap=ListedColormap(["#2166AC", "#D1E5F0", "#FEE0B6", "#FDAE61", "#B2182B"]), norm=BoundaryNorm([0, .05, .5, .8, .95, 1.0001], 5), ax=axs[1, j], edgecolor="#a0a0a0", linewidth=0.05); base(axs[1, j])
    fig.colorbar(plt.cm.ScalarMappable(cmap="RdBu_r", norm=TwoSlopeNorm(vmin=math.exp(-x), vcenter=1, vmax=math.exp(x))), ax=axs[0, :], orientation="horizontal", fraction=0.04, pad=0.02, label="Spatial odds ratio (BYM2, refitted per window); solid = hotspot, dashed = coldspot")
    fig.colorbar(plt.cm.ScalarMappable(cmap=ListedColormap(["#2166AC", "#D1E5F0", "#FEE0B6", "#FDAE61", "#B2182B"]), norm=BoundaryNorm([0, .05, .5, .8, .95, 1.0], 5)), ax=axs[1, :], orientation="horizontal", fraction=0.04, pad=0.02, label="Posterior probability that the spatial OR exceeds 1")
    fig.savefig(os.path.join(FIG, "FigS_winmap_window_spatialOR_exceedance.png"), dpi=300, facecolor="white"); plt.close(fig)
YRS = sorted(P.year.unique()); obs = FO[FO.year.isin(YRS)].merge(key, on="fbd_region")
sm_ = pd.read_csv(ef) if os.path.exists(ef) else None
def yearly(col_df, col, fn, label):
    vmin, vmax = np.quantile(col_df[col], [.02, .98]); fig, axs = plt.subplots(3, 5, figsize=(7.1, 5.4))
    for k, ax in enumerate(axs.flat):
        if k >= len(YRS): ax.set_axis_off(); continue
        y = YRS[k]; g = cjoin(col_df[col_df.year == y][["sigungu_cd", col]], f"{fn}_{y}")
        g.plot(column=col, cmap="YlOrRd", vmin=vmin, vmax=vmax, ax=ax, edgecolor="#a0a0a0", linewidth=0.04); base(ax); ax.set_title(str(y), fontsize=7.5)
    fig.colorbar(plt.cm.ScalarMappable(cmap="YlOrRd", norm=plt.Normalize(vmin, vmax)), ax=axs, orientation="horizontal", fraction=0.03, pad=0.02, label=label)
    fig.savefig(os.path.join(FIG, fn + ".png"), dpi=300, facecolor="white"); plt.close(fig)
yearly(obs, "SR", "FigS_yearmap_yearly_observed_SR", "Self-reported obesity, age-sex-standardized % (direct survey estimate; 2019 not surveyed for this item)")
if sm_ is not None: yearly(sm_.assign(sigungu_cd=sm_.sigungu_cd.astype(int)), "fit_SR", "FigS_yearmap_yearly_smoothed_SR", "Model-smoothed prevalence, % (main model fitted values without the cell term)")
import pymannkendall as pmk
from statsmodels.stats.multitest import multipletests
gl = open(os.path.join(RUN, "04_분석/OBES/graph_593.graph")).read().split("\n"); nbg = {}
for line in gl[1:]:
    t = line.split()
    if len(t) >= 2: nbg[int(t[0])] = [int(q) for q in t[2:]]
area = {r: i + 1 for i, r in enumerate(G.region)}
def gistar(x, nbo):
    nn = len(x); xb = x.mean(); s = math.sqrt((x ** 2).sum() / nn - xb ** 2); z = np.zeros(nn)
    for i in range(nn):
        nbr = [i] + nbo[i]; k = len(nbr); w = 1.0 / k
        z[i] = (x[nbr].sum() * w - xb) / (s * math.sqrt((nn * (k * w * w) - 1.0) / (nn - 1)))
    return z
def last_run(b):
    c = 0
    for v in b[::-1]:
        if v: c += 1
        else: break
    return c
thr = math.ceil(0.9 * len(YRS))
def esri(zv, tau, p):
    hot = zv > 1.96; cold = zv < -1.96; nh, ncl = hot.sum(), cold.sum(); up = (p < 0.05) and tau > 0; dn = (p < 0.05) and tau < 0
    if hot[-1]:
        if nh >= thr: return "Intensifying hot" if up else ("Diminishing hot" if dn else "Persistent hot")
        if nh == 1: return "New hot"
        if last_run(hot) >= 2 and last_run(hot) == nh: return "Consecutive hot"
        return "Oscillating hot" if cold.any() else "Sporadic hot"
    if cold[-1]:
        if ncl >= thr: return "Intensifying cold" if dn else ("Diminishing cold" if up else "Persistent cold")
        if ncl == 1: return "New cold"
        if last_run(cold) >= 2 and last_run(cold) == ncl: return "Consecutive cold"
        return "Oscillating cold" if hot.any() else "Sporadic cold"
    if nh >= thr: return "Historical hot"
    if ncl >= thr: return "Historical cold"
    return "No pattern"
ehs = []
for src, df, col in [("observed_logit", FO[FO.year.isin(YRS)].assign(v=lambda q: q.y_logit), "v")] + ([("smoothed_logit", sm_.assign(v=sm_.fit_logit_nocell), "v")] if sm_ is not None else []):
    df = df.copy(); df["fbd_region"] = df.fbd_region.map(NFC); M_ = df.pivot(index="year", columns="fbd_region", values=col)
    regs = sorted(M_.columns, key=lambda r: area[r]); pos = {area[r]: k for k, r in enumerate(regs)}; nbo = [[pos[a] for a in nbg[area[r]] if a in pos] for r in regs]
    Z = np.vstack([gistar(M_.loc[y, regs].values.astype(float), nbo) for y in YRS]); yrs = np.array(YRS, float)
    orig = [kendalltau(yrs, Z[:, k]) for k in range(len(regs))]; p_o = np.array([o[1] for o in orig]); p_f = multipletests(p_o, method="fdr_bh")[1]
    hr = [pmk.hamed_rao_modification_test(Z[:, k]) for k in range(len(regs))]
    for k, r in enumerate(regs):
        ehs.append(dict(source=src, fbd_region=r, EHSA_original=esri(Z[:, k], orig[k][0], p_o[k]), EHSA_hamed_rao=esri(Z[:, k], hr[k].Tau, hr[k].p), EHSA_mk_bh_fdr=esri(Z[:, k], orig[k][0], p_f[k]),
                        n_hot_years=int((Z[:, k] > 1.96).sum()), n_cold_years=int((Z[:, k] < -1.96).sum()), mk_tau=round(orig[k][0], 3), mk_p=round(p_o[k], 4), mk_q_bh=round(p_f[k], 4), hr_p=round(hr[k].p, 4)))
EH = pd.DataFrame(ehs); EH.to_csv(os.path.join(T, "ehsa_EHSA_district.csv"), index=False, encoding="utf-8-sig")
cats = sorted(set(EH.EHSA_original) | set(EH.EHSA_hamed_rao) | set(EH.EHSA_mk_bh_fdr)); sm_rows = []
for src, d in EH.groupby("source"):
    for c in cats: sm_rows.append(dict(source=src, category=c, original_MK=int((d.EHSA_original == c).sum()), hamed_rao_MK=int((d.EHSA_hamed_rao == c).sum()), MK_BH_FDR=int((d.EHSA_mk_bh_fdr == c).sum())))
    sm_rows.append(dict(source=src, category="changed_vs_original", original_MK=0, hamed_rao_MK=int((d.EHSA_hamed_rao != d.EHSA_original).sum()), MK_BH_FDR=int((d.EHSA_mk_bh_fdr != d.EHSA_original).sum())))
ES = pd.DataFrame(sm_rows); ES.to_csv(os.path.join(T, "ehsa_EHSA_summary.csv"), index=False); print(ES[(ES.original_MK + ES.hamed_rao_MK + ES.MK_BH_FDR) > 0].to_string())
pal = {"Persistent hot": "#b2182b", "Intensifying hot": "#67001f", "Diminishing hot": "#ef8a62", "Consecutive hot": "#f4a582", "Oscillating hot": "#fddbc7", "Sporadic hot": "#fee8d8", "New hot": "#d6604d", "Historical hot": "#c51b7d",
       "Persistent cold": "#2166ac", "Intensifying cold": "#053061", "Diminishing cold": "#67a9cf", "Consecutive cold": "#92c5de", "Oscillating cold": "#d1e5f0", "Sporadic cold": "#e0ecf4", "New cold": "#4393c3", "Historical cold": "#542788", "No pattern": "#f0f0f0"}
srcs = list(EH.source.unique()); fig, axs = plt.subplots(1, len(srcs), figsize=(3.6 * len(srcs), 4.0)); axs = np.atleast_1d(axs)
for ax, src in zip(axs, srcs):
    g = cjoin(EH[EH.source == src][["fbd_region", "EHSA_mk_bh_fdr"]], f"EHSA_{src}"); g["c"] = g.EHSA_mk_bh_fdr.map(pal)
    g.plot(color=g.c, ax=ax, edgecolor="#a0a0a0", linewidth=0.05); base(ax); ax.set_title(f"{'AB'[srcs.index(src)]}  Emerging hot spots, {src.replace('_', ' ')} (MK with BH-FDR)", fontsize=7, loc="left", fontweight="bold")
used = [c for c in pal if c in set(EH.EHSA_mk_bh_fdr)]; fig.legend(handles=[Patch(facecolor=pal[c], edgecolor="#888", label=f"{c}") for c in used], loc="lower center", ncol=4, fontsize=6.5, frameon=False)
plt.subplots_adjust(bottom=0.16); fig.savefig(os.path.join(FIG, "FigS_ehsa_EHSA.png"), dpi=300, facecolor="white"); plt.close(fig)
ROB6 = ["walking_practice_rate", "tertiary_educated_share", "perceived_stress_rate", "high_risk_drinking_rate_drinkers", "ext__restaurant_per10k", "current_smoking_rate"]
D = P.merge(FO[["year", "fbd_region", "y_logit"]], on=["year", "fbd_region"]); tw = lambda s: s - D.groupby("fbd_region")[s.name].transform("mean") - D.groupby("year")[s.name].transform("mean") + s.mean()
D["s_dev"] = tw(D.y_logit); D["period"] = np.where(D.year <= 2018, "2011–2018", "2020–2024"); rows = []
fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.8), sharey=True)
for ax, v in zip(axs.flat, ROB6):
    D["x_dev"] = tw(D[v + "_raw"].rename(v + "_raw")) if False else (D[v + "_raw"] - D.groupby("fbd_region")[v + "_raw"].transform("mean") - D.groupby("year")[v + "_raw"].transform("mean") + D[v + "_raw"].mean())
    for per, col in [("2011–2018", "#4d4d4d"), ("2020–2024", "#d6604d")]:
        d = D[D.period == per]; q = pd.qcut(d.x_dev, 10, labels=False, duplicates="drop"); gq = d.groupby(q).agg(x=("x_dev", "mean"), s=("s_dev", "mean"), se=("s_dev", lambda z: z.std(ddof=1) / np.sqrt(len(z))))
        sl = np.polyfit(d.x_dev, d.s_dev, 1)[0]; ax.errorbar(gq.x, gq.s, yerr=1.96 * gq.se, fmt="o", ms=3, color=col, capsize=1.5, lw=0.7, label=per)
        xs = np.linspace(d.x_dev.quantile(.02), d.x_dev.quantile(.98), 10); ax.plot(xs, sl * xs, color=col, lw=0.8)
        for i_, r in gq.iterrows(): rows.append(dict(cov=v, period=per, decile=int(i_) + 1, x_dev_mean=round(r.x, 4), logit_dev_mean=round(r.s, 4), se=round(r.se, 4), ols_slope=round(sl, 5)))
    ax.axhline(0, color="grey", lw=0.4); ax.axvline(0, color="grey", lw=0.4); ax.set_title(LAB.get(v, v), fontsize=7)
axs[0, 0].legend(fontsize=6.5, frameon=False); fig.supylabel("Within-district, within-year deviation in logit prevalence", fontsize=7.5); fig.supxlabel("Deviation in covariate (natural units after transform; district and year means removed)", fontsize=7.5)
plt.tight_layout(); fig.savefig(os.path.join(FIG, "FigS_ehsa_twoway_centred_robust6.png"), dpi=300, facecolor="white"); plt.close(fig)
pd.DataFrame(rows).to_csv(os.path.join(T, "ehsa_twoway_centred_bins.csv"), index=False)
pd.DataFrame(chk).to_csv(os.path.join(T, "mapjoin_map_join_checks.csv"), index=False); print("[mapjoin] checks", len(chk), "all set_equal", all(c["set_equal"] for c in chk), "all 1:1", all(c["one_to_one"] for c in chk), "names agree", all(c["code_join_region_name_agrees"] for c in chk))
