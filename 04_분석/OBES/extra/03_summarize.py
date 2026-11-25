#!/usr/bin/env python3
# Summary of the additional fits A-F: hypothesis OR and credibility table, WAIC, two-way decomposition against S9, Rubin pooling. Environment: OBES_RES, OBES_EXTRA.
import os, glob, numpy as np, pandas as pd
from decimal import Decimal as _D, ROUND_HALF_UP as _RHU
import builtins as _bi
def round(x, nd=None):
    try:
        q = _D(repr(float(x))).quantize(_D(1).scaleb(-(nd or 0)), rounding=_RHU)
    except Exception:
        return _bi.round(x, nd) if nd is not None else _bi.round(x)
    return float(q) if nd else int(q)

from scipy import stats
OB = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
R3 = os.environ.get("OBES_RES", os.path.join(OB, "results_std3")); OUT = os.environ.get("OBES_EXTRA", os.path.join(OB, "results_extra")); FD = os.path.join(OUT, "fits")
S3 = pd.read_csv(os.path.join(R3, "S3_main_OR.csv")); HYP = list(S3[S3.role == "hypothesis"]["cov"]); FORCED = list(S3[S3.role == "forced"]["cov"])
main = S3.set_index("cov")
def fmt(r): return f"{r.OR:.3f} ({r.lo:.3f}–{r.hi:.3f}){'*' if r.cred else ''}"
FIT = ["A1_icar", "A2_nospatial", "C1_km2", "C2_both", "D1_bnd_terr", "D2_bnd_any", "E1_drop_mvpa", "E1cc_mvpa", "E2_drop_fac", "E2cc_fac", "E3_drop_tert", "E3cc_tert", "E4_drop_bl"]
long = []; stat = []
for f in ["orig_rerun"] + FIT:
    a = pd.read_csv(os.path.join(FD, f"fit_{f}_OR.csv")); a["fit"] = f; long.append(a)
    s = pd.read_csv(os.path.join(FD, f"fit_{f}_stat.csv")); stat.append(s)
L = pd.concat(long, ignore_index=True); ST = pd.concat(stat, ignore_index=True)
ST["dWAIC_vs_main"] = (ST.WAIC - ST.loc[ST.fit == "orig_rerun", "WAIC"].iloc[0]).round(1)
ST.to_csv(os.path.join(OUT, "ABCDE_fitstats.csv"), index=False)
rows = []; changes = []
for v in HYP + ["pub_km2_z", "vil_km2_z"]:
    row = {"cov": v, "main": fmt(main.loc[v]) if v in main.index else "—"}
    for f in FIT:
        r = L[(L.fit == f) & (L["cov"] == v)]
        row[f] = fmt(r.iloc[0]) if len(r) else "—"
        if len(r) and v in main.index and bool(r.iloc[0].cred) != bool(main.loc[v, "cred"]):
            changes.append(dict(analysis=f, cov=v, main=fmt(main.loc[v]), alt=fmt(r.iloc[0]), change=("credible→not" if main.loc[v, "cred"] else "not→credible")))
        if len(r) and v in main.index and np.sign(np.log(r.iloc[0].OR)) != np.sign(np.log(main.loc[v, "OR"])) and (r.iloc[0].cred or main.loc[v, "cred"]):
            changes.append(dict(analysis=f, cov=v, main=fmt(main.loc[v]), alt=fmt(r.iloc[0]), change="sign flip"))
    rows.append(row)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "ABCDE_OR_compare.csv"), index=False)
L.to_csv(os.path.join(OUT, "ABCDE_OR_long.csv"), index=False)
B2 = pd.read_csv(os.path.join(FD, "fit_B_twoway_OR.csv")); S9 = pd.read_csv(os.path.join(R3, "S9_between_within.csv"))
br = []
for v in HYP:
    s9b = S9[(S9["cov"] == v) & (S9.component == "between")].iloc[0]; s9w = S9[(S9["cov"] == v) & (S9.component == "within")].iloc[0]
    tb = B2[(B2["cov"] == v) & (B2.component == "between")].iloc[0]; tw = B2[(B2["cov"] == v) & (B2.component == "within_twoway")].iloc[0]
    sdw = tw.sd_component
    br.append(dict(cov=v, main=fmt(main.loc[v]), S9_between=fmt(s9b), S9_within_oneway=fmt(s9w), B_between=fmt(tb), B_within_twoway=fmt(tw),
                   within_twoway_SD=round(sdw, 4), B_within_per_withinSD=f"{np.exp(np.log(tw.OR)*sdw):.4f} ({np.exp(np.log(tw.lo)*sdw):.4f}–{np.exp(np.log(tw.hi)*sdw):.4f})",
                   within_cred_change=("same" if bool(s9w.cred) == bool(tw.cred) else f"S9 {s9w.cred} → two-way {tw.cred}"),
                   between_cred_change=("same" if bool(s9b.cred) == bool(tb.cred) else f"S9 {s9b.cred} → two-way {tb.cred}")))
pd.DataFrame(br).to_csv(os.path.join(OUT, "B_twoway_vs_S9.csv"), index=False)
mi = sorted(glob.glob(os.path.join(FD, "fit_F_MI_*_OR.csv"))); Mn = len(mi)
A = pd.concat([pd.read_csv(p).assign(m=i) for i, p in enumerate(mi)], ignore_index=True)
nu_com = 2977 - 19; rr = []
for v in HYP + FORCED:
    x = A[A["cov"] == v]; q = x.beta_mean.values; u = x.beta_sd.values ** 2
    Q = q.mean(); U = u.mean(); Bv = q.var(ddof=1); T = U + (1 + 1 / Mn) * Bv; lam = (1 + 1 / Mn) * Bv / T
    nu_old = (Mn - 1) / lam ** 2 if lam > 0 else np.inf; nu_obs = (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lam); nu = 1 / (1 / nu_old + 1 / nu_obs)
    t = stats.t.ppf(0.975, nu); lo, hi = Q - t * np.sqrt(T), Q + t * np.sqrt(T)
    rr.append(dict(cov=v, role="hypothesis" if v in HYP else "forced", M=Mn, OR=round(np.exp(Q), 4), lo=round(np.exp(lo), 4), hi=round(np.exp(hi), 4), cred=bool(lo > 0 or hi < 0),
                   main=fmt(main.loc[v]), main_cred=bool(main.loc[v, "cred"]), FMI=round(lam, 3), df=round(nu, 1), n_cred_of_M=int(x.cred.sum()), OR_min_m=round(x.OR.min(), 4), OR_max_m=round(x.OR.max(), 4)))
RB = pd.DataFrame(rr); RB.to_csv(os.path.join(OUT, "F_MI_rubin.csv"), index=False)
for _, r in RB[RB.role == "hypothesis"].iterrows():
    if r.cred != r.main_cred: changes.append(dict(analysis="F_MI_rubin", cov=r["cov"], main=r.main, alt=f"{r.OR:.3f} ({r.lo:.3f}–{r.hi:.3f}){'*' if r.cred else ''}", change=("credible→not" if r.main_cred else "not→credible")))
for b in br:
    pass
pd.DataFrame(changes).to_csv(os.path.join(OUT, "judgement_changes.csv"), index=False)
print(ST.to_string()); print(pd.DataFrame(changes).to_string()); print(RB.to_string()); print(pd.DataFrame(br)[["cov", "S9_within_oneway", "B_within_twoway", "B_between", "within_cred_change", "between_cred_change"]].to_string())
