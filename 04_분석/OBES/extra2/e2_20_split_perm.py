#!/usr/bin/env python3
# Selection-corrected test of the 2018 split: two-way fixed-effects weighted least squares interaction t for each candidate split year, with year-label permutation (B = 999). Environment: E2_RUN (run root with 02_패널/ and 04_분석/OBES/ model inputs), E2_OUT (results folder, default 04_분석/OBES/results_extra2).
import os, json, numpy as np, pandas as pd, unicodedata
NFC = lambda s: unicodedata.normalize("NFC", str(s))
OBD = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.environ.get("E2_OUT", os.path.join(OBD, "results_extra2")); RUN = os.environ["E2_RUN"]
P = pd.read_csv(os.path.join(RUN, "04_분석/OBES/panel_OBES_std3.csv"), encoding="utf-8-sig"); P["fbd_region"] = P.fbd_region.map(NFC)
FO = pd.read_csv(os.path.join(RUN, "02_패널/OBES_fh_outcomes.csv"), encoding="utf-8-sig"); FO["fbd_region"] = FO.fbd_region.map(NFC)
CFG = json.load(open(os.path.join(RUN, "04_분석/OBES/config_OBES_std3.json"))); HYP = CFG["sets_after_vif"][CFG["rules"]["main_set"]]; MAINV = list(dict.fromkeys(HYP + CFG["forced"]))
ROB6 = ["walking_practice_rate", "tertiary_educated_share", "perceived_stress_rate", "high_risk_drinking_rate_drinkers", "ext__restaurant_per10k", "current_smoking_rate"]
D = P.merge(FO[["year", "fbd_region", "y_logit", "v_logit"]], on=["year", "fbd_region"], validate="one_to_one").sort_values(["year", "fbd_region"]).reset_index(drop=True)
YRS = sorted(D.year.unique()); K = [k for k in YRS if sum(y < k for y in YRS) >= 2 and sum(y >= k for y in YRS) >= 2]; assert 2018 in K
w = 1 / D.v_logit.values; sw = np.sqrt(w); y = D.y_logit.values
reg = pd.get_dummies(D.fbd_region, drop_first=True).values.astype(float); yd = pd.get_dummies(D.year, drop_first=True).values.astype(float)
X0 = np.column_stack([np.ones(len(D)), reg, yd, D[MAINV].values]); H = D[HYP].values; ir = [HYP.index(v) for v in ROB6]
def tstats(post):
    X = np.column_stack([X0, H * post[:, None]]); Xw = X * sw[:, None]; yw = y * sw
    b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); r = yw - Xw @ b; dof = len(y) - np.linalg.matrix_rank(Xw); s2 = r @ r / dof
    XtX = np.linalg.pinv(Xw.T @ Xw); se = np.sqrt(np.diag(XtX) * s2); j = X0.shape[1] + np.array(ir)
    return b[j] / se[j], b[j], se[j]
yr = D.year.values
def all_k(ymap):
    py = np.vectorize(ymap.get)(yr); return np.vstack([tstats((py >= k).astype(float))[0] for k in K])
ident = {t: t for t in YRS}; T_obs = all_k(ident); ik = K.index(2018); t18 = T_obs[ik]
_, b18, se18 = tstats((yr >= 2018).astype(float))
rng = np.random.default_rng(20261129); B = 999; Tp = np.empty((B, len(K), len(ROB6)))
for bb in range(B):
    perm = rng.permutation(YRS); Tp[bb] = all_k(dict(zip(YRS, perm)))
from scipy.stats import norm
rows = []
for c, v in enumerate(ROB6):
    a = abs(t18[c]); rows.append(dict(cov=v, split_year=2018, beta_interaction_FE=round(b18[c], 5), t_2018=round(t18[c], 3), p_normal=round(2 * norm.sf(a), 4),
        p_perm_fixed2018=round((1 + np.sum(np.abs(Tp[:, ik, c]) >= a)) / (B + 1), 4), p_perm_max_over_splits=round((1 + np.sum(np.abs(Tp[:, :, c]).max(1) >= a)) / (B + 1), 4),
        rank_of_2018_among_splits=int((np.abs(T_obs[:, c]) > a).sum() + 1), n_candidate_splits=len(K), max_abs_t_any_split=round(np.abs(T_obs[:, c]).max(), 3), split_with_max=int(K[int(np.abs(T_obs[:, c]).argmax())])))
R = pd.DataFrame(rows); fam = (1 + np.sum(np.abs(Tp).max((1, 2)) >= np.abs(t18).max())) / (B + 1)
R["p_family_max_over_splits_and_6covs"] = round(fam, 4); R["B"] = B; R["seed"] = 20261129; R["candidate_splits"] = ",".join(map(str, K))
R.to_csv(os.path.join(OUT, "splitperm_split_selection_permutation.csv"), index=False)
pd.DataFrame(T_obs, index=K, columns=ROB6).round(3).to_csv(os.path.join(OUT, "splitperm_FE_t_by_split.csv"))
print(R.drop(columns=["candidate_splits"]).to_string()); print("family p", fam)
