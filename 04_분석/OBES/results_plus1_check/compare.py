# Compares the stored main model with the rerun and the alternative green-zone transforms. Usage: python3 compare.py <this folder>
import sys, pandas as pd, numpy as np
D = sys.argv[1]; R = D + "/../results_std3"
orig = pd.read_csv(R + "/S3_main_OR.csv"); of = pd.read_csv(R + "/S3_main_fit.csv")
fits = ["orig_rerun", "A_halfmin", "B_zeroNA_fill", "B2_zeroNA_drop"]
T = {"stored": orig[["cov", "OR", "lo", "hi", "cred", "label_en", "role"]]}
S = {"stored": dict(WAIC=of.WAIC[0], DIC=of.DIC[0], LCPO=of.LCPO[0], phi=of.phi[0], n_obs=of.n_obs[0])}
for f in fits:
    T[f] = pd.read_csv(f"{D}/fit_{f}_OR.csv"); s = pd.read_csv(f"{D}/fit_{f}_stat.csv").iloc[0]; S[f] = s.to_dict()
rows = []
for _, r in orig.iterrows():
    row = dict(cov=r["cov"], label=r.label_en, role=r.role)
    for k, t in T.items():
        x = t[t["cov"] == r["cov"]].iloc[0]; row[k] = f"{x.OR:.4f} ({x.lo:.4f}–{x.hi:.4f}){'*' if x.cred else ''}"; row[k + "_cred"] = bool(x.cred)
        for q in ["OR", "lo", "hi"]: row[f"{k}_{q}"] = x[q]
    rows.append(row)
W = pd.DataFrame(rows); W.to_csv(D + "/plus1_compare_all.csv", index=False)
logd = lambda a, b: max(abs(np.log(W[f"{a}_{q}"]) - np.log(W[f"{b}_{q}"])).max() for q in ["OR", "lo", "hi"])
logd_ex = lambda a, b, ex: max(abs(np.log(W.loc[W["cov"] != ex, f"{a}_{q}"]) - np.log(W.loc[W["cov"] != ex, f"{b}_{q}"])).max() for q in ["OR", "lo", "hi"])
G = "eloc__dt_1yl20421e__1___1"
print("max |dlog| stored vs rerun:", logd("stored", "orig_rerun"))
for f in fits[1:]:
    print(f, "max|dlog| vs rerun all", round(logd("orig_rerun", f), 5), "excl green", round(logd_ex("orig_rerun", f, G), 5),
          "cred changes", W.loc[W[f + "_cred"] != W["orig_rerun_cred"], "cov"].tolist())
print("stored vs rerun cred changes", W.loc[W["stored_cred"] != W["orig_rerun_cred"], "cov"].tolist())
print(pd.DataFrame(S).T)
pd.DataFrame(S).T.to_csv(D + "/plus1_fitstats.csv")
print(W[["cov", "stored", "orig_rerun", "A_halfmin", "B_zeroNA_fill", "B2_zeroNA_drop"]].to_string())
W["min_dist_bound_to1"] = W[[f"{k}_{q}" for k in T for q in ["lo", "hi"]]].sub(1).abs().min(axis=1)
print(W[["cov", "min_dist_bound_to1"]].sort_values("min_dist_bound_to1").head(6).to_string())
