# Writes the job lists for runner.sh (phases A, B, C) for the current panel.
# Usage: python3 make_jobs.py <run root (parent of 01_스크립트)> <boundary file> <contact graph> <out folder for .jobs>
# Results: <run root>/04_분석/OBES/results_std3, results_extra (needs extra_inputs.csv from 82_extra_inputs.py), results_plus1_check (needs green_variants.csv).
import sys, os
ROOT, SHP, GRAPH, JD = sys.argv[1:5]
S = f"{ROOT}/01_스크립트"; OB = f"{ROOT}/04_분석/OBES"
BASE = f"ST_SHP={SHP} ST_GRAPH={GRAPH} ST_THREADS=1"
E = f"{BASE} ST_RES={OB}/results_std3 ST_OLD65_COL=pop65_share_raw_fbd5"
R = f"Rscript {S}/81fh_obes.R OBES"
EX = f"{BASE} EX_RES={OB}/results_extra"; PC = f"{BASE} PC_RES={OB}/results_plus1_check"
J = {"A": [], "B": [], "C": []}
add = lambda ph, nm, to, cmd: J[ph].append(f"{nm}\t{to}\t{cmd}")
add("A", "S0_S2_S3", 90, f"env ST_FH_ONLY=S0,S2,S3 {E} {R}")
add("A", "S7_specs", 30, f"env ST_FH_ONLY=S7 ST_MV_STOP=1 {E} {R}")
for r in ["1-50", "51-100", "101-150", "151-200"]:
    add("B", f"S7_{r}", 360, f"env ST_FH_ONLY=S7 ST_MV_RANGE={r} {E} {R}")
add("B", "sensA_S4_S5_S6_S8_S9", 150, f"env ST_FH_ONLY=S4,S5,S6,S8,S9 {E} {R}")
add("B", "sensB_S10_S11_S12_S13", 150, f"env ST_FH_ONLY=S10,S11,S12,S13 {E} {R}")
add("B", "sensC_S14_S16_S18_S19", 150, f"env ST_FH_ONLY=S14,S16,S18,S19 {E} {R}")
add("B", "newA_S20_S21_S22_S23", 150, f"env ST_FH_ONLY=S20,S21,S22,S23 {E} {R}")
add("B", "newB_S24_S25_S26", 150, f"env ST_FH_ONLY=S24,S25,S26 {E} {R}")
add("B", "newC_S27_S28_S29_S31", 150, f"env ST_FH_ONLY=S27,S28,S29,S31 {E} {R}")
add("B", "extra_ABCDE", 150, f"env EX_ONLY=orig_rerun,A1_icar,A2_nospatial,B_twoway,C1_km2,C2_both,D1_bnd_terr,D2_bnd_any,E1_drop_mvpa,E1cc_mvpa,E2_drop_fac,E2cc_fac,E3_drop_tert,E3cc_tert,E4_drop_bl {EX} Rscript {S}/83_extra_fits.R OBES")
add("B", "extra_F_MI", 150, f"env EX_ONLY='F_MI_*' {EX} Rscript {S}/83_extra_fits.R OBES")
add("B", "plus1", 90, f"env {PC} Rscript {S}/84_refit_green_plus1.R OBES")
add("C", "S7merge_S17_S44", 120, f"env ST_FH_ONLY=S7,S17,S44 {E} {R}")
os.makedirs(JD, exist_ok=True)
for ph, L in J.items():
    open(f"{JD}/{ph}.jobs", "w").write("\n".join(L) + "\n"); print(ph, len(L))
