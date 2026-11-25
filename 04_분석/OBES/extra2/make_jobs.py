# Writes the round-2 job list (stages J1-J7) for extra2/runner.sh.
# Usage: python3 make_jobs.py <run root (parent of 01_스크립트)> <boundary file> <contact graph> <out folder for .jobs>
import sys, os
ROOT, SHP, GRAPH, JD = sys.argv[1:5]
OB = f"{ROOT}/04_분석/OBES"; E2 = f"{OB}/extra2"
ENV = (f"E2_RES={OB}/results_extra2 E2_MAINRES={OB}/results_std3 E2_RUNSRC={ROOT}/01_스크립트/81fh_obes.R "
       f"ST_SHP={SHP} ST_GRAPH={GRAPH} ST_THREADS=1 ST_OLD65_COL=pop65_share_raw_fbd5")
L = [f"e2_{j}\t90\tenv {ENV} E2_STAGE={j} Rscript {E2}/e2_fits.R OBES" for j in ["J1", "J2", "J3", "J4", "J5", "J6", "J7"]]
os.makedirs(JD, exist_ok=True); open(f"{JD}/E2.jobs", "w").write("\n".join(L) + "\n"); print("E2", len(L))
