#!/usr/bin/env python3
# Pre-model data gate for a district-year panel: whole-year zero/missing/constant columns, district-wise structural zeros, columns explained by year (R2 >= 0.8), 2x year-mean jumps, rates above 100. Usage: python3 data_gate.py <panel.csv> [--year year] [--region region]
"""분석 전 데이터 게이트: 지역-연도 패널의 공변량을 모형에 넣기 전에 점검한다.
사용: python3 data_gate.py <panel.csv> [--year year] [--region region]
점검: (1) 특정 연도에 전 지역이 0/결측/상수  (2) 매년 같은 지역들이 0 (지역형 구조적 0)
      (3) 연도만으로 분산이 설명되는 열(R^2>=0.8)  (4) 연도 간 평균 급변(2배 이상)  (5) 비율형인데 100 초과
"""
import sys, argparse, numpy as np, pandas as pd
ap = argparse.ArgumentParser(); ap.add_argument("panel"); ap.add_argument("--year", default="year"); ap.add_argument("--region", default="region")
a = ap.parse_args()
try: d = pd.read_csv(a.panel, low_memory=False, encoding="utf-8-sig")
except UnicodeDecodeError: d = pd.read_csv(a.panel, low_memory=False, encoding="cp949")
Y, R = a.year, a.region
d[Y] = pd.to_numeric(d[Y], errors="coerce"); yrs = sorted(d[Y].dropna().unique()); flags = []
for c in d.columns:
    if c in (Y, R): continue
    x = pd.to_numeric(d[c], errors="coerce")
    if x.notna().sum() < 30 or x.nunique(dropna=True) < 2: continue
    g = pd.DataFrame({"y": d[Y], "r": d[R] if R in d else 0, "x": x})
    by = g.groupby("y")["x"]
    zero = [int(y) for y, s in by if (s.fillna(0) == 0).mean() >= 0.98]
    const = [int(y) for y, s in by if s.nunique(dropna=True) <= 1 and int(y) not in zero]
    if zero: flags.append((c, "연도 전체 0/결측", zero))
    if const: flags.append((c, "연도 전체 상수", const))
    if R in d:
        zr = g.assign(z=(g.x.fillna(0) == 0)).groupby("r")["z"].mean()
        n = int((zr == 1).sum())
        if 0 < n < zr.size * 0.6 and n >= 5 and not zero: flags.append((c, "전 연도 0인 지역", f"{n}곳"))
    m = by.mean(); tot = ((x - x.mean()) ** 2).sum()
    if tot > 0:
        r2 = float(((g["y"].map(m) - x.mean()) ** 2).sum() / tot)
        if r2 >= 0.8: flags.append((c, "연도가 분산 대부분 설명", f"R2={r2:.2f}"))
    mm = m[m > 0]
    if len(mm) >= 2 and (mm.max() / mm.min()) >= 2: flags.append((c, "연도 간 평균 2배 이상 차이", f"{mm.idxmin():.0f}:{mm.min():.3g} vs {mm.idxmax():.0f}:{mm.max():.3g}"))
    if any(k in c for k in ("율", "률", "rate", "pct", "ratio", "%")) and x.max() > 100: flags.append((c, "비율형인데 100 초과", f"max={x.max():.1f}"))
print(f"패널 {len(d):,}행, 연도 {yrs[0]:.0f}-{yrs[-1]:.0f}, 점검 열 {d.shape[1]}")
if not flags: print("통과: 경고 없음")
else:
    print(f"경고 {len(flags)}건 — 모형 적합 전에 원자료·출처에서 확인할 것")
    for c, k, v in flags: print(f"  [{k}] {c}: {v}")
    sys.exit(1)
