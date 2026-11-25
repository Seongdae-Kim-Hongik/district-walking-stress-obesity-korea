# Year-anomaly gate on the district-year master table: years where a variable's district median changes by <0.5x or >2x, or the number of districts with data halves or falls below 100. Environment: ST_DATA_ROOT (folder with 02_패널/ and 03_문서/).
"""연도 이상 게이트(FBD 563 패널 year_anomaly_gate 와 같은 취지): 변수마다 연도별 시군구 중앙값 시계열을 보고
(a) 전년 대비 중앙값 비율이 0.5 미만 또는 2 초과인 해, (b) 비결측 시군구 수가 전년 대비 절반 이하로 줄거나 100 미만인 해를 표시.
산출: 02_패널/MASTER_연도이상_게이트.csv (변수, 연도, 유형, 값), 03_문서/03_연도이상_게이트_요약.md"""
import os, numpy as np, pandas as pd
ROOT = os.environ.get("ST_DATA_ROOT", ".")
m = pd.read_parquet(os.path.join(ROOT, "02_패널", "MASTER_시군구_연도_wide.parquet"))
m = m[(m.year >= 1998) & (m.year <= 2026)]
meta = {"year", "key", "sido_std", "sigungu", "fbd_region", "population", "boundary_flag"}
vars_ = [c for c in m.columns if c not in meta and pd.api.types.is_numeric_dtype(m[c])]
med = m.groupby("year")[vars_].median(); cnt = m.groupby("year")[vars_].count()
rows = []
for v in vars_:
    s = med[v].dropna(); n = cnt[v]
    yrs = s.index.tolist()
    for i in range(1, len(yrs)):
        y0, y1 = yrs[i - 1], yrs[i]
        if y1 - y0 != 1: continue
        a, b = s[y0], s[y1]
        if a and b and a != 0 and (abs(b / a) > 2 or abs(b / a) < 0.5):
            rows.append({"column": v, "year": y1, "type": "중앙값 급변", "prev": a, "value": b, "ratio": round(b / a, 3)})
        if n[y0] >= 100 and n[y1] < n[y0] * 0.5:
            rows.append({"column": v, "year": y1, "type": "시군구 수 급감", "prev": n[y0], "value": n[y1], "ratio": round(n[y1] / n[y0], 3)})
    for y in yrs:
        if 0 < n[y] < 100: rows.append({"column": v, "year": y, "type": "시군구 수 100 미만", "prev": None, "value": n[y], "ratio": None})
g = pd.DataFrame(rows).sort_values(["column", "year"])
g.to_csv(os.path.join(ROOT, "02_패널", "MASTER_연도이상_게이트.csv"), index=False, encoding="utf-8-sig")
by = g.groupby("type").column.nunique().to_dict()
jump = g[g.type == "중앙값 급변"].groupby("year").column.nunique().sort_values(ascending=False).head(12)
lines = ["# 연도 이상 게이트 요약(1998~2026)", f"- 점검 변수 {len(vars_)}개, 표시 {len(g)}건, 유형별 변수 수 {by}",
         "- 중앙값 급변이 몰린 연도(변수 수): " + ", ".join(f"{y}({c})" for y, c in jump.items()),
         "- 해석: 같은 해에 여러 변수가 함께 급변하면 원천 정의·표 교체(예: 검진 2018 판정기준, 2010 판정 중복집계, 대사증후군 2018 수검자수)일 가능성이 큼. 모형 전 해당 연도·변수를 제외하거나 민감도로 둘 것.",
         "\n## 변수별 표시(상위 60)"]
for r in g.head(60).itertuples(): lines.append(f"- {r.column} | {r.year} | {r.type} | {r.prev} → {r.value} (×{r.ratio})")
open(os.path.join(ROOT, "03_문서", "03_연도이상_게이트_요약.md"), "w", encoding="utf-8").write("\n".join(lines))
print(lines[1]); print(lines[2])
