#!/usr/bin/env python3
# Standard sensitivity sets s577, s2 and s8 written after the main results (marked post hoc): copies of the main configuration with only the changed keys. Environment: E2_RUN (run root with 02_패널/ and 04_분석/OBES/ model inputs), E2_OUT (results folder, default 04_분석/OBES/results_extra2).
import os, json, hashlib, datetime
OBD = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.environ.get("E2_OUT", os.path.join(OBD, "results_extra2")); RUN = os.environ["E2_RUN"]
base = json.load(open(os.path.join(RUN, "04_분석/OBES/sets_OBES_std3.json")))
cfg = json.load(open(os.path.join(RUN, "04_분석/OBES/config_OBES_std3.json")))
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
common = {"posthoc": True, "posthoc_note": f"Written after the main results were seen ({now}): standard sensitivity sets s577, s2 and s8 adapted to this outcome. Not fixed before fitting. The main sets file sets_OBES_std3.json is unchanged.",
          "base_sets_file": "sets_OBES_std3.json", "base_sets_sha256": hashlib.sha256(open(os.path.join(RUN, "04_분석/OBES/sets_OBES_std3.json"), "rb").read()).hexdigest(), "main_set": cfg["rules"]["main_set"], "forced": cfg["forced"]}
V = {"s577": dict(sens_note="그래프만 queen 인접(셰이프 재생성 queen + 섬 최근접 연결, 577 간선); 공변량·사전분포·우도 그대로", graph="queen+island nearest (577 edges)", fit="fits/fit_residspec_S13_queen_s577_*.csv (= S13 queen 행과 같은 모형)"),
     "s2": dict(sens_note="80세 이상 비율을 강제통제에서 뺌(성·연령 직접표준화 결과와 겹치는 축)", forced=[f for f in cfg["forced"] if f != "pop80_share"], fit="fits/fit_s2_no_pop80_forced_*.csv"),
     "s8": dict(sens_note="시리즈 원정의 = 연령만 표준화 기대건수(E_std_age). 이 결과변수는 KCHS 직접표준화 유병률(성·연령, 2005 추계인구)이며 원표에 연령만 표준화율이 없음(기준 SQLite region_year 비만 열 = SR·CR·SE·N뿐). 대체 = 조율(CR) + 65세 이상 비율 통제 FH(S12 crude_rate_FH_plus_old65)",
                applicable=False, substitute="results_std3/S12_likelihood.csv likelihood=crude_rate_FH_plus_old65", fit="(대체 적합 = S12)")}
rows = []
for k, v in V.items():
    j = dict(base); j.update(common); j.update(v); j["version_tag"] = k
    p = os.path.join(OUT, "sets", f"sets_OBES_{k}.json"); json.dump(j, open(p, "w"), ensure_ascii=False, indent=1)
    rows.append(f"{k} sha256 {hashlib.sha256(open(p, 'rb').read()).hexdigest()}")
open(os.path.join(OUT, "sets", "sets_OBES_s577_s2_s8_sha256_posthoc.txt"), "w").write("\n".join(rows) + "\n"); print("\n".join(rows))
