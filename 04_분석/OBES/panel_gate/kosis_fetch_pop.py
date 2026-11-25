# KOSIS open API: total resident registration population by district, December 2011-2024 (DT_1B04005N, item T2), for the denominator check. Key from KOSIS_API_KEY. Usage: python3 kosis_fetch_pop.py <out folder>
import json, urllib.request, urllib.parse, os, sys, time, hashlib
KEY = os.environ["KOSIS_API_KEY"]; OUT = sys.argv[1]
for y in range(2011, 2025):
    fp = os.path.join(OUT, f"DT_1B04005N_T2_{y}12.json")
    if os.path.exists(fp): continue
    q = dict(method="getList", apiKey=KEY, itmId="T2+", format="json", jsonVD="Y", prdSe="M", startPrdDe=f"{y}12", endPrdDe=f"{y}12", orgId="101", tblId="DT_1B04005N", objL1="ALL", objL2="0+")
    for k in range(3):
        try:
            raw = urllib.request.urlopen("https://kosis.kr/openapi/Param/statisticsParameterData.do?" + urllib.parse.urlencode(q), timeout=300).read(); d = json.loads(raw.decode()); break
        except Exception as e: print("retry", y, e); time.sleep(5)
    open(fp, "wb").write(raw); print(y, len(d) if isinstance(d, list) else d, hashlib.sha256(raw).hexdigest()[:16])
