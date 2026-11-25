# KOSIS open API: item metadata of table DT_1B04005N (resident registration population by district). Key from KOSIS_API_KEY. Usage: python3 kosis_fetch_meta.py <out folder>
import json, urllib.request, urllib.parse, os, sys
KEY = os.environ["KOSIS_API_KEY"]
for typ in ["ITM"]:
    q = dict(method="getMeta", type=typ, apiKey=KEY, format="json", jsonVD="Y", orgId="101", tblId="DT_1B04005N")
    d = json.loads(urllib.request.urlopen("https://kosis.kr/openapi/statisticsData.do?" + urllib.parse.urlencode(q), timeout=120).read().decode())
    json.dump(d, open(os.path.join(sys.argv[1], f"meta_{typ}.json"), "w"), ensure_ascii=False)
    items = [(x.get("OBJ_ID"), x.get("ITM_ID"), x.get("ITM_NM")) for x in d]
    print(len(items)); print([i for i in items if i[0] != "A"][:10] if items else d)
    print(sorted(set(i[0] for i in items)))
    print([i for i in items if i[0] not in ("A",) and str(i[1]) in ("0", "T2", "T3", "T4", "T1")][:10])
    print([i for i in items if i[0]=="ITEM"][:10])
