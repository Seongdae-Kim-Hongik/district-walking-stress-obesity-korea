# Reference metadata from PubMed and Crossref.
import os, re, json, time, urllib.request, urllib.parse, glob
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); MS = os.path.dirname(HERE)
man = pd.read_csv(os.path.join(MS, "refs", "manifest", "ALL_refs.csv"), dtype=str).fillna("")
man["pmid"] = man.pmid.str.replace(r"\.0$", "", regex=True); M = man.set_index("key")
src = open(os.path.join(HERE, "content_obes.py"), encoding="utf-8").read()
src = re.sub(r"<!--.*?-->", "", src, flags=re.S); src = re.sub(r"(?m)^\s*#.*$", "", src)
keys = sorted({k.strip().lstrip("@") for g in re.findall(r"\[(@[^\]]+)\]", src) for k in g.split(";")})
keys = [k for k in keys if re.fullmatch(r"(web:)?[a-z][a-z0-9]*", k)]
WEB = {
 "web:kdca2026kchs": dict(type="web", authors="Korea Disease Control and Prevention Agency", title="Korea Community Health Survey: district obesity prevalence (self-reported, table DT_H_OBE_OBE; measured, table DT_117075_MEA_OBE) and health behavior indicators, crude and standardized rates with standard errors, 2008-2025. Korean Statistical Information Service", url="https://kosis.kr", year=2026, accessed="November 10, 2026"),
 "web:mois2026licence": dict(type="web", authors="Ministry of the Interior and Safety", title="Local administrative licensing data: general restaurants and instant food manufacturing and processing businesses, 2011-2024. Public Data Portal", url="https://www.data.go.kr", year=2026, accessed="October 20, 2026"),
 "web:mcst2026facilities": dict(type="web", authors="Ministry of Culture, Sports and Tourism", title="National public sports facility status, 2007-2024", url="https://www.mcst.go.kr", year=2026, accessed="October 18, 2026"),
}
out = {}
miss = [k for k in keys if not k.startswith("web:") and k not in M.index]; assert not miss, miss
pm = [k for k in keys if not k.startswith("web:") and M.loc[k, "pmid"]]
for i in range(0, len(pm), 150):
    ids = ",".join(M.loc[k, "pmid"] for k in pm[i:i + 150])
    d = json.load(urllib.request.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=" + ids, timeout=60))["result"]
    for k in pm[i:i + 150]:
        r = d[M.loc[k, "pmid"]]
        au = [a["name"] for a in r["authors"] if a.get("authtype", "Author") == "Author"] or [a["name"] for a in r["authors"]]
        doi = next((a["value"] for a in r["articleids"] if a["idtype"] == "doi"), "") or M.loc[k, "doi"]
        pages = r.get("pages", "") or (r.get("elocationid", "").replace("pii: ", "").split(" doi")[0] if r.get("elocationid", "").startswith("pii:") else "")
        out[k] = dict(type="journal", authors=au, title=r["title"].rstrip("."), journal=r["source"], year=int(r["pubdate"][:4]), vol=r.get("volume", ""), issue=r.get("issue", ""), pages=pages, doi=doi, pmid=M.loc[k, "pmid"], verified="PubMed esummary " + time.strftime("%Y-%m-%d"))
    time.sleep(0.4)
JABBR = {"Journal of the American Statistical Association": "J Am Stat Assoc", "The American Statistician": "Am Stat", "Political Science Research and Methods": "Polit Sci Res Methods", "Geographical Analysis": "Geogr Anal", "Econometrica": "Econometrica", "Statistical Science": "Stat Sci"}
for k in keys:
    if k.startswith("web:"): out[k] = WEB[k]; continue
    if k in out: continue
    doi = M.loc[k, "doi"]; assert doi, k
    m = json.load(urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(doi), headers={"User-Agent": "research (mailto:kimsd@hongik.ac.kr)"}), timeout=60))["message"]
    au = [f"{a.get('family','')} {''.join(x[0] for x in re.split(r'[ -]', a.get('given','')) if x)}".strip() for a in m.get("author", [])]
    ct = (m.get("container-title") or [""])[0]; jt = JABBR.get(ct, (m.get("short-container-title") or [ct])[0])
    yr = (m.get("published-print") or m.get("published-online") or m.get("issued"))["date-parts"][0][0]
    out[k] = dict(type="journal", authors=au, title=re.sub(r"<[^>]+>", "", (m.get("title") or [""])[0]), journal=jt, year=int(yr), vol=m.get("volume", ""), issue=m.get("issue", ""), pages=m.get("page", "") or m.get("article-number", ""), doi=doi, verified="Crossref " + time.strftime("%Y-%m-%d"))
    time.sleep(0.3)
JN = {"Polit. Sci. Res. Methods": "Polit Sci Res Methods", "Geogr. Anal.": "Geogr Anal", "Stat. Sci.": "Stat Sci", "J R Stat Soc Series B Stat Methodol": "J R Stat Soc Series B Stat Methodol",
      "Journal of the Royal Statistical Society Series B: Statistical Methodology": "J R Stat Soc Series B Stat Methodol"}
def fullp(p):
    m = re.match(r"^(\d+)-(\d+)$", str(p))
    if not m: return p
    a, b = m.groups(); return p if len(b) >= len(a) else f"{a}-{a[:len(a)-len(b)]+b}"
FIXP = {"choi2020walkability": "2060", "kim2017seoul": "898", "netherlands2020clustering": "839", "walkability2019adults": "3135"}
AHEAD = set()
for k, v in out.items():
    if v.get("type") != "journal": continue
    if k in FIXP: v["pages"] = FIXP[k]
    if k == "fay1979estimates": v["issue"] = "366"
    if k in AHEAD: v["pages"] = ""; v["vol"] = ""; v["ahead"] = True
    v["journal"] = JN.get(v["journal"], v["journal"]).replace(".", "")
    v["pages"] = fullp(v.get("pages", "")) if k not in AHEAD else ""
    if v["title"].isupper(): v["title"] = v["title"].capitalize()
    v["authors"] = [a.title() if a.isupper() else a for a in v["authors"]]
SENT = {
 "rue2009approximate": "Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations",
 "knorrheld2000bayesian": "Bayesian modelling of inseparable space-time variation in disease risk",
 "kang2015korea": "Korea Community Health Survey data profiles",
 "simpson2017penalising": "Penalising model component complexity: a principled, practical approach to constructing priors",
 "fay1979estimates": "Estimates of income for small places: an application of James-Stein procedures to census data",
 "hodges2010adding": "Adding spatially-correlated errors can mess up the fixed effect you love",
 "mercer2015space": "Space-time smoothing of complex survey data: small area estimation for child mortality",
 "riebler2016intuitive": "An intuitive Bayesian spatial model for disease mapping that accounts for scaling",
 "swinburn2019global": "The global syndemic of obesity, undernutrition, and climate change: the Lancet Commission report",
 "walkability2019adults": "Walkability, overweight, and obesity in adults: a systematic review of observational studies",
 "kim2017seoul": "The relationship between obesity and urban environment in Seoul",
 "kim2026obesity": "2025 obesity fact sheet for Korea: prevalence of obesity, abdominal obesity from 2014 to 2023 and prevalence of chronic diseases according to obesity status",
 "steegen2016increasing": "Increasing transparency through a multiverse analysis",
 "traversy2015alcohol": "Alcohol consumption and obesity: an update",
 "choi2020walkability": "Do the walkability and urban leisure amenities of neighborhoods affect the body mass index of individuals? Based on a case study in Seoul, South Korea",
 "namgung2019built": "The role of built environment on health of older adults in Korea: obesity and gender differences",
}
for k, t in SENT.items():
    if k in out: out[k]["title"] = t
for k, v in out.items():
    if v.get("type") == "journal":
        v["title"] = re.sub(r":\s+([A-Z][a-z]*)\b", lambda m: ": " + (m.group(1) if m.group(1) in ("Korea", "Korean", "South", "Bayesian", "Spatial", "E", "GOLD") else m.group(1).lower()), v["title"])
AUTHOR_FIX = {"Paulo Dos Anjos Souza Barbosa J": "Barbosa JPAS", "Henrique Guerra P": "Guerra PH", "de Oliveira Santos C": "Santos CO", "de Oliveira Barbosa Nunes AP": "Nunes APOB", "Antonio Florindo A": "Florindo AA", "Cheol Seong S": "Seong SC", "Heon Park J": "Park JH", "Hyon Bang J": "Bang JH", "Ae Shin S": "Shin SA"}
for _k, _v in out.items():
    if isinstance(_v.get("authors"), list): _v["authors"] = [AUTHOR_FIX.get(a, a) for a in _v["authors"]]
assert not any("Cheol Seong" in a for v in out.values() if isinstance(v.get("authors"), list) for a in v["authors"])
json.dump(out, open(os.path.join(HERE, "refs_obes.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
import re as _re
for _k, _v in out.items():
    if isinstance(_v.get("authors"), list): assert not [a for a in _v["authors"] if a not in ("Tchetgen Tchetgen E", "De Schutter O", "Paulo Dos Anjos Souza Barbosa J") and _re.fullmatch(r"[A-Z][a-z]+ [A-Z][a-z]+ [A-Z]", a)], (_k, _v["authors"])
print(len(keys), "keys;", sum(v["type"] == "journal" for v in out.values()), "journal,", sum(v["type"] == "web" for v in out.values()), "web")
