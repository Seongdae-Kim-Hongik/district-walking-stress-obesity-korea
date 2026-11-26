# Reference metadata from a Zotero collection, formatted with a PubMed/Crossref cache (writes refs_obes.json).
# Caches: refs_cache_zotero.json (--refresh-zotero, local Zotero API; ZOTERO_USER, ZOTERO_COLLECTION) and refs_cache_pubmed.json (--refresh-pubmed, runs build_refs_pubmed.py).
# REFS_SOURCE = manuscript source files with the citation keys, separated by the path separator.
import os, re, sys, json, subprocess, unicodedata, urllib.request
HERE = os.path.dirname(os.path.abspath(__file__))
ZC, PC = os.path.join(HERE, "refs_cache_zotero.json"), os.path.join(HERE, "refs_cache_pubmed.json")
LIB, COL = os.environ.get("ZOTERO_USER", "0"), os.environ.get("ZOTERO_COLLECTION", "")
if "--refresh-zotero" in sys.argv:
    out = []
    for st in range(0, 2000, 100):
        d = json.load(urllib.request.urlopen(f"http://localhost:23119/api/users/{LIB}/collections/{COL}/items/top?limit=100&start={st}&format=json", timeout=60)); out += d
        if len(d) < 100: break
    json.dump(out, open(ZC, "w", encoding="utf-8"), ensure_ascii=False)
if "--refresh-pubmed" in sys.argv:
    bak = open(os.path.join(HERE, "refs_obes.json"), encoding="utf-8").read()
    subprocess.run([sys.executable, os.path.join(HERE, "build_refs_pubmed.py")], check=True)
    P = json.load(open(os.path.join(HERE, "refs_obes.json"), encoding="utf-8"))
    json.dump({k: v for k, v in P.items() if v.get("type") == "journal"}, open(PC, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    open(os.path.join(HERE, "refs_obes.json"), "w", encoding="utf-8").write(bak)
Z = json.load(open(ZC, encoding="utf-8")); P = json.load(open(PC, encoding="utf-8"))
src = "".join(open(p, encoding="utf-8").read() for p in os.environ["REFS_SOURCE"].split(os.pathsep))
src = re.sub(r"(?m)^\s*#.*$", "", src)
keys = sorted({k.strip().lstrip("@") for g in re.findall(r"\[(@[^\]]+)\]", src) for k in g.split(";")})
keys = [k for k in keys if re.fullmatch(r"(web:)?[a-z][a-z0-9]*", k)]
byTag = {}
for it in Z:
    for t in it["data"].get("tags", []): byTag.setdefault(t["tag"], []).append(it)
miss = [k for k in keys if len(byTag.get(k, [])) != 1]; assert not miss, ("Zotero 태그 매핑 실패", miss)
ZTAG = r"</?(?:i|b|sup|scp)>"
def norm(s): s = unicodedata.normalize("NFKC", re.sub(ZTAG, "", s)).replace("‐", "-").replace("‑", "-"); return re.sub(r"[^0-9a-z]+", "", s.casefold())
PART = {"de", "da", "do", "dos", "das", "van", "von", "der", "den", "del", "la", "le"}
def initials(first): return "".join(x[0].upper() for x in re.split(r"[ .\-\u2010]+", first) if x and x not in PART)
def sname(n): n = n.replace("\u2010", "-"); return n.title() if n.isupper() else n
def fullp(p):
    m = re.match(r"^(\d+)-(\d+)$", str(p));
    if not m: return str(p)
    a, b = m.groups(); return p if len(b) >= len(a) else f"{a}-{a[:len(a)-len(b)] + b}"
out, diff = {}, []
for k in keys:
    it = byTag[k][0]; d = it["data"]; zk = it["key"]
    if d["itemType"] in ("webpage", "computerProgram", "report", "dataset"):
        ex = dict(re.findall(r"(?m)^(\w+):\s*(.+)$", d.get("extra", "")))
        acc = d.get("accessDate", "")[:10] or (lambda m: f"{m[0]}-{int(m[1]):02d}-{int(m[2]):02d}" if m else "")(re.findall(r"\d+", ex.get("accessed", "")))
        assert acc, (k, "접근일 없음")
        import calendar; y, mo, dd = acc.split("-")
        title = re.sub(ZTAG, "", d["title"])
        if d["itemType"] == "computerProgram":
            title = re.sub(r"^R: A ", "R: a ", title) + ". " + ex["company"] + (f", {d['place']}" if d.get("place") else "")
        out[k] = dict(type="web", authors=d["creators"][0].get("name") or d["creators"][0]["lastName"], title=title, url=d["url"], year=int(re.findall(r"\d{4}", d.get("date", ""))[0]),
                      accessed=f"{calendar.month_name[int(mo)]} {int(dd)}, {y}", zotero=zk, verified="Zotero SPPE_OBES_2026 " + d.get("dateModified", "")[:10])
        continue
    au = [f"{sname(c.get('lastName', c.get('name', '')))} {initials(c.get('firstName', ''))}".strip() for c in d["creators"] if c.get("creatorType", "author") == "author"]
    p = P.get(k, {}); zt = re.sub(ZTAG, "", d["title"]).replace("‐", "-").rstrip(".")
    title = p["title"] if p and norm(p["title"]) == norm(zt) else zt
    yr = int(re.findall(r"\d{4}", d.get("date", ""))[0])
    zau = au
    if p.get("authors") and [a.split(" ")[0].casefold() for a in p["authors"]] == [a.split(" ")[0].casefold() for a in au]: au = p["authors"]
    e = dict(type="journal", authors=au, title=title, journal=p.get("journal") or d.get("journalAbbreviation") or d.get("publicationTitle"), year=yr,
             vol=p.get("vol", d.get("volume", "")), issue=p.get("issue", d.get("issue", "")), pages=p.get("pages", fullp(d.get("pages", "").replace("–", "-"))), doi=(p["doi"] if p.get("doi") and p["doi"].lower() == d.get("DOI", "").lower() else d.get("DOI", "") or p.get("doi", "")), zotero=zk,
             verified="Zotero SPPE_OBES_2026 " + d.get("dateModified", "")[:10] + " (format: " + p.get("verified", "no PubMed/Crossref cache") + ")")
    if p.get("pmid"): e["pmid"] = p["pmid"]
    zv = dict(authors=zau, title=title, year=yr, vol=d.get("volume", ""), issue=d.get("issue", ""), pages=fullp(d.get("pages", "").replace("–", "-")), doi=d.get("DOI", ""))
    for f in ("authors", "title", "year", "vol", "issue", "pages", "doi"):
        a, b = zv[f], p.get(f)
        if f == "title" and b is not None and norm(a) == norm(b): continue
        if f == "doi" and b and a.lower() == b.lower(): continue
        if f == "pages" and b and str(a) == str(b): continue
        if b is not None and str(a) != str(b): diff.append((k, zk, f, b, a))
    out[k] = e
json.dump(out, open(os.path.join(HERE, "refs_obes.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
with open(os.path.join(HERE, "refs_zotero_diff.md"), "w", encoding="utf-8") as f:
    f.write(f"# Zotero item metadata vs PubMed/Crossref cache: {len(out)} keys, {len(diff)} differing fields\n\nRule: key mapping, title text, year, DOI and web sources (URL, access date) from Zotero; volume, issue, pages and initials (same surnames) from PubMed/Crossref; title case and journal abbreviation in PubMed/Crossref format.\n\n| key | Zotero item | field | PubMed/Crossref (adopted) | Zotero value |\n|---|---|---|---|---|\n")
    for r in diff: f.write("| " + " | ".join(str(x).replace("|", "/") for x in r) + " |\n")
print(len(keys), "keys from Zotero;", sum(v["type"] == "journal" for v in out.values()), "journal,", sum(v["type"] == "web" for v in out.values()), "web; PubMed/Crossref 차이", len(diff))
