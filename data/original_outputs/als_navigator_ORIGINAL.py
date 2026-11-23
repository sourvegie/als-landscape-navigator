#!/usr/bin/env python3
"""
ALS Landscape Navigator - minimal demo
======================================

Turns public data into a decision-ready view of the ALS therapeutic landscape:

  * ClinicalTrials.gov API v2  -> every ALS study (phase, status, sponsor, sites, dates)
  * Open Targets Platform      -> genetic / overall evidence for ALS-associated targets
  * PubMed E-utilities         -> literature volume per target (all-time and last 5 years)

Each drug intervention is classified by modality (ASO, AAV, antibody, ...) and by
mechanism (SOD1, C9orf72, TDP-43 biology, neuroinflammation, ...) with transparent,
auditable rules. The output is ONE self-contained HTML dashboard plus CSV tables.

Standard library only - no pip install needed.

Usage
-----
  python3 als_navigator.py                   # fetch (or reuse cache) and build the dashboard
  python3 als_navigator.py --refresh         # force a fresh download
  python3 als_navigator.py --evaluate        # score your manual labels in data/label_sample.csv

Outputs
-------
  output/als_landscape.html          the dashboard (open in any browser; host on GitHub Pages)
  data/trials.csv                    one row per trial, with classifications
  data/targets.csv                   Open Targets + PubMed evidence per target
  data/unclassified_interventions.csv  interventions the rules could not place (fill in -> overrides)
  data/label_sample.csv              50 random drug trials for manual validation (never overwritten)
  data/raw/*.json                    cached raw API responses (reproducibility)
"""

import argparse
import collections
import csv
import datetime as dt
import html
import json
import os
import random
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

try:  # macOS python.org builds often lack root certificates; certifi (bundled with conda) fixes that
    import ssl
    import certifi
    _SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except Exception:
    _SSL_CTX = None

# --------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------
CONDITION = "amyotrophic lateral sclerosis"
CTGOV_URL = "https://clinicaltrials.gov/api/v2/studies"
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
USER_AGENT = "ALS-Landscape-Navigator/0.1 (research demo; python urllib)"

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
OUT_DIR = os.path.join(HERE, "output")

N_OT_TARGETS = 30          # top Open Targets associations to pull
N_PUBMED_TARGETS = 20      # how many of those get PubMed counts (rate-limit friendly)
LABEL_SAMPLE_SIZE = 50
RANDOM_SEED = 42

ACTIVE_STATUSES = {"RECRUITING", "NOT_YET_RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION"}
STOPPED_STATUSES = {"TERMINATED", "WITHDRAWN", "SUSPENDED"}
DRUG_TYPES = {"DRUG", "BIOLOGICAL", "GENETIC", "COMBINATION_PRODUCT", "DIETARY_SUPPLEMENT"}
SOC_RE = re.compile(r"^\s*(riluzole|edaravone|standard of care|best supportive care)\b", re.I)
PLACEBO_RE = re.compile(r"\b(placebo|sham|vehicle|matching control|saline control)\b", re.I)

MODALITIES = ["ASO", "siRNA / RNAi", "AAV gene therapy", "Vectorized antibody", "Antibody",
              "Cell therapy", "Other biologic", "Small molecule", "Supplement"]
MECHANISMS = ["SOD1", "C9orf72", "FUS", "ATXN2", "TDP-43 biology",
              "Neuroinflammation & immune", "Oxidative stress", "Excitotoxicity",
              "Mitochondria & cell stress", "Integrated stress response", "Proteostasis & autophagy",
              "Neurotrophic & cell-based", "Muscle & NMJ", "Other neuroprotective",
              "Symptomatic", "Unclassified"]
# mechanism bucket -> gene symbols used to join with Open Targets / PubMed
MECH_TO_GENES = {"SOD1": ["SOD1"], "C9orf72": ["C9orf72", "C9ORF72"], "FUS": ["FUS"],
                 "ATXN2": ["ATXN2"], "TDP-43 biology": ["TARDBP", "STMN2"]}
PUBMED_ALIASES = {"TARDBP": '(TARDBP[tiab] OR "TDP-43"[tiab] OR TDP43[tiab])',
                  "C9orf72": '(C9orf72[tiab] OR "C9ORF72"[tiab])'}
ASIA = {"Japan", "China", "Korea, Republic of", "South Korea", "Taiwan", "Singapore", "Hong Kong",
        "India", "Thailand", "Malaysia", "Vietnam", "Viet Nam", "Philippines", "Indonesia"}

# --------------------------------------------------------------------------------------
# Classification rules (ordered: first match wins for each axis).
# Each rule: (rule_id, regex, mechanism or None, modality or None)
# Named-asset rules come first (high confidence), keyword rules after (lower confidence).
# --------------------------------------------------------------------------------------
NAMED_RULES = [
    ("tofersen", r"tofersen|BIIB067|IONIS[- ]?SOD1", "SOD1", "ASO"),
    ("ulefnersen", r"ulefnersen|jacifusen|ION[- ]?363", "FUS", "ASO"),
    ("BIIB078", r"BIIB078|IONIS[- ]?C9", "C9orf72", "ASO"),
    ("WVE-004", r"WVE[- ]?004", "C9orf72", "ASO"),
    ("BIIB105", r"BIIB105|ION[- ]?541", "ATXN2", "ASO"),
    ("AMT-162", r"AMT[- ]?162", "SOD1", "AAV gene therapy"),
    ("AP-101", r"\bAP[- ]?101\b", "SOD1", "Antibody"),
    ("VTx-002", r"VTx[- ]?002", "TDP-43 biology", "Vectorized antibody"),
    ("QRL-201", r"QRL[- ]?201", "TDP-43 biology", "ASO"),
    ("PIKfyve", r"apilimod|\bAS[- ]?202\b|PIKfyve", "TDP-43 biology", "Small molecule"),
    ("edaravone", r"edaravone|MT[- ]?1186|radicava", "Oxidative stress", "Small molecule"),
    ("riluzole", r"riluzole|troriluzole", "Excitotoxicity", "Small molecule"),
    ("AMX0035", r"AMX0035|AMX[- ]0035|phenylbutyrate|taurursodiol|TUDCA|tauroursodeoxychol",
     "Mitochondria & cell stress", "Small molecule"),
    ("masitinib", r"masitinib", "Neuroinflammation & immune", "Small molecule"),
    ("ibudilast", r"ibudilast|MN[- ]?166", "Neuroinflammation & immune", "Small molecule"),
    ("verdiperstat", r"verdiperstat|BHV[- ]?3241", "Neuroinflammation & immune", "Small molecule"),
    ("zilucoplan", r"zilucoplan", "Neuroinflammation & immune", "Other biologic"),
    ("complement-mab", r"ravulizumab|eculizumab", "Neuroinflammation & immune", "Antibody"),
    ("Treg", r"regulatory T[- ]?cell|\bTregs?\b", "Neuroinflammation & immune", "Cell therapy"),
    ("IL-2", r"interleukin[- ]?2\b|\bIL[- ]?2\b|aldesleukin", "Neuroinflammation & immune", "Other biologic"),
    ("pridopidine", r"pridopidine", "Other neuroprotective", "Small molecule"),
    ("ISR", r"DNL[- ]?343|ABBV[- ]?CLS[- ]?7262|eIF2B", "Integrated stress response", "Small molecule"),
    ("CNM-Au8", r"CNM[- ]?Au8", "Mitochondria & cell stress", "Other biologic"),
    ("troponin", r"reldesemtiv|tirasemtiv|CK[- ]?2127107|troponin activator", "Muscle & NMJ", "Small molecule"),
    ("proteostasis-sm", r"trehalose|rapamycin|sirolimus|arimoclomol|colchicine", "Proteostasis & autophagy",
     "Small molecule"),
    ("NurOwn", r"NurOwn|MSC[- ]?NTF", "Neurotrophic & cell-based", "Cell therapy"),
    ("methylcobalamin", r"methyl[- ]?cobalamin|mecobalamin", "Other neuroprotective", "Small molecule"),
    ("symptomatic", r"dextromethorphan|quinidine|nuedexta|baclofen|mexiletine|botulinum|glycopyrr|"
                    r"scopolamine|atropine|morphine|cannabi|nabiximols", "Symptomatic", None),
]
KEYWORD_RULES = [
    # modality keywords
    ("kw-ASO", r"antisense|oligonucleotide|\b\w+nersen\b", None, "ASO"),
    ("kw-siRNA", r"\bsiRNA\b|\bRNAi\b|\b\w+siran\b", None, "siRNA / RNAi"),
    ("kw-AAV", r"\bAAV\d*\b|adeno[- ]associated|gene therapy|gene transfer", None, "AAV gene therapy"),
    ("kw-cell", r"stem cell|mesenchymal|\bMSCs?\b|cell therapy|autologous .*cells|bone marrow", None, "Cell therapy"),
    ("kw-mab", r"monoclonal|antibod|\b\w+mab\b|immunoglobulin|\bIVIG\b", None, "Antibody"),
    # mechanism keywords
    ("kw-SOD1", r"\bSOD[- ]?1\b|superoxide dismutase", "SOD1", None),
    ("kw-C9", r"C9orf72|C9ORF72|\bC9\b", "C9orf72", None),
    ("kw-FUS", r"\bFUS\b", "FUS", None),
    ("kw-ATXN2", r"ATXN2|ataxin", "ATXN2", None),
    ("kw-TDP43", r"TDP[- ]?43|TARDBP|STMN2|stathmin", "TDP-43 biology", None),
    ("kw-inflam", r"inflamm|microglia|complement|immun|cytokine|\bTNF\b", "Neuroinflammation & immune", None),
    ("kw-oxid", r"oxidat|antioxid|free radical|glutathione|N-acetylcysteine", "Oxidative stress", None),
    ("kw-glut", r"glutamat|excitotox|AMPA|NMDA|perampanel|memantine", "Excitotoxicity", None),
    ("kw-mito", r"mitochond|bioenerget|coenzyme ?Q|CoQ10|creatine|ketogenic", "Mitochondria & cell stress", None),
    ("kw-ISR", r"integrated stress|\bISR\b", "Integrated stress response", None),
    ("kw-proteo", r"autophag|proteasom|heat shock|\bHSP\d*\b|chaperone|aggregat", "Proteostasis & autophagy", None),
    ("kw-neurotroph", r"neurotroph|GDNF|BDNF|VEGF|IGF[- ]?1|growth factor", "Neurotrophic & cell-based", None),
    ("kw-symptom", r"spasticity|cramp|sialorrh|saliva|pain|pseudobulbar|dyspn", "Symptomatic", None),
]
TYPE_FALLBACK_MODALITY = {"DRUG": "Small molecule", "BIOLOGICAL": "Other biologic", "GENETIC": "AAV gene therapy",
                          "COMBINATION_PRODUCT": "Small molecule", "DIETARY_SUPPLEMENT": "Supplement"}
SOC_RULE_IDS = {"riluzole", "edaravone"}
_COMPILED = [(rid, re.compile(rx, re.I), mech, mod) for rid, rx, mech, mod in NAMED_RULES + KEYWORD_RULES]

# Known-answer checks: (NCT ID, expected mechanism, expected modality, label)
KNOWN_TRIALS = [
    ("NCT02623699", "SOD1", "ASO", "VALOR - tofersen, SOD1-ALS"),
    ("NCT04856982", "SOD1", "ASO", "ATLAS - tofersen, presymptomatic SOD1 carriers"),
]
# Known-answer checks on intervention names: (regex, expected mechanism, expected modality)
KNOWN_ASSETS = [
    (r"tofersen", "SOD1", "ASO"),
    (r"ulefnersen|jacifusen|ION[- ]?363", "FUS", "ASO"),
    (r"VTx[- ]?002", "TDP-43 biology", "Vectorized antibody"),
    (r"edaravone", "Oxidative stress", "Small molecule"),
    (r"riluzole", "Excitotoxicity", "Small molecule"),
]


# --------------------------------------------------------------------------------------
# HTTP helpers (stdlib only, polite retries)
# --------------------------------------------------------------------------------------
def _request(url, data=None, headers=None, retries=4, timeout=60):
    hdrs = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    hdrs.update(headers or {})
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers=hdrs)
            with urllib.request.urlopen(req, timeout=timeout, context=_SSL_CTX) as r:
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (400, 404):
                raise
            time.sleep(1.5 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            last = e
            time.sleep(1.5 * (attempt + 1))
    hint = ""
    if "CERTIFICATE_VERIFY_FAILED" in str(last):
        hint = ("\n  hint: run 'pip install certifi', or on python.org Python run "
                "'/Applications/Python 3.x/Install Certificates.command', or use your conda Python.")
    raise RuntimeError(f"Request failed after {retries} attempts: {url}\n  last error: {last}{hint}")


def get_json(url, params=None):
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    return _request(url)


def post_json(url, payload):
    return _request(url, data=json.dumps(payload).encode("utf-8"), headers={"Content-Type": "application/json"})


def now_iso():
    return dt.datetime.now().replace(microsecond=0).isoformat()


# --------------------------------------------------------------------------------------
# Fetchers (each caches its raw response in data/raw/)
# --------------------------------------------------------------------------------------
def fetch_ctgov():
    print(f"[ClinicalTrials.gov] querying condition = '{CONDITION}' ...")
    params = {"query.cond": CONDITION, "pageSize": 1000, "countTotal": "true", "format": "json"}
    studies, total, token, page = [], None, None, 0
    while True:
        p = dict(params)
        if token:
            p["pageToken"] = token
        resp = get_json(CTGOV_URL, p)
        if total is None:
            total = resp.get("totalCount")
        studies.extend(resp.get("studies", []))
        page += 1
        print(f"   page {page}: {len(studies)} / {total} studies")
        token = resp.get("nextPageToken")
        if not token:
            break
        time.sleep(0.3)
    return {"fetched_at": now_iso(), "query": CONDITION, "totalCount": total, "studies": studies}


OT_SEARCH_Q = """
query Search($q: String!) {
  search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 5}) {
    hits { id name entity }
  }
}"""
OT_ASSOC_Q = """
query Assoc($id: String!, $size: Int!) {
  disease(efoId: $id) {
    id name
    associatedTargets(page: {index: 0, size: $size}) {
      count
      rows { score target { id approvedSymbol approvedName } %s }
    }
  }
}"""


def fetch_opentargets():
    print("[Open Targets] resolving disease ID ...")
    hits = post_json(OT_URL, {"query": OT_SEARCH_Q, "variables": {"q": CONDITION}})
    hits = (((hits.get("data") or {}).get("search") or {}).get("hits")) or []
    disease = next((h for h in hits if h.get("name", "").lower() == CONDITION), hits[0] if hits else None)
    if not disease:
        raise RuntimeError("Open Targets: disease not found")
    print(f"   -> {disease['name']} ({disease['id']}); pulling top {N_OT_TARGETS} targets ...")
    try:
        q = OT_ASSOC_Q % "datatypeScores { id score }"
        resp = post_json(OT_URL, {"query": q, "variables": {"id": disease["id"], "size": N_OT_TARGETS}})
        if resp.get("errors"):
            raise RuntimeError(resp["errors"])
    except Exception as e:  # schema drift fallback: overall score only
        print(f"   datatypeScores unavailable ({str(e)[:80]}); falling back to overall score")
        q = OT_ASSOC_Q % ""
        resp = post_json(OT_URL, {"query": q, "variables": {"id": disease["id"], "size": N_OT_TARGETS}})
    d = (resp.get("data") or {}).get("disease") or {}
    rows = (d.get("associatedTargets") or {}).get("rows") or []
    return {"fetched_at": now_iso(), "disease": {"id": disease["id"], "name": disease["name"]},
            "count": (d.get("associatedTargets") or {}).get("count"), "rows": rows}


def pubmed_count(term, mindate=None, maxdate=None):
    p = {"db": "pubmed", "term": term, "retmode": "json", "rettype": "count"}
    if os.environ.get("NCBI_API_KEY"):
        p["api_key"] = os.environ["NCBI_API_KEY"]
    if mindate:
        p.update({"datetype": "pdat", "mindate": str(mindate), "maxdate": str(maxdate)})
    r = get_json(EUTILS_URL, p)
    time.sleep(0.12 if os.environ.get("NCBI_API_KEY") else 0.4)  # NCBI: 3 req/s without key
    return int(r["esearchresult"]["count"])


def fetch_pubmed(symbols):
    yr = dt.date.today().year
    base = f'"{CONDITION}"[tiab]'
    print(f"[PubMed] counting papers for {len(symbols)} targets ...")
    out = {"fetched_at": now_iso(), "als_total": pubmed_count(base), "window": [yr - 4, yr], "counts": {}}
    for s in symbols:
        t = f"{base} AND {PUBMED_ALIASES.get(s, s + '[tiab]')}"
        out["counts"][s] = {"total": pubmed_count(t), "recent": pubmed_count(t, yr - 4, yr), "term": t}
    return out


def load_or_fetch(name, fetch_fn, refresh):
    path = os.path.join(RAW_DIR, name)
    if os.path.exists(path) and not refresh:
        with open(path) as f:
            print(f"[cache] using {os.path.relpath(path, HERE)}")
            return json.load(f)
    data = fetch_fn()
    os.makedirs(RAW_DIR, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f)
    return data


# --------------------------------------------------------------------------------------
# Parsing & classification
# --------------------------------------------------------------------------------------
def parse_date(s):
    if not s:
        return None
    for fmt in ("%Y-%m-%d", "%Y-%m"):
        try:
            return dt.datetime.strptime(s, fmt).date()
        except ValueError:
            pass
    return None


PHASE_RANK = {"EARLY_PHASE1": 0.5, "PHASE1": 1, "PHASE2": 2, "PHASE3": 3, "PHASE4": 4}


def phase_label(phases):
    ph = [p for p in (phases or []) if p in PHASE_RANK]
    if not ph:
        return "N/A", 0
    ph.sort(key=PHASE_RANK.get)
    nums = [("Early 1" if p == "EARLY_PHASE1" else p.replace("PHASE", "")) for p in ph]
    return "Ph " + "/".join(nums), PHASE_RANK[ph[-1]]


def load_overrides():
    """data/overrides.csv: intervention_name, mechanism, modality - your manual fixes, applied first."""
    path = os.path.join(DATA_DIR, "overrides.csv")
    ov = {}
    if os.path.exists(path):
        with open(path, newline="") as f:
            for r in csv.DictReader(f):
                name = (r.get("intervention_name") or "").strip().lower()
                if name and ((r.get("mechanism") or "").strip() or (r.get("modality") or "").strip()):
                    ov[name] = ((r.get("mechanism") or "").strip() or None, (r.get("modality") or "").strip() or None)
    return ov


def classify_intervention(iv, overrides):
    name = iv.get("name") or ""
    name_text = " | ".join([name, " ".join(iv.get("otherNames") or [])])
    text = " | ".join([name_text, iv.get("description") or ""])
    mech = mod = None
    mech_rule = mod_rule = None
    if name.strip().lower() in overrides:
        mech, mod = overrides[name.strip().lower()]
        mech_rule = "override" if mech else None
        mod_rule = "override" if mod else None
    for rid, rx, r_mech, r_mod in _COMPILED:
        if mech and mod:
            break
        # Standard-of-care drugs (riluzole, edaravone) are often named in the description only as
        # background therapy or a double-dummy comparator, so match them on the name alone.
        hay = name_text if rid in SOC_RULE_IDS else text
        if rx.search(hay):
            if r_mech and not mech:
                mech, mech_rule = r_mech, rid
            if r_mod and not mod:
                mod, mod_rule = r_mod, rid
    if not mod:
        mod = TYPE_FALLBACK_MODALITY.get(iv.get("type"), "Small molecule")
        mod_rule = "type-fallback"
    if not mech:
        mech, mech_rule = "Unclassified", None
    return {"name": name, "type": iv.get("type"), "mechanism": mech, "modality": mod,
            "mech_rule": mech_rule, "mod_rule": mod_rule}


def normalise_trials(ct, overrides):
    trials = []
    for s in ct["studies"]:
        ps = s.get("protocolSection", {})
        idm, stm = ps.get("identificationModule", {}), ps.get("statusModule", {})
        dsm, spm = ps.get("designModule", {}), ps.get("sponsorCollaboratorsModule", {})
        aim, clm = ps.get("armsInterventionsModule", {}), ps.get("contactsLocationsModule", {})
        if dsm.get("studyType") != "INTERVENTIONAL":
            continue
        ivs_all = aim.get("interventions") or []
        ivs = [iv for iv in ivs_all if not PLACEBO_RE.search(iv.get("name") or "")]
        drug_ivs = [iv for iv in ivs if iv.get("type") in DRUG_TYPES]
        classified = [classify_intervention(iv, overrides) for iv in drug_ivs]
        # Standard-of-care drugs listed alongside an investigational agent are background therapy,
        # not the mechanism being tested - drop them from the trial's classification.
        background = []
        if len(classified) > 1:
            inv = [c for c in classified if not SOC_RE.search(c["name"])]
            if inv:
                background = [c["name"] for c in classified if SOC_RE.search(c["name"])]
                classified = inv
        plab, prank = phase_label(dsm.get("phases"))
        countries = sorted({loc.get("country") for loc in (clm.get("locations") or []) if loc.get("country")})
        status = stm.get("overallStatus", "UNKNOWN")
        trials.append({
            "nct_id": idm.get("nctId", ""),
            "title": idm.get("briefTitle", ""),
            "status": status,
            "status_group": ("Active" if status in ACTIVE_STATUSES else "Completed" if status == "COMPLETED"
                             else "Stopped" if status in STOPPED_STATUSES else "Unknown"),
            "phase": plab, "phase_rank": prank,
            "enrollment": (dsm.get("enrollmentInfo") or {}).get("count"),
            "sponsor": (spm.get("leadSponsor") or {}).get("name", ""),
            "sponsor_class": (spm.get("leadSponsor") or {}).get("class", "UNKNOWN"),
            "start": parse_date((stm.get("startDateStruct") or {}).get("date")),
            "primary_completion": parse_date((stm.get("primaryCompletionDateStruct") or {}).get("date")),
            "primary_completion_raw": (stm.get("primaryCompletionDateStruct") or {}).get("date"),
            "countries": countries,
            "is_drug": bool(drug_ivs),
            "interventions": classified,
            "nondrug_interventions": [iv.get("name") for iv in ivs if iv.get("type") not in DRUG_TYPES],
            "background": background,
            "mechanisms": sorted({c["mechanism"] for c in classified}, key=MECHANISMS.index),
            "modalities": sorted({c["modality"] for c in classified}, key=MODALITIES.index),
            "pairs": sorted({(c["mechanism"], c["modality"]) for c in classified}),
        })
    return trials


def sponsor_group(c):
    return {"INDUSTRY": "Industry", "NIH": "Government", "FED": "Government", "OTHER_GOV": "Government",
            "NETWORK": "Network", "OTHER": "Academic / hospital", "INDIV": "Individual"}.get(c, "Other")


# --------------------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------------------
def analyse(ct, ot, pm, trials):
    today = dt.date.today()
    drug = [t for t in trials if t["is_drug"]]
    active = [t for t in drug if t["status_group"] == "Active"]

    # gap map: mechanism x modality, active drug trials (and all-time totals for context)
    grid_active, grid_all, grid_ids = collections.Counter(), collections.Counter(), collections.defaultdict(list)
    for t in drug:
        for pair in t["pairs"]:
            grid_all[pair] += 1
            if t["status_group"] == "Active":
                grid_active[pair] += 1
                grid_ids[pair].append(t["nct_id"])

    # phase profile among active trials, per mechanism (trial counted once per mechanism)
    phase_cols = ["Ph 1", "Ph 2", "Ph 3", "Ph 4", "N/A"]
    phase_prof = collections.defaultdict(collections.Counter)
    for t in active:
        col = {0: "N/A", 0.5: "Ph 1", 1: "Ph 1", 2: "Ph 2", 3: "Ph 3", 4: "Ph 4"}[t["phase_rank"]]
        for m in t["mechanisms"]:
            phase_prof[m][col] += 1

    # catalyst calendar: active drug trials, Phase >= 2, primary completion between today and end of year+2
    horizon = dt.date(today.year + 2, 12, 31)
    catalysts = sorted([t for t in active if t["phase_rank"] >= 2 and t["primary_completion"]
                        and today <= t["primary_completion"] <= horizon], key=lambda t: t["primary_completion"])
    overdue = [t for t in active if t["phase_rank"] >= 2 and t["primary_completion"] and t["primary_completion"] < today]

    # target evidence table: Open Targets + PubMed + trial counts mapped through MECH_TO_GENES
    gene_to_mech = {g.upper(): m for m, gs in MECH_TO_GENES.items() for g in gs}
    mech_active = collections.Counter(m for t in active for m in t["mechanisms"])
    mech_all = collections.Counter(m for t in drug for m in t["mechanisms"])
    targets = []
    for r in ot.get("rows", []):
        sym = r["target"]["approvedSymbol"]
        dts = {d["id"]: d["score"] for d in (r.get("datatypeScores") or [])}
        mech = gene_to_mech.get(sym.upper())
        pmc = (pm.get("counts") or {}).get(sym) or {}
        targets.append({
            "symbol": sym, "name": r["target"].get("approvedName", ""), "ensembl": r["target"].get("id", ""),
            "ot_score": r.get("score"), "genetic_score": dts.get("genetic_association"),
            "literature_score": dts.get("literature"),
            "mechanism_bucket": mech or "",
            "active_trials": mech_active.get(mech, 0) if mech else 0,
            "all_trials": mech_all.get(mech, 0) if mech else 0,
            "pubmed_total": pmc.get("total"), "pubmed_recent": pmc.get("recent"),
        })

    # sponsors & geography (active drug trials)
    spons = collections.Counter(sponsor_group(t["sponsor_class"]) for t in active)
    countries = collections.Counter(c for t in active for c in t["countries"])
    asia_trials = [t for t in active if any(c in ASIA for c in t["countries"])]
    top_sponsors = collections.Counter(t["sponsor"] for t in active if t["sponsor_class"] == "INDUSTRY")

    tdp = [t for t in drug if "TDP-43 biology" in t["mechanisms"]]

    # unclassified interventions (for the human-in-the-loop overrides file)
    unclassified = collections.Counter()
    for t in drug:
        for c in t["interventions"]:
            if c["mechanism"] == "Unclassified":
                unclassified[c["name"].strip()] += 1

    return {"today": today, "horizon": horizon, "drug": drug, "active": active,
            "grid_active": grid_active, "grid_all": grid_all, "grid_ids": grid_ids,
            "phase_cols": phase_cols, "phase_prof": phase_prof,
            "catalysts": catalysts, "overdue": overdue, "targets": targets,
            "mech_active": mech_active, "mech_all": mech_all,
            "sponsors": spons, "countries": countries, "asia_trials": asia_trials, "top_sponsors": top_sponsors,
            "tdp": tdp, "unclassified": unclassified}


# --------------------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------------------
def validate(ct, trials, a):
    checks = []

    def add(name, status, detail):
        checks.append({"name": name, "status": status, "detail": detail})

    n_fetched, n_total = len(ct["studies"]), ct.get("totalCount")
    add("Record count reconciles with API total", "PASS" if n_total == n_fetched else "FAIL",
        f"fetched {n_fetched} studies; API reported totalCount = {n_total}")

    ids = [t["nct_id"] for t in trials]
    bad = [i for i in ids if not re.fullmatch(r"NCT\d{8}", i)]
    dup = [i for i, c in collections.Counter(ids).items() if c > 1]
    add("Trial IDs well-formed and unique", "PASS" if not bad and not dup else "FAIL",
        f"{len(ids)} interventional trials; {len(bad)} malformed, {len(dup)} duplicated. "
        "Every trial cited on this page links to its own ClinicalTrials.gov record.")

    ivs = [c for t in a["drug"] for c in t["interventions"]]
    n_mech = sum(c["mechanism"] != "Unclassified" for c in ivs)
    n_named = sum(c["mod_rule"] not in (None, "type-fallback") for c in ivs)
    cov = n_mech / len(ivs) if ivs else 0
    add("Mechanism coverage of drug interventions", "PASS" if cov >= 0.6 else "WARN",
        f"{n_mech}/{len(ivs)} ({cov:.0%}) assigned a mechanism; the rest are listed in "
        "data/unclassified_interventions.csv for manual review. "
        f"Modality set by an explicit rule for {n_named}/{len(ivs)}; others fall back to the registry's intervention type.")

    by_id = {t["nct_id"]: t for t in trials}
    for nct, mech, mod, label in KNOWN_TRIALS:
        t = by_id.get(nct)
        if not t:
            add(f"Known trial present: {label}", "WARN", f"{nct} not in this snapshot")
        else:
            ok = (mech, mod) in t["pairs"]
            add(f"Known trial classified correctly: {label}", "PASS" if ok else "FAIL",
                f"{nct}: expected {mech} / {mod}; got {', '.join(f'{m} / {d}' for m, d in t['pairs']) or 'none'}")

    for rx, mech, mod in KNOWN_ASSETS:
        hits = [c for c in ivs if re.search(rx, c["name"], re.I)]
        if not hits:
            add(f"Known asset /{rx}/", "WARN", "no intervention with this name in snapshot")
            continue
        wrong = [c for c in hits if (c["mechanism"], c["modality"]) != (mech, mod)]
        asset_label = re.sub(r"\[[^]]*\]\??", "-", rx.split("|")[0])  # e.g. VTx[- ]?002 -> VTx-002
        add(f"Known asset classified: {asset_label}", "PASS" if not wrong else "FAIL",
            f"{len(hits)} intervention(s); expected {mech} / {mod}; {len(wrong)} mismatched")

    no_date = sum(1 for t in a["active"] if not t["primary_completion"])
    add("Active drug trials with a parseable primary-completion date",
        "PASS" if no_date == 0 else "WARN", f"{len(a['active']) - no_date}/{len(a['active'])}")

    ev_path = os.path.join(DATA_DIR, "validation_results.json")
    manual = None
    if os.path.exists(ev_path):
        with open(ev_path) as f:
            manual = json.load(f)
        add("Manual label agreement (mechanism)", "PASS" if manual["mechanism_acc"] >= 0.85 else "WARN",
            f"{manual['mechanism_acc']:.0%} on {manual['n']} hand-labelled trials")
        add("Manual label agreement (modality)", "PASS" if manual["modality_acc"] >= 0.85 else "WARN",
            f"{manual['modality_acc']:.0%} on {manual['n']} hand-labelled trials")
    else:
        add("Manual label agreement", "TODO",
            f"label {LABEL_SAMPLE_SIZE} random trials in data/label_sample.csv, then run --evaluate")
    return checks, manual


def write_label_sample(a):
    path = os.path.join(DATA_DIR, "label_sample.csv")
    if os.path.exists(path):
        return  # never overwrite manual work
    rng = random.Random(RANDOM_SEED)
    pool = sorted(a["drug"], key=lambda t: t["nct_id"])
    sample = rng.sample(pool, min(LABEL_SAMPLE_SIZE, len(pool)))
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["nct_id", "title", "interventions", "auto_mechanism", "auto_modality",
                    "manual_mechanism", "manual_modality", "notes"])
        for t in sample:
            w.writerow([t["nct_id"], t["title"], "; ".join(c["name"] for c in t["interventions"]),
                        "; ".join(t["mechanisms"]), "; ".join(t["modalities"]), "", "", ""])


def evaluate_labels():
    path = os.path.join(DATA_DIR, "label_sample.csv")
    if not os.path.exists(path):
        sys.exit("No data/label_sample.csv yet - run the pipeline once first.")
    n = mech_ok = mod_ok = 0
    errors = []
    norm = lambda s: {x.strip().lower() for x in (s or "").split(";") if x.strip()}
    # Score against the CURRENT classifier output (data/trials.csv), so rule fixes and overrides
    # are reflected; the auto_* columns in label_sample.csv are only a snapshot from the first run.
    current = {}
    trials_path = os.path.join(DATA_DIR, "trials.csv")
    if os.path.exists(trials_path):
        with open(trials_path, newline="") as f:
            for t in csv.DictReader(f):
                current[t["nct_id"]] = t
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if not (r["manual_mechanism"].strip() or r["manual_modality"].strip()):
                continue
            cur = current.get(r["nct_id"])
            if cur:
                r["auto_mechanism"], r["auto_modality"] = cur["mechanisms"], cur["modalities"]
            n += 1
            m_ok = norm(r["manual_mechanism"]) == norm(r["auto_mechanism"])
            d_ok = norm(r["manual_modality"]) == norm(r["auto_modality"])
            mech_ok += m_ok
            mod_ok += d_ok
            if not (m_ok and d_ok):
                errors.append({"nct_id": r["nct_id"], "auto": [r["auto_mechanism"], r["auto_modality"]],
                               "manual": [r["manual_mechanism"], r["manual_modality"]]})
    if n == 0:
        sys.exit("No rows labelled yet: fill manual_mechanism / manual_modality in data/label_sample.csv")
    res = {"n": n, "mechanism_acc": mech_ok / n, "modality_acc": mod_ok / n, "errors": errors,
           "evaluated_at": now_iso()}
    with open(os.path.join(DATA_DIR, "validation_results.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(f"Labelled rows: {n}\n  mechanism agreement: {res['mechanism_acc']:.0%}\n"
          f"  modality agreement:  {res['modality_acc']:.0%}\n  disagreements: {len(errors)}")
    for e in errors[:15]:
        print(f"   {e['nct_id']}: auto={e['auto']}  manual={e['manual']}")
    print("Re-run the pipeline to show these results on the dashboard.")


# --------------------------------------------------------------------------------------
# Observations (deterministic, data-derived; no LLM)
# --------------------------------------------------------------------------------------
def observations(a, pm):
    obs = []
    ma = {m: c for m, c in a["mech_active"].items() if m != "Unclassified"}
    if ma:
        top = sorted(ma.items(), key=lambda x: -x[1])[:3]
        obs.append("Most active mechanisms (active drug trials): "
                   + ", ".join(f"<b>{esc(m)}</b> ({c})" for m, c in top) + ".")
    gen = ["SOD1", "C9orf72", "FUS", "ATXN2", "TDP-43 biology"]
    gen_active = sum(1 for t in a["active"] if set(t["mechanisms"]) & set(gen))
    if a["active"]:
        share = gen_active / len(a["active"])
        tail = ("most activity still targets broad downstream pathways." if share < 0.5
                else "the field has shifted towards genetically defined targets.")
        obs.append(f"Genetically defined or target-specific programmes (SOD1, C9orf72, FUS, ATXN2, TDP-43 biology) "
                   f"account for <b>{gen_active}</b> of <b>{len(a['active'])}</b> active drug trials "
                   f"({share:.0%}); {tail}")
    nm = collections.Counter(m for t in a["active"] for m in t["modalities"])
    adv = sum(nm[m] for m in ["ASO", "siRNA / RNAi", "AAV gene therapy", "Vectorized antibody", "Antibody", "Cell therapy"])
    if a["active"]:
        obs.append(f"Advanced modalities (oligonucleotides, gene therapy, antibodies, cells) appear in <b>{adv}</b> "
                   f"active trials; small molecules in <b>{nm['Small molecule']}</b>.")
    tdp_active = [t for t in a["tdp"] if t["status_group"] == "Active"]
    obs.append("TDP-43 pathology is present in the large majority of ALS cases; "
               f"<b>{len(a['tdp'])}</b> drug trials all-time target this biology, <b>{len(tdp_active)}</b> of them active "
               "(see deep-dive).")
    gaps = [t for t in a["targets"] if (t["genetic_score"] or 0) >= 0.5 and t["active_trials"] == 0]
    if gaps:
        obs.append("Strong genetic evidence but no active trial mapped in this snapshot: "
                   + ", ".join(f"<b>{esc(t['symbol'])}</b>" for t in gaps[:8])
                   + " - candidate white space (or programmes still preclinical).")
    if a["catalysts"]:
        by_year = collections.Counter(t["primary_completion"].year for t in a["catalysts"])
        obs.append("Phase 2+ primary-completion dates in the window: "
                   + ", ".join(f"{y}: <b>{c}</b>" for y, c in sorted(by_year.items()))
                   + ". These are potential data readouts to track.")
    if a["active"]:
        obs.append(f"<b>{len(a['asia_trials'])}</b> active drug trials ({len(a['asia_trials']) / len(a['active']):.0%}) "
                   "have at least one site in Asia.")
    return obs


# --------------------------------------------------------------------------------------
# HTML rendering
# --------------------------------------------------------------------------------------
def esc(x):
    return html.escape("" if x is None else str(x))


def nct_link(i):
    return f'<a href="https://clinicaltrials.gov/study/{esc(i)}" target="_blank" rel="noopener">{esc(i)}</a>'


SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
       "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]


def seq_color(v, vmax):
    if not v:
        return None
    i = min(len(SEQ) - 1, int(round((v / vmax) ** 0.6 * (len(SEQ) - 1)))) if vmax else 0
    return SEQ[i]


def render_heatmap(a):
    mechs = [m for m in MECHANISMS if any(a["grid_all"][(m, d)] for d in MODALITIES)]
    mods = [d for d in MODALITIES if any(a["grid_all"][(m, d)] for m in MECHANISMS)]
    vmax = max([a["grid_active"][(m, d)] for m in mechs for d in mods] + [1])
    rows = []
    for m in mechs:
        cells = []
        for d in mods:
            v, tot = a["grid_active"][(m, d)], a["grid_all"][(m, d)]
            col = seq_color(v, vmax)
            ink = "#fff" if col and SEQ.index(col) >= 6 else "var(--text-primary)"
            tip = f"{m} × {d}: {v} active / {tot} all-time"
            ids = a["grid_ids"][(m, d)][:12]
            if ids:
                tip += " | " + ", ".join(ids) + (" …" if len(a["grid_ids"][(m, d)]) > 12 else "")
            style = f"background:{col};color:{ink}" if col else ""
            label = str(v) if v else ("·" if tot else "")
            sub = f'<span class="sub">{tot}</span>' if tot else ""
            cells.append(f'<td class="hm" style="{style}" title="{esc(tip)}">{label}{sub}</td>')
        tot_row = a["mech_active"].get(m, 0)
        rows.append(f'<tr><th class="rowh">{esc(m)}</th>{"".join(cells)}<td class="tot">{tot_row}</td></tr>')
    head = "".join(f"<th class='colh'>{esc(d)}</th>" for d in mods)
    return (f'<table class="heat"><thead><tr><th></th>{head}<th class="colh">Active trials</th></tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table>'
            '<p class="note">Big number = active trials; small number = all-time trials. "·" = past activity only. '
            'Hover a cell for trial IDs. A trial testing several interventions counts in each cell it touches.</p>')


def render_phase(a):
    mechs = [m for m in MECHANISMS if a["phase_prof"].get(m)]
    head = "".join(f"<th>{c}</th>" for c in a["phase_cols"])
    rows = []
    for m in mechs:
        pr = a["phase_prof"][m]
        cells = "".join(f"<td class='num'>{pr.get(c, 0) or ''}</td>" for c in a["phase_cols"])
        rows.append(f"<tr><th class='rowh'>{esc(m)}</th>{cells}</tr>")
    return f'<table class="data compact"><thead><tr><th>Mechanism</th>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table>'


def iv_summary(t, n=3):
    names = [c["name"] for c in t["interventions"]]
    s = "; ".join(names[:n]) + (" …" if len(names) > n else "")
    return s


def render_trial_table(trials, date_col=True, empty="No trials in this snapshot."):
    if not trials:
        return f"<p class='note'>{empty}</p>"
    rows = []
    for t in trials:
        rows.append(
            "<tr>"
            + (f"<td class='nowrap'>{esc(t['primary_completion_raw'] or '')}</td>" if date_col else "")
            + f"<td class='nowrap'>{nct_link(t['nct_id'])}</td>"
            f"<td>{esc(iv_summary(t))}</td>"
            f"<td>{esc(', '.join(t['mechanisms']))}</td>"
            f"<td>{esc(', '.join(t['modalities']))}</td>"
            f"<td class='nowrap'>{esc(t['phase'])}</td>"
            f"<td>{esc(t['sponsor'])}</td>"
            f"<td class='num'>{esc(t['enrollment'] or '')}</td>"
            f"<td class='nowrap'>{esc(t['status'].replace('_', ' ').title())}</td>"
            "</tr>")
    head = ("<th>Primary completion</th>" if date_col else "") + (
        "<th>Trial</th><th>Intervention(s)</th><th>Mechanism</th><th>Modality</th><th>Phase</th>"
        "<th>Lead sponsor</th><th>n</th><th>Status</th>")
    return f'<div class="scroll"><table class="data"><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>'


def bar(v, vmax, w=120):
    px = 0 if not vmax or v is None else max(2, round(v / vmax * w))
    return f'<span class="bar" style="width:{px}px"></span>'


def render_targets(a):
    if not a["targets"]:
        return "<p class='note'>Open Targets data unavailable.</p>"
    tmax = max((t["ot_score"] or 0) for t in a["targets"]) or 1
    rows = []
    for t in a["targets"]:
        gs = t["genetic_score"]
        flag = ""
        if (gs or 0) >= 0.5 and t["active_trials"] == 0:
            flag = "<span class='tag'>open</span>"
        rows.append(
            f"<tr><td><b>{esc(t['symbol'])}</b><div class='sub2'>{esc(t['name'][:48])}</div></td>"
            f"<td class='nowrap'>{bar(t['ot_score'], tmax)} {(t['ot_score'] or 0):.2f}</td>"
            f"<td class='num'>{'' if gs is None else f'{gs:.2f}'}</td>"
            f"<td class='num'>{'' if t['pubmed_total'] is None else t['pubmed_total']}</td>"
            f"<td class='num'>{'' if t['pubmed_recent'] is None else t['pubmed_recent']}</td>"
            f"<td>{esc(t['mechanism_bucket'])}</td>"
            f"<td class='num'>{t['active_trials'] if t['mechanism_bucket'] else '–'}</td>"
            f"<td>{flag}</td></tr>")
    return ('<div class="scroll"><table class="data"><thead><tr><th>Target</th><th>Open Targets overall</th>'
            '<th>Genetic evidence</th><th>PubMed (all)</th><th>PubMed (5 y)</th><th>Mapped bucket</th>'
            '<th>Active trials</th><th></th></tr></thead><tbody>' + "".join(rows) + "</tbody></table></div>"
            '<p class="note">Scores 0–1 from Open Targets. PubMed counts = title/abstract co-mentions with ALS. '
            '"open" = genetic evidence ≥ 0.5 and no active trial mapped to that gene. That means "check this", '
            'not "it\'s a gap": the trial-to-gene mapping covers only the buckets listed.</p>')


def render_bars(counter, n=10):
    items = counter.most_common(n)
    if not items:
        return "<p class='note'>No data.</p>"
    vmax = items[0][1]
    return "<table class='bars'>" + "".join(
        f"<tr><th>{esc(k)}</th><td>{bar(v, vmax, 160)}</td><td class='num'>{v}</td></tr>" for k, v in items) + "</table>"


def render_checks(checks):
    icon = {"PASS": "✔", "FAIL": "✖", "WARN": "!", "TODO": "…"}
    return "<table class='data checks'><thead><tr><th></th><th>Check</th><th>Result</th></tr></thead><tbody>" + "".join(
        f"<tr><td><span class='st st-{c['status'].lower()}'>{icon[c['status']]} {c['status']}</span></td>"
        f"<td>{esc(c['name'])}</td><td>{esc(c['detail'])}</td></tr>" for c in checks) + "</tbody></table>"


CSS = """
.viz-root{color-scheme:light;--surface-0:#f6f5f2;--surface-1:#fcfcfb;--text-primary:#0b0b0b;--text-secondary:#52514e;
--muted:#7a7974;--rule:#e4e2dc;--accent:#2a78d6;--good:#0ca30c;--warn:#b07800;--crit:#d03b3b}
@media (prefers-color-scheme:dark){.viz-root{color-scheme:dark;--surface-0:#121211;--surface-1:#1a1a19;
--text-primary:#f0efec;--text-secondary:#c3c2b7;--muted:#9a998f;--rule:#33332f;--accent:#3987e5}}
*{box-sizing:border-box}body{margin:0}
.viz-root{background:var(--surface-0);color:var(--text-primary);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",
Inter,Helvetica,Arial,sans-serif;padding:32px 20px 64px}
.wrap{max-width:1180px;margin:0 auto}
h1{font-size:26px;margin:0 0 4px;letter-spacing:-.01em}h2{font-size:18px;margin:0 0 4px}
.lede{color:var(--text-secondary);margin:0 0 20px;max-width:820px}
.meta{color:var(--muted);font-size:12px;margin-bottom:20px}
.card{background:var(--surface-1);border:1px solid var(--rule);border-radius:10px;padding:20px 22px;margin:0 0 18px}
.card .why{color:var(--text-secondary);margin:0 0 14px;font-size:13px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;margin-bottom:18px}
.kpi{background:var(--surface-1);border:1px solid var(--rule);border-radius:10px;padding:14px 16px}
.kpi .v{font-size:28px;font-weight:650;line-height:1.1}.kpi .l{color:var(--text-secondary);font-size:12px}
ul.obs{margin:0;padding-left:18px}ul.obs li{margin:6px 0}
table{border-collapse:collapse;width:100%}
table.data th,table.data td{border-bottom:1px solid var(--rule);padding:6px 8px;text-align:left;vertical-align:top;font-size:13px}
table.data thead th{color:var(--text-secondary);font-weight:600;font-size:12px;position:sticky;top:0;background:var(--surface-1)}
table.compact{width:auto}table.compact td,table.compact th{padding:4px 12px}
.num{text-align:right;font-variant-numeric:tabular-nums}.nowrap{white-space:nowrap}
.scroll{max-height:520px;overflow:auto;border:1px solid var(--rule);border-radius:8px}
table.heat{width:auto;border-collapse:separate;border-spacing:2px}
table.heat td.hm{width:86px;height:40px;text-align:center;font-weight:650;font-size:15px;border-radius:4px;
background:transparent;outline:1px solid var(--rule);outline-offset:-1px;cursor:default;font-variant-numeric:tabular-nums}
table.heat td.hm .sub{display:block;font-size:10px;font-weight:400;opacity:.75}
table.heat th.colh{font-size:11px;color:var(--text-secondary);font-weight:600;vertical-align:bottom;padding:0 2px 6px;max-width:86px}
table.heat th.rowh{text-align:right;font-size:12px;padding-right:10px;font-weight:500;white-space:nowrap}
table.heat td.tot{padding-left:10px;font-weight:650;color:var(--text-secondary);font-variant-numeric:tabular-nums}
.rowh{font-weight:500}
.note{color:var(--muted);font-size:12px;margin:10px 0 0}
.sub2{color:var(--muted);font-size:11px}
.bar{display:inline-block;height:10px;background:var(--accent);border-radius:0 4px 4px 0;vertical-align:middle}
table.bars{width:auto}table.bars th{text-align:right;font-weight:500;padding:3px 10px 3px 0;font-size:13px}
table.bars td{padding:3px 6px}
.two{display:grid;grid-template-columns:1fr 1fr;gap:18px}@media(max-width:820px){.two{grid-template-columns:1fr}}
.tag{font-size:11px;border:1px solid var(--accent);color:var(--accent);border-radius:10px;padding:0 7px}
.st{font-size:11px;font-weight:650;white-space:nowrap}.st-pass{color:var(--good)}.st-fail{color:var(--crit)}
.st-warn{color:var(--warn)}.st-todo{color:var(--muted)}
.banner{background:#fff3cd;color:#5c4400;border:1px solid #e9c46a;border-radius:8px;padding:10px 14px;margin-bottom:16px}
a{color:var(--accent)}
footer{color:var(--muted);font-size:12px;margin-top:28px}
"""


def render_html(ct, ot, pm, trials, a, checks, obs, fixture):
    today = a["today"]
    n_inter = len(trials)
    kpis = [
        (n_inter, "interventional ALS trials (all-time)"),
        (len(a["drug"]), "testing a drug, biologic or gene therapy"),
        (len(a["active"]), "of those, currently active"),
        (len(a["catalysts"]), f"Phase 2+ primary completions to {a['horizon'].year}"),
    ]
    kpi_html = "".join(f"<div class='kpi'><div class='v'>{v}</div><div class='l'>{esc(l)}</div></div>" for v, l in kpis)
    banner = ("<div class='banner'><b>TEST FIXTURE DATA.</b> This dashboard was built from a synthetic sample "
              "to test the pipeline. Run the script with internet access to build it from live data.</div>") if fixture else ""
    tdp_text = ("TDP-43 mislocalisation and aggregation appear in the large majority of ALS cases, both sporadic "
                "and familial. That makes it the broadest disease-biology target in ALS, but it is hard to drug. "
                "Programmes here include direct approaches (antibodies, including AAV-vectorised antibodies) and "
                "approaches that restore TDP-43 function downstream (e.g. STMN2). "
                "This view lists every trial the rules place in this bucket.")
    sections = f"""
<div class="card"><h2>Key observations</h2>
<p class="why">Generated by fixed rules from the data below (no language model). Check each one against the linked trials before using it.</p>
<ul class="obs">{''.join(f'<li>{o}</li>' for o in obs)}</ul></div>

<div class="card"><h2>1 · Gap map: mechanism × modality</h2>
<p class="why">Where is the field crowded and where is it open? Darker cells = more active drug trials.</p>
{render_heatmap(a)}</div>

<div class="two">
<div class="card"><h2>2 · Phase profile of active trials</h2>
<p class="why">How mature is each mechanism? Highest phase per active trial.</p>{render_phase(a)}</div>
<div class="card"><h2>3 · Who is running trials, and where</h2>
<p class="why">Sponsor type, top industry sponsors and top countries for active drug trials.
{len(a['asia_trials'])} of {len(a['active'])} have an Asian site.</p>
{render_bars(a['sponsors'], 6)}<div style="height:12px"></div>{render_bars(a['top_sponsors'], 8)}
<div style="height:12px"></div>{render_bars(a['countries'], 10)}</div>
</div>

<div class="card"><h2>4 · Upcoming results: Phase 2+ primary completions, {today:%b %Y} – {a['horizon']:%b %Y}</h2>
<p class="why">These dates are when data could read out. Registry dates often slip, so treat them as estimates.
A further {len(a['overdue'])} active Phase 2+ trials are past their registered primary-completion date and show no update.</p>
{render_trial_table(a['catalysts'])}</div>

<div class="card"><h2>5 · Target evidence vs. clinical activity</h2>
<p class="why">Is clinical activity going where the genetic evidence points? Top {len(a['targets'])} ALS-associated targets
from Open Targets ({esc(ot.get('disease', {}).get('id', ''))}), with literature volume and mapped trial activity.</p>
{render_targets(a)}</div>

<div class="card"><h2>6 · Deep-dive: TDP-43 biology</h2>
<p class="why">{tdp_text}</p>
{render_trial_table(sorted(a['tdp'], key=lambda t: (t['status_group'] != 'Active', -(t['phase_rank'] or 0))), date_col=True)}</div>

<div class="card"><h2>7 · Method & validation</h2>
<p class="why">Each drug intervention is classified by an ordered list of rules: {len(NAMED_RULES)} rules for named assets,
then {len(KEYWORD_RULES)} keyword rules, then a fallback to the registry's intervention type. The rule that fired is
saved in data/trials.csv, so every label can be traced. Placebo/sham arms and non-drug interventions are excluded from the landscape.</p>
{render_checks(checks)}
<p class="note">Limitations: ClinicalTrials.gov under-represents trials registered only in other registries (e.g. jRCT in Japan, ChiCTR in China, EU CTIS).
Preclinical programmes are not visible. Rule-based labels are approximate until the manual-validation step is completed.
Gene mapping to Open Targets covers only SOD1, C9orf72, FUS, ATXN2 and TARDBP/STMN2. This is not investment advice.</p></div>
"""
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>ALS Landscape Navigator</title><style>{CSS}</style></head>
<body><div class="viz-root"><div class="wrap">
{banner}
<h1>ALS Landscape Navigator</h1>
<p class="lede">The ALS therapeutic landscape as seen in public data: which mechanisms and modalities are crowded or open,
which results are due soon, and whether clinical activity follows the genetic evidence.</p>
<div class="meta">Data snapshot: ClinicalTrials.gov {esc(ct.get('fetched_at', '')[:10])} · Open Targets {esc(ot.get('fetched_at', '')[:10])}
· PubMed {esc(pm.get('fetched_at', '')[:10])} · built {today.isoformat()} · query: "{esc(CONDITION)}"</div>
<div class="kpis">{kpi_html}</div>
{sections}
<footer>Sources: ClinicalTrials.gov API v2 · Open Targets Platform GraphQL API · NCBI PubMed E-utilities.
Built with a standard-library Python pipeline (als_navigator.py); raw responses are cached for reproducibility.</footer>
</div></div></body></html>"""


# --------------------------------------------------------------------------------------
# CSV outputs
# --------------------------------------------------------------------------------------
def write_csvs(trials, a):
    with open(os.path.join(DATA_DIR, "trials.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["nct_id", "title", "status", "status_group", "phase", "enrollment", "lead_sponsor",
                    "sponsor_class", "start", "primary_completion", "countries", "is_drug", "interventions",
                    "mechanisms", "modalities", "rules_fired", "nondrug_interventions"])
        for t in trials:
            w.writerow([t["nct_id"], t["title"], t["status"], t["status_group"], t["phase"], t["enrollment"],
                        t["sponsor"], t["sponsor_class"], t["start"] or "", t["primary_completion_raw"] or "",
                        "; ".join(t["countries"]), t["is_drug"],
                        "; ".join(c["name"] for c in t["interventions"]),
                        "; ".join(t["mechanisms"]), "; ".join(t["modalities"]),
                        "; ".join(f"{c['mech_rule']}/{c['mod_rule']}" for c in t["interventions"]),
                        "; ".join(n or "" for n in t["nondrug_interventions"])])
    with open(os.path.join(DATA_DIR, "targets.csv"), "w", newline="") as f:
        cols = ["symbol", "name", "ensembl", "ot_score", "genetic_score", "literature_score",
                "mechanism_bucket", "active_trials", "all_trials", "pubmed_total", "pubmed_recent"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for t in a["targets"]:
            w.writerow({k: t.get(k) for k in cols})
    with open(os.path.join(DATA_DIR, "unclassified_interventions.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["intervention_name", "n_trials", "mechanism", "modality"])
        for name, n in a["unclassified"].most_common():
            w.writerow([name, n, "", ""])


# --------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true", help="ignore cache and re-download everything")
    ap.add_argument("--evaluate", action="store_true", help="score manual labels in data/label_sample.csv")
    args = ap.parse_args()
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    if args.evaluate:
        evaluate_labels()
        return

    try:
        ct = load_or_fetch("ctgov_studies.json", fetch_ctgov, args.refresh)
    except Exception as e:
        sys.exit(f"\nCould not download from ClinicalTrials.gov:\n  {e}\n"
                 "Check your internet connection (a VPN or corporate proxy can block it) and try again.")
    try:
        ot = load_or_fetch("opentargets.json", fetch_opentargets, args.refresh)
    except Exception as e:
        print(f"[Open Targets] skipped: {e}")
        ot = {"fetched_at": "", "rows": []}
    symbols = [r["target"]["approvedSymbol"] for r in ot.get("rows", [])][:N_PUBMED_TARGETS]
    for must in ["SOD1", "C9orf72", "FUS", "TARDBP", "ATXN2"]:
        if must not in symbols:
            symbols.append(must)
    try:
        pm = load_or_fetch("pubmed_counts.json", lambda: fetch_pubmed(symbols), args.refresh)
    except Exception as e:
        print(f"[PubMed] skipped: {e}")
        pm = {"fetched_at": "", "counts": {}}

    fixture = bool(ct.get("_fixture"))
    trials = normalise_trials(ct, load_overrides())
    a = analyse(ct, ot, pm, trials)
    checks, _ = validate(ct, trials, a)
    obs = observations(a, pm)
    write_csvs(trials, a)
    write_label_sample(a)

    out = os.path.join(OUT_DIR, "als_landscape.html")
    with open(out, "w") as f:
        f.write(render_html(ct, ot, pm, trials, a, checks, obs, fixture))

    print("\nDone.")
    print(f"  interventional trials: {len(trials)} | drug trials: {len(a['drug'])} | active: {len(a['active'])}")
    for c in checks:
        print(f"  [{c['status']}] {c['name']}: {c['detail']}")
    print(f"\n  Dashboard: {out}")
    print(f"  Tables:    {DATA_DIR}/trials.csv, targets.csv, unclassified_interventions.csv, label_sample.csv")


if __name__ == "__main__":
    main()
