#!/usr/bin/env python3
"""
ALS Landscape Navigator - exploratory scientific-intelligence prototype (v2, repaired)
====================================================================================

Combines three public sources into one self-contained HTML report plus CSV tables:

  * ClinicalTrials.gov API v2  -> ALS interventional studies (phase, status, arms, interventions, dates)
  * Open Targets Platform      -> ALS-associated targets with association scores (prioritisation evidence)
  * PubMed E-utilities         -> title/abstract co-mention counts per target (literature volume only)

Every derived label keeps its provenance: the rule or curated source behind it, the evidence tier
(curated asset / keyword / registry-type fallback / unresolved), and the intervention's role in the
trial (experimental, background, comparator, control-arm, placebo, non-therapeutic, unresolved).

This is an exploratory prototype. Classification accuracy has NOT been measured against a
human-reviewed reference set; see data/validation/ and README.md.

Standard library only.

Usage
-----
  python3 als_navigator.py                          # build from the cached snapshot in data/raw/
  python3 als_navigator.py --analysis-date 2026-11-18
  python3 als_navigator.py --refresh                # re-download (changes the inputs; new snapshot)
  python3 als_navigator.py --evaluate               # score HUMAN-approved labels in the review pack
  python3 -m unittest discover -s tests -v          # regression checks
"""

import argparse
import calendar
import collections
import csv
import datetime as dt
import hashlib
import html
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

try:  # macOS python.org builds often lack root certificates; certifi fixes that
    import ssl
    import certifi
    _SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except Exception:
    _SSL_CTX = None

# --------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------
CLASSIFICATION_VERSION = "2.0.0"
CONDITION = "amyotrophic lateral sclerosis"
CTGOV_URL = "https://clinicaltrials.gov/api/v2/studies"
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
USER_AGENT = "ALS-Landscape-Navigator/2.0 (research prototype; python urllib)"

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
OUT_DIR = os.path.join(HERE, "output")
CURATION_PATH = os.path.join(DATA_DIR, "curation", "asset_curation.csv")
VALID_DIR = os.path.join(DATA_DIR, "validation")

N_OT_TARGETS = 30          # top Open Targets associations pulled (selection limit)
N_PUBMED_TARGETS = 20      # of those, how many get PubMed counts (+ REQUIRED_PUBMED_GENES)
REQUIRED_PUBMED_GENES = ["SOD1", "C9orf72", "FUS", "TARDBP", "ATXN2"]
MILESTONE_HORIZON_YEARS = 2  # window = analysis date .. 31 Dec of (analysis year + 2)

QUERY_DESCRIPTION = {
    "ClinicalTrials.gov": f'API v2 /studies, query.cond="{CONDITION}", pageSize=1000, all pages; '
                          "analysis keeps studyType=INTERVENTIONAL only",
    "Open Targets": f'search "{CONDITION}" -> disease ID; associatedTargets top {N_OT_TARGETS} by overall score',
    "PubMed": f'esearch counts, "{CONDITION}"[tiab] AND <gene>[tiab] (aliases for TARDBP, C9orf72); '
              f"top {N_PUBMED_TARGETS} Open Targets genes plus {', '.join(REQUIRED_PUBMED_GENES)}; "
              "recent window = last 5 calendar years by publication date",
}

ACTIVE_STATUSES = {"RECRUITING", "NOT_YET_RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION"}
STOPPED_STATUSES = {"TERMINATED", "WITHDRAWN", "SUSPENDED"}
DRUG_TYPES = {"DRUG", "BIOLOGICAL", "GENETIC", "COMBINATION_PRODUCT", "DIETARY_SUPPLEMENT"}

# Names that identify standard-of-care background agents (anchored at the start of the name).
SOC_RE = re.compile(r"^\s*(riluzole|rilutek|edaravone|radicava|standard[- ]of[- ]care|standard care|"
                    r"best supportive care|usual care|background therapy)\b", re.I)
PLACEBO_RE = re.compile(r"\b(placebos?|sham|vehicle|matching control|saline control|dummy)\b", re.I)
# Interventions that are registered with a drug-like type but are sampling, imaging or enabling agents.
NONTHERAPEUTIC_RE = re.compile(
    r"blood sample|biological sample|sample collection|genetic (analysis|testing)|\bbiopsy\b|survey|"
    r"\bSPECT\b|\bPET\b|radiotracer|\btracer\b|\[\d+[A-Z][a-z]?\]|\b(89Zr|18F|11C|123I|99mTc)\b|DATSCAN|"
    r"contrast agent|microspheres?\b", re.I)
_PLACEBO_FILLER = set("""tablet tablets capsule capsules injection injections matching match matched mg ml daily oral
orally solution dose doses dosing comparator comparators group arm for of the to or and with subcutaneous intravenous iv
infusion infusions control saline normal sugar pill pills placebo placebos sham vehicle dummy identical look alike
lookalike once twice weekly per day film suspension powder drink administration administered bid tid qd sc it
intrathecal procedure""".split())


def is_placebo_name(name):
    """True if the name describes only a placebo/sham/vehicle. 'Drug X or placebo' (blinded naming) is NOT placebo."""
    if not PLACEBO_RE.search(name or ""):
        return False
    if re.match(r"^\W*(matching\s+|identical\s+)?(placebos?|sham|vehicle|dummy)\b", name, re.I):
        return True  # e.g. "Placebo to match riluzole", "Sham procedure"
    tokens = re.findall(r"[A-Za-z][A-Za-z0-9\-]*", name)
    substantive = [t for t in tokens if t.lower() not in _PLACEBO_FILLER and len(t) >= 3]
    return not substantive


SUPPORTIVE_RE = re.compile(r"prophyla|premedicat|pre-medicat|immunosuppress|immune suppression regimen|"
                           r"anti-?emetic|rescue medication", re.I)
COMBO_SPLIT_RE = re.compile(r"\s*(?:\+|\bplus\b|\bcombined with\b|\bin combination with\b)\s*", re.I)

ROLE_EXPERIMENTAL = "experimental"
ROLE_EXPERIMENTAL_INFERRED = "experimental (inferred)"
ROLE_UNRESOLVED = "unresolved role"
ROLE_BACKGROUND = "background"
ROLE_COMPARATOR = "comparator"
ROLE_CONTROL = "control-arm treatment"
ROLE_PLACEBO = "placebo/sham"
ROLE_NONTHERAPEUTIC = "non-therapeutic"
TEST_ROLES = {ROLE_EXPERIMENTAL, ROLE_EXPERIMENTAL_INFERRED, ROLE_UNRESOLVED}

MODALITIES = ["ASO", "siRNA / RNAi", "AAV gene therapy", "Gene therapy (non-AAV/unspecified)",
              "Vectorized antibody", "Antibody", "Cell therapy", "Other biologic", "Inorganic nanocrystal",
              "Small molecule", "Supplement", "Unresolved"]
MECHANISMS = ["SOD1", "C9orf72", "FUS", "ATXN2", "TDP-43 biology", "Other gene-directed",
              "Neuroinflammation & immune", "Oxidative stress", "Excitotoxicity",
              "Mitochondria & cell stress", "Integrated stress response", "Proteostasis & autophagy",
              "Neurotrophic & cell-based", "Muscle & NMJ", "Other neuroprotective",
              "Symptomatic", "Unclassified"]
GENE_DIRECTED_BUCKETS = ["SOD1", "C9orf72", "FUS", "ATXN2", "TDP-43 biology", "Other gene-directed"]
# Open Targets gene -> mechanism bucket used ONLY for the labelled bucket-level proxy column.
BUCKET_FOR_GENE = {"SOD1": "SOD1", "C9ORF72": "C9orf72", "FUS": "FUS", "ATXN2": "ATXN2", "TARDBP": "TDP-43 biology"}

PUBMED_ALIASES = {"TARDBP": '(TARDBP[tiab] OR "TDP-43"[tiab] OR TDP43[tiab])',
                  "C9orf72": '(C9orf72[tiab] OR "C9ORF72"[tiab])'}
ASIA = {"Japan", "China", "Korea, Republic of", "South Korea", "Taiwan", "Singapore", "Hong Kong",
        "India", "Thailand", "Malaysia", "Vietnam", "Viet Nam", "Philippines", "Indonesia"}

# Keyword rules: (rule_id, regex, mechanism or None, modality or None).
# Applied to the intervention NAME + otherNames. Modality keywords may also be applied to the
# intervention's own description (weaker tier). Mechanism keywords found only in the description are
# stored as review CANDIDATES and never assigned, because descriptions often mention other arms,
# background drugs or disease biology (e.g. "SOD1 mouse model").
KEYWORD_RULES = [
    ("kw-ASO", r"antisense|oligonucleotide|\b\w+nersen\b", None, "ASO"),
    ("kw-siRNA", r"\bsiRNA\b|\bRNAi\b|\b\w+siran\b", None, "siRNA / RNAi"),
    ("kw-AAV", r"\br?AAV\d*\b|adeno[- ]associated", None, "AAV gene therapy"),
    ("kw-gene-therapy", r"gene therapy|gene transfer|plasmid", None, "Gene therapy (non-AAV/unspecified)"),
    ("kw-cell", r"stem cell|mesenchymal|\bMSCs?\b|cell therapy|autologous .*cells|bone marrow|"
                r"\bT[- ]cells?\b|cells? (injection|infusion|transplant)", None, "Cell therapy"),
    ("kw-mab", r"monoclonal|antibod|\b\w+mab\b|immunoglobulin|\bIVIG\b", None, "Antibody"),
    ("kw-SOD1", r"\bSOD[- ]?1\b|superoxide dismutase", "SOD1", None),
    ("kw-C9", r"C9orf72|\bC9\b", "C9orf72", None),
    ("kw-FUS", r"\bFUS\b", "FUS", None),
    ("kw-ATXN2", r"ATXN2|ataxin", "ATXN2", None),
    ("kw-TDP43", r"TDP[- ]?43|TARDBP|STMN2|stathmin", "TDP-43 biology", None),
    ("kw-inflam", r"inflamm|microglia|complement|immun|cytokine|\bTNF\b", "Neuroinflammation & immune", None),
    ("kw-oxid", r"oxidat|antioxid|free radical|glutathione|acetylcysteine", "Oxidative stress", None),
    ("kw-glut", r"glutamat|excitotox|\bAMPA\b|\bNMDA\b|perampanel|fycompa|memantine|talampanel",
     "Excitotoxicity", None),
    ("kw-mito", r"mitochond|bioenerget|coenzyme ?Q|CoQ10|creatine|ketogenic", "Mitochondria & cell stress", None),
    ("kw-ISR", r"integrated stress|\bISR\b|eIF2B", "Integrated stress response", None),
    ("kw-proteo", r"autophag|proteasom|heat shock|\bHSP\d*\b|chaperone", "Proteostasis & autophagy", None),
    ("kw-neurotroph", r"neurotroph|GDNF|BDNF|VEGF|IGF[- ]?1|growth factor", "Neurotrophic & cell-based", None),
    ("kw-symptom", r"spasticity|\bcramps?\b|sialorrh|drool|\bpain\b|pseudobulbar|dyspn", "Symptomatic", None),
]
_KW = [(rid, re.compile(rx, re.I), mech, mod) for rid, rx, mech, mod in KEYWORD_RULES]
TYPE_FALLBACK_MODALITY = {"DRUG": "Small molecule", "BIOLOGICAL": "Other biologic",
                          "DIETARY_SUPPLEMENT": "Supplement"}
# GENETIC and COMBINATION_PRODUCT deliberately have no fallback: they stay "Unresolved".

# Evidence tiers (strongest first)
TIER_CURATED_VERIFIED = "curated asset (source verified)"
TIER_CURATED_REGISTRY = "curated asset (registry text only)"
TIER_CURATED_UNVERIFIED = "curated rule carried over from v1 (source not re-verified)"
TIER_OVERRIDE = "manual override"
TIER_KW_NAME = "keyword in intervention name"
TIER_KW_DESC = "keyword in intervention description"
TIER_TYPE = "registry-type fallback"
TIER_NONE = "unresolved"
EXPLICIT_MOD_TIERS = [TIER_CURATED_VERIFIED, TIER_CURATED_REGISTRY, TIER_CURATED_UNVERIFIED, TIER_OVERRIDE,
                      TIER_KW_NAME, TIER_KW_DESC]

# Known-answer REGRESSION checks (not accuracy evidence).
KNOWN_TRIALS = [
    ("NCT02623699", "SOD1", "ASO", "VALOR - tofersen, SOD1-ALS"),
    ("NCT04856982", "SOD1", "ASO", "ATLAS - tofersen, presymptomatic SOD1 carriers"),
]
KNOWN_ASSETS = [
    ("tofersen", r"tofersen", "SOD1", "ASO"),
    ("ulefnersen (ION363)", r"ulefnersen|ION[- ]?363", "FUS", "ASO"),
    ("VTx-002", r"VTx[- ]?002", "TDP-43 biology", "Vectorized antibody"),
    ("CNM-Au8", r"CNM[- ]?Au[- ]?8|gold nanocrystal", "Mitochondria & cell stress", "Inorganic nanocrystal"),
    ("ARO-SOD1", r"ARO[- ]?SOD1", "SOD1", "siRNA / RNAi"),
]


# --------------------------------------------------------------------------------------
# HTTP helpers (only used with --refresh)
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
    raise RuntimeError(f"Request failed after {retries} attempts: {url}\n  last error: {last}")


def get_json(url, params=None):
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    return _request(url)


def post_json(url, payload):
    return _request(url, data=json.dumps(payload).encode("utf-8"), headers={"Content-Type": "application/json"})


def now_iso():
    return dt.datetime.now().replace(microsecond=0).isoformat()


def fetch_ctgov():
    print(f"[ClinicalTrials.gov] querying condition = '{CONDITION}' ...")
    params = {"query.cond": CONDITION, "pageSize": 1000, "countTotal": "true", "format": "json"}
    studies, total, token = [], None, None
    while True:
        p = dict(params)
        if token:
            p["pageToken"] = token
        resp = get_json(CTGOV_URL, p)
        total = resp.get("totalCount") if total is None else total
        studies.extend(resp.get("studies", []))
        token = resp.get("nextPageToken")
        if not token:
            break
        time.sleep(0.3)
    return {"fetched_at": now_iso(), "query": CONDITION, "query_params": params, "totalCount": total,
            "studies": studies}


OT_SEARCH_Q = """
query Search($q: String!) {
  search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 5}) { hits { id name entity } }
}"""
OT_ASSOC_Q = """
query Assoc($id: String!, $size: Int!) {
  disease(efoId: $id) { id name
    associatedTargets(page: {index: 0, size: $size}) {
      count rows { score target { id approvedSymbol approvedName } datatypeScores { id score } } } }
}"""


def fetch_opentargets():
    hits = post_json(OT_URL, {"query": OT_SEARCH_Q, "variables": {"q": CONDITION}})
    hits = (((hits.get("data") or {}).get("search") or {}).get("hits")) or []
    disease = next((h for h in hits if h.get("name", "").lower() == CONDITION), hits[0] if hits else None)
    if not disease:
        raise RuntimeError("Open Targets: disease not found")
    resp = post_json(OT_URL, {"query": OT_ASSOC_Q, "variables": {"id": disease["id"], "size": N_OT_TARGETS}})
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
    time.sleep(0.12 if os.environ.get("NCBI_API_KEY") else 0.4)
    return int(r["esearchresult"]["count"])


def fetch_pubmed(symbols, year):
    base = f'"{CONDITION}"[tiab]'
    out = {"fetched_at": now_iso(), "als_total": pubmed_count(base), "window": [year - 4, year], "counts": {}}
    for s in symbols:
        t = f"{base} AND {PUBMED_ALIASES.get(s, s + '[tiab]')}"
        out["counts"][s] = {"total": pubmed_count(t), "recent": pubmed_count(t, year - 4, year), "term": t}
    return out


def load_or_fetch(path, fetch_fn, refresh):
    if os.path.exists(path) and not refresh:
        with open(path) as f:
            return json.load(f)
    data = fetch_fn()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f)
    return data


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------------------
# Curation and overrides
# --------------------------------------------------------------------------------------
CURATION_COLS = ["asset_id", "rule_kind", "match_regex", "mechanism", "modality", "direct_targets",
                 "target_relationship", "source_status", "source_urls", "rationale"]


def load_curation(path=CURATION_PATH):
    rules = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != CURATION_COLS:
            raise ValueError(f"{path}: unexpected columns {reader.fieldnames}")
        for r in reader:
            if None in r:
                raise ValueError(f"{path}: row for {r.get('asset_id')} has too many fields")
            st = r["source_status"]
            tier = (TIER_CURATED_VERIFIED if st.startswith("verified") else
                    TIER_CURATED_REGISTRY if st == "registry_only" else TIER_CURATED_UNVERIFIED)
            mech = r["mechanism"] or None
            mod = r["modality"] or None
            if mech and mech not in MECHANISMS:
                raise ValueError(f"unknown mechanism {mech!r} for {r['asset_id']}")
            if mod and mod not in MODALITIES:
                raise ValueError(f"unknown modality {mod!r} for {r['asset_id']}")
            rules.append({**r, "rx": re.compile(r["match_regex"], re.I), "mech": mech, "mod": mod, "tier": tier,
                          "targets": [g.strip() for g in r["direct_targets"].split(";") if g.strip()]})
    return rules


def load_overrides(data_dir):
    """data/overrides.csv: intervention_name, mechanism, modality, source, rationale (optional)."""
    path = os.path.join(data_dir, "overrides.csv")
    ov = {}
    if os.path.exists(path):
        with open(path, newline="") as f:
            for r in csv.DictReader(f):
                name = (r.get("intervention_name") or "").strip().lower()
                mech = (r.get("mechanism") or "").strip() or None
                mod = (r.get("modality") or "").strip() or None
                if name and (mech or mod):
                    ov[name] = {"mech": mech, "mod": mod, "source": (r.get("source") or "").strip(),
                                "rationale": (r.get("rationale") or "").strip()}
    return ov


# --------------------------------------------------------------------------------------
# Dates
# --------------------------------------------------------------------------------------
def parse_registry_date(s):
    """Registry dates are YYYY-MM-DD or YYYY-MM (month precision). Keep the precision explicitly."""
    if not s:
        return None
    m = re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})", s)
    if m:
        d = dt.date(int(m[1]), int(m[2]), int(m[3]))
        return {"raw": s, "precision": "day", "start": d, "end": d}
    m = re.fullmatch(r"(\d{4})-(\d{2})", s)
    if m:
        y, mo = int(m[1]), int(m[2])
        return {"raw": s, "precision": "month", "start": dt.date(y, mo, 1),
                "end": dt.date(y, mo, calendar.monthrange(y, mo)[1])}
    m = re.fullmatch(r"(\d{4})", s)
    if m:
        y = int(m[1])
        return {"raw": s, "precision": "year", "start": dt.date(y, 1, 1), "end": dt.date(y, 12, 31)}
    return None


def milestone_position(d, analysis_date, horizon):
    """'passed' only if the WHOLE registered period (day or month) ended before the analysis date."""
    if not d:
        return "no date"
    if d["end"] < analysis_date:
        return "passed"
    if d["start"] > horizon:
        return "beyond window"
    return "in window"


# --------------------------------------------------------------------------------------
# Roles: what each intervention is doing in its trial
# --------------------------------------------------------------------------------------
def resolve_roles(interventions, arm_groups):
    arm_type = {a.get("label"): a.get("type") or "UNKNOWN_ARM" for a in (arm_groups or [])}
    trial_has_exp = "EXPERIMENTAL" in arm_type.values()
    out = []
    for iv in interventions:
        name = iv.get("name") or ""
        labels = iv.get("armGroupLabels") or []
        types = sorted({arm_type.get(lb, "UNKNOWN_ARM") for lb in labels})
        soc = bool(SOC_RE.search(name))
        supportive = bool(SUPPORTIVE_RE.search(name))
        flags = []
        if PLACEBO_RE.search(name) and not is_placebo_name(name):
            flags.append("name combines an active agent with placebo (blinded naming)")
        if is_placebo_name(name):
            role, basis = ROLE_PLACEBO, "name describes only a placebo/sham/vehicle"
        elif NONTHERAPEUTIC_RE.search(name):
            role, basis = ROLE_NONTHERAPEUTIC, "name matches sampling/imaging/enabling-agent pattern"
        elif not arm_type or not labels:
            role, basis = ROLE_EXPERIMENTAL_INFERRED, "no arm-group linkage in registry; treated as test intervention"
            flags.append("no arm data")
        else:
            ts = set(types)
            exp = "EXPERIMENTAL" in ts
            ctrl = bool(ts & {"PLACEBO_COMPARATOR", "SHAM_COMPARATOR", "NO_INTERVENTION"})
            act = "ACTIVE_COMPARATOR" in ts
            other = bool(ts & {"OTHER", "UNKNOWN_ARM"})
            if exp and (ctrl or act):
                if soc:
                    role, basis = ROLE_BACKGROUND, "standard-of-care agent linked to experimental and control/comparator arms"
                elif ctrl and not act:
                    role, basis = ROLE_EXPERIMENTAL, "linked to experimental arm (also linked to placebo/sham arm)"
                    flags.append("also linked to placebo/sham arm")
                else:
                    role, basis = ROLE_UNRESOLVED, "linked to both experimental and active-comparator arms"
            elif exp:
                role, basis = ROLE_EXPERIMENTAL, "linked only to arm(s) typed EXPERIMENTAL"
            elif act and not ctrl:
                if soc or trial_has_exp:
                    role, basis = ROLE_COMPARATOR, "linked only to ACTIVE_COMPARATOR arm(s)"
                else:
                    role, basis = ROLE_EXPERIMENTAL_INFERRED, "no arm typed EXPERIMENTAL; active-comparator arm treated as test arm"
                    flags.append("no EXPERIMENTAL arm in registry")
            elif ctrl and not act:
                role, basis = ROLE_CONTROL, "linked only to placebo/sham/no-intervention arm(s)"
            elif act and ctrl:
                if soc:
                    role, basis = ROLE_BACKGROUND, "standard-of-care agent linked to comparator and placebo arms"
                elif trial_has_exp:
                    role, basis = ROLE_UNRESOLVED, "linked to comparator and placebo arms"
                else:
                    role, basis = ROLE_EXPERIMENTAL_INFERRED, "no arm typed EXPERIMENTAL; linked to comparator and placebo arms"
                    flags.append("no EXPERIMENTAL arm in registry")
            elif other:
                if trial_has_exp:
                    role, basis = ROLE_UNRESOLVED, "linked only to arm typed OTHER"
                else:
                    role, basis = ROLE_EXPERIMENTAL_INFERRED, "no arm typed EXPERIMENTAL; arm typed OTHER treated as test arm"
                    flags.append("no EXPERIMENTAL arm in registry")
            else:
                role, basis = ROLE_UNRESOLVED, f"arm types {types}"
        out.append({"role": role, "role_basis": basis, "arm_labels": labels, "arm_types": types,
                    "soc_name": soc, "supportive_name": supportive, "flags": flags})
    # Post-pass: a standard-of-care agent co-listed with an investigational drug in the same arms
    # (or with no arm data) is background therapy, not the mechanism being tested.
    drug_idx = [i for i, iv in enumerate(interventions) if iv.get("type") in DRUG_TYPES]
    for i in drug_idx:
        r = out[i]
        if not (r["soc_name"] or r["supportive_name"]) or r["role"] not in (ROLE_EXPERIMENTAL, ROLE_EXPERIMENTAL_INFERRED):
            continue
        for j in drug_idx:
            if j == i or out[j]["soc_name"] or out[j]["supportive_name"] or out[j]["role"] not in TEST_ROLES:
                continue
            if not r["arm_labels"] or not out[j]["arm_labels"] or set(r["arm_labels"]) & set(out[j]["arm_labels"]):
                r["role"] = ROLE_BACKGROUND
                r["role_basis"] = ("supportive/prophylactic co-medication given alongside an investigational agent"
                                   if r["supportive_name"] and not r["soc_name"] else
                                   "standard-of-care agent co-listed with an investigational agent in the same arm(s)")
                r["flags"].append("background inferred from co-listing")
                break
    return out


# --------------------------------------------------------------------------------------
# Intervention classification
# --------------------------------------------------------------------------------------
def split_background_components(name):
    """'Plasma ... + Riluzole' -> ('Plasma ...', ['Riluzole']). Only standard-of-care parts are split off."""
    parts = [p for p in COMBO_SPLIT_RE.split(name) if p and p.strip()]
    if len(parts) < 2:
        return name, []
    soc = [p.strip() for p in parts if SOC_RE.search(p)]
    rest = [p.strip() for p in parts if not SOC_RE.search(p)]
    if soc and rest:
        return " + ".join(rest), soc
    return name, []


def classify_intervention(iv, curation, overrides):
    raw_name = iv.get("name") or ""
    core_name, bg_parts = split_background_components(raw_name)
    other_names = [n for n in (iv.get("otherNames") or []) if not SOC_RE.search(n)]
    name_text = " | ".join([core_name] + other_names)
    desc = iv.get("description") or ""
    reg_type = iv.get("type") or ""
    res = {"name": raw_name, "classified_text": core_name, "registry_type": reg_type,
           "background_components": bg_parts,
           "mechanism": None, "mechanism_tier": None, "mechanism_rule": None,
           "modality": None, "modality_tier": None, "modality_rule": None,
           "asset_id": None, "direct_targets": [], "target_relationship": "",
           "source_urls": "", "mechanism_candidate": None, "candidate_rule": None, "flags": []}
    if bg_parts:
        res["flags"].append("combined registry intervention; standard-of-care component(s) treated as background: "
                            + ", ".join(bg_parts))

    ov = overrides.get(raw_name.strip().lower())
    if ov:
        if ov["mech"]:
            res.update(mechanism=ov["mech"], mechanism_tier=TIER_OVERRIDE, mechanism_rule="override")
        if ov["mod"]:
            res.update(modality=ov["mod"], modality_tier=TIER_OVERRIDE, modality_rule="override")
        res["source_urls"] = ov["source"]

    for r in curation:
        if r["rx"].search(name_text):
            res["asset_id"] = r["asset_id"]
            res["asset_kind"] = r["rule_kind"]
            res["direct_targets"] = r["targets"]
            res["target_relationship"] = r["target_relationship"]
            res["source_urls"] = res["source_urls"] or r["source_urls"]
            if r["mech"] and not res["mechanism"]:
                res.update(mechanism=r["mech"], mechanism_tier=r["tier"], mechanism_rule="asset:" + r["asset_id"])
            if r["mod"] and not res["modality"]:
                res.update(modality=r["mod"], modality_tier=r["tier"], modality_rule="asset:" + r["asset_id"])
            break
    else:
        for r in curation:
            if r["rx"].search(desc) and r["asset_id"] not in ("riluzole", "edaravone"):
                res["mechanism_candidate"] = r["mech"]
                res["candidate_rule"] = "asset:" + r["asset_id"] + " (description only)"
                res["flags"].append(f"curated asset '{r['asset_id']}' mentioned only in description")
                break

    for rid, rx, mech, mod in _KW:
        if mech and not res["mechanism"] and rx.search(name_text):
            res.update(mechanism=mech, mechanism_tier=TIER_KW_NAME, mechanism_rule=rid)
        if mod and not res["modality"] and rx.search(name_text):
            res.update(modality=mod, modality_tier=TIER_KW_NAME, modality_rule=rid)
    if not res["modality"]:
        for rid, rx, mech, mod in _KW:
            if mod and rx.search(desc):
                res.update(modality=mod, modality_tier=TIER_KW_DESC, modality_rule=rid)
                break
    if not res["mechanism"] and not res["mechanism_candidate"]:
        for rid, rx, mech, mod in _KW:
            m = rx.search(desc) if mech else None
            if m:
                res["mechanism_candidate"] = mech
                res["candidate_rule"] = f"{rid} (description: '{m.group(0)}')"
                break

    if not res["modality"]:
        fb = TYPE_FALLBACK_MODALITY.get(reg_type)
        if fb:
            res.update(modality=fb, modality_tier=TIER_TYPE, modality_rule=f"type:{reg_type}")
        else:
            res.update(modality="Unresolved", modality_tier=TIER_NONE, modality_rule=f"type:{reg_type} (no fallback)")
    if not res["mechanism"]:
        res.update(mechanism="Unclassified", mechanism_tier=TIER_NONE, mechanism_rule="")
    return res


# --------------------------------------------------------------------------------------
# Trials
# --------------------------------------------------------------------------------------
PHASE_RANK = {"EARLY_PHASE1": 0.5, "PHASE1": 1, "PHASE2": 2, "PHASE3": 3, "PHASE4": 4}


def phase_label(phases):
    ph = [p for p in (phases or []) if p in PHASE_RANK]
    if not ph:
        return "N/A", 0
    ph.sort(key=PHASE_RANK.get)
    nums = [("Early 1" if p == "EARLY_PHASE1" else p.replace("PHASE", "")) for p in ph]
    return "Ph " + "/".join(nums), PHASE_RANK[ph[-1]]


def _sorted_unique(items, order):
    return sorted(set(items), key=lambda x: (order.index(x) if x in order else len(order), x))


def normalise_trials(ct, curation, overrides):
    trials = []
    for s in ct["studies"]:
        ps = s.get("protocolSection", {})
        idm, stm = ps.get("identificationModule", {}), ps.get("statusModule", {})
        dsm, spm = ps.get("designModule", {}), ps.get("sponsorCollaboratorsModule", {})
        aim, clm = ps.get("armsInterventionsModule", {}), ps.get("contactsLocationsModule", {})
        if dsm.get("studyType") != "INTERVENTIONAL":
            continue
        nct = idm.get("nctId", "")
        ivs = aim.get("interventions") or []
        arms = aim.get("armGroups") or []
        roles = resolve_roles(ivs, arms)
        records = []
        for k, (iv, rl) in enumerate(zip(ivs, roles)):
            rec = {"nct_id": nct, "iv_index": k, "name": iv.get("name") or "", "registry_type": iv.get("type") or "",
                   "drug_type": iv.get("type") in DRUG_TYPES, **rl}
            if rec["drug_type"] and rl["role"] in TEST_ROLES:
                c = classify_intervention(iv, curation, overrides)
                rec.update(c)
                rec["flags"] = rl["flags"] + c["flags"]
            records.append(rec)
        test = [r for r in records if r["drug_type"] and r["role"] in TEST_ROLES]
        background = [r["name"] for r in records if r["drug_type"] and r["role"] == ROLE_BACKGROUND]
        background += [b for r in test for b in r.get("background_components", [])]
        status = stm.get("overallStatus", "UNKNOWN")
        plab, prank = phase_label(dsm.get("phases"))
        pcd_struct = stm.get("primaryCompletionDateStruct") or {}
        mechs = _sorted_unique([r["mechanism"] for r in test], MECHANISMS)
        classified_mechs = [m for m in mechs if m != "Unclassified"]
        n_unc = sum(r["mechanism"] == "Unclassified" for r in test)
        trials.append({
            "nct_id": nct,
            "title": idm.get("briefTitle", ""),
            "status": status,
            "status_group": ("Active" if status in ACTIVE_STATUSES else "Completed" if status == "COMPLETED"
                             else "Stopped" if status in STOPPED_STATUSES else "Other/unknown"),
            "phase": plab, "phase_rank": prank,
            "enrollment": (dsm.get("enrollmentInfo") or {}).get("count"),
            "sponsor": (spm.get("leadSponsor") or {}).get("name", ""),
            "sponsor_class": (spm.get("leadSponsor") or {}).get("class", "UNKNOWN"),
            "start_raw": (stm.get("startDateStruct") or {}).get("date") or "",
            "pcd": parse_registry_date(pcd_struct.get("date")),
            "pcd_raw": pcd_struct.get("date") or "",
            "pcd_type": pcd_struct.get("type") or "",
            "last_update_posted": (stm.get("lastUpdatePostDateStruct") or {}).get("date") or "",
            "status_verified": stm.get("statusVerifiedDate") or "",
            "has_results": bool(s.get("hasResults")),
            "countries": sorted({loc.get("country") for loc in (clm.get("locations") or []) if loc.get("country")}),
            "records": records,
            "test": test,
            "has_drug_type_any_role": any(r["drug_type"] for r in records),
            "is_drug": bool(test),
            "background": sorted(set(background)),
            "comparators": sorted({r["name"] for r in records if r["drug_type"] and r["role"] in (ROLE_COMPARATOR, ROLE_CONTROL)}),
            "mechanisms": mechs,
            "classified_mechanisms": classified_mechs,
            "mech_status": ("none" if not test else "all classified" if n_unc == 0
                            else "none classified" if n_unc == len(test) else "partially classified"),
            "modalities": _sorted_unique([r["modality"] for r in test], MODALITIES),
            "pairs": sorted({(r["mechanism"], r["modality"]) for r in test if r["mechanism"] != "Unclassified"}),
            "assets": sorted({r["asset_id"] for r in test if r.get("asset_id") and r.get("asset_kind") != "drug_class"}),
            "direct_targets": sorted({g for r in test for g in r.get("direct_targets", [])}),
        })
    trials.sort(key=lambda t: t["nct_id"])
    return trials


def sponsor_group(c):
    return {"INDUSTRY": "Industry", "NIH": "Government", "FED": "Government", "OTHER_GOV": "Government",
            "NETWORK": "Network", "OTHER": "Academic / hospital", "INDIV": "Individual"}.get(c, "Other")


def most_common_sorted(counter, n=None):
    items = sorted(counter.items(), key=lambda kv: (-kv[1], str(kv[0])))
    return items[:n] if n else items


# --------------------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------------------
def analyse(ot, pm, trials, curation, analysis_date):
    horizon = dt.date(analysis_date.year + MILESTONE_HORIZON_YEARS, 12, 31)
    drug = [t for t in trials if t["is_drug"]]
    active = [t for t in drug if t["status_group"] == "Active"]
    active_classified = [t for t in active if t["classified_mechanisms"]]

    # --- coverage ------------------------------------------------------------------
    def iv_cov(ts):
        ivs = [r for t in ts for r in t["test"]]
        mech_t = collections.Counter(r["mechanism_tier"] for r in ivs)
        mod_t = collections.Counter(r["modality_tier"] for r in ivs)
        return {"n": len(ivs), "mech_tiers": mech_t, "mod_tiers": mod_t,
                "mech_classified": sum(r["mechanism"] != "Unclassified" for r in ivs),
                "mech_candidates": sum(1 for r in ivs if r["mechanism"] == "Unclassified" and r["mechanism_candidate"]),
                "mod_explicit": sum(mod_t[x] for x in EXPLICIT_MOD_TIERS),
                "mod_fallback": mod_t[TIER_TYPE], "mod_unresolved": mod_t[TIER_NONE]}

    def trial_cov(ts):
        c = collections.Counter(t["mech_status"] for t in ts)
        return {"n": len(ts), "all": c["all classified"], "partial": c["partially classified"],
                "none": c["none classified"], "any": c["all classified"] + c["partially classified"]}

    role_counts = collections.Counter(r["role"] for t in trials for r in t["records"] if r["drug_type"])
    coverage = {
        "iv_all": iv_cov(drug), "iv_active": iv_cov(active),
        "trial_all": trial_cov(drug), "trial_active": trial_cov(active),
        "roles": role_counts,
        "n_interventional": len(trials),
        "n_drug_any_role": sum(t["has_drug_type_any_role"] for t in trials),
        "n_drug": len(drug), "n_active": len(active),
        "unique_assets_all": sorted({a for t in drug for a in t["assets"]}),
        "unique_assets_active": sorted({a for t in active for a in t["assets"]}),
    }

    # --- classified activity grid (unique trials per cell) ---------------------------------
    grid_active, grid_all = collections.Counter(), collections.Counter()
    grid_ids = collections.defaultdict(list)
    grid_fallback = collections.Counter()
    for t in drug:
        fb_pairs = {(r["mechanism"], r["modality"]) for r in t["test"] if r["modality_tier"] in (TIER_TYPE, TIER_NONE)}
        for pair in t["pairs"]:
            grid_all[pair] += 1
            if t["status_group"] == "Active":
                grid_active[pair] += 1
                grid_ids[pair].append(t["nct_id"])
                if pair in fb_pairs:
                    grid_fallback[pair] += 1
    mech_active = collections.Counter(m for t in active for m in t["classified_mechanisms"])
    mech_all = collections.Counter(m for t in drug for m in t["classified_mechanisms"])
    multi_mech_active = [t for t in active if len(t["classified_mechanisms"]) > 1]
    multi_cell_active = [t for t in active if len(t["pairs"]) > 1]

    # --- phase profile (classified subset; unclassified reported separately) ------------
    phase_cols = ["Ph 1", "Ph 2", "Ph 3", "Ph 4", "N/A"]
    colmap = {0: "N/A", 0.5: "Ph 1", 1: "Ph 1", 2: "Ph 2", 3: "Ph 3", 4: "Ph 4"}
    phase_prof = collections.defaultdict(collections.Counter)
    phase_unc = collections.Counter()
    for t in active:
        col = colmap[t["phase_rank"]]
        for m in t["classified_mechanisms"]:
            phase_prof[m][col] += 1
        if not t["classified_mechanisms"]:
            phase_unc[col] += 1

    # --- registered primary-completion milestones ---------------------------------------
    ms_pool = [t for t in active if t["phase_rank"] >= 2]
    for t in ms_pool:
        t["_ms_pos"] = milestone_position(t["pcd"], analysis_date, horizon)
    milestones = sorted([t for t in ms_pool if t["_ms_pos"] == "in window"],
                        key=lambda t: (t["pcd"]["start"], t["nct_id"]))
    passed = sorted([t for t in ms_pool if t["_ms_pos"] == "passed"], key=lambda t: (t["pcd"]["start"], t["nct_id"]))
    ms_summary = {
        "pool": len(ms_pool), "in_window": len(milestones), "passed": len(passed),
        "beyond": sum(t["_ms_pos"] == "beyond window" for t in ms_pool),
        "no_date": sum(t["_ms_pos"] == "no date" for t in ms_pool),
        "in_window_month_precision": sum(t["pcd"]["precision"] == "month" for t in milestones),
        "in_window_current_month": sum(t["pcd"]["start"] <= analysis_date <= t["pcd"]["end"] for t in milestones),
        "passed_actual": sum(t["pcd_type"] == "ACTUAL" for t in passed),
        "passed_estimated": sum(t["pcd_type"] == "ESTIMATED" for t in passed),
        "phase12_in_window": sum(t["phase"] == "Ph 1/2" for t in milestones),
    }

    # --- target evidence vs curated clinical mapping -----------------------------------------
    assessed, display = {}, {}
    for r in curation:
        for g in r["targets"]:
            assessed.setdefault(g.upper(), []).append(r["asset_id"])
            display[g.upper()] = g
    targets = []
    ot_rows = ot.get("rows") or []
    pm_counts = (pm or {}).get("counts")
    for row in ot_rows:
        sym = row["target"]["approvedSymbol"]
        dts = {d["id"]: d["score"] for d in (row.get("datatypeScores") or [])}
        pmc = (pm_counts or {}).get(sym)
        rec = {"symbol": sym, "name": row["target"].get("approvedName", ""), "ensembl": row["target"].get("id", ""),
               "ot_score": row.get("score"), "genetic_score": dts.get("genetic_association"),
               "literature_score": dts.get("literature"),
               "pubmed_total": pmc.get("total") if pmc else None, "pubmed_recent": pmc.get("recent") if pmc else None,
               "pubmed_status": ("unavailable" if pm_counts is None else "retrieved" if pmc else "not retrieved"),
               "curated_assets": "", "direct_active": None, "direct_all": None, "direct_ids_active": [],
               "bucket": "", "bucket_active": None, "bucket_all": None, "flag": ""}
        key = sym.upper()
        if key in assessed:
            hit = [t for t in drug if key in {g.upper() for g in t["direct_targets"]}]
            rec["curated_assets"] = "; ".join(sorted(set(assessed[key])))
            rec["direct_all"] = len(hit)
            rec["direct_active"] = sum(t["status_group"] == "Active" for t in hit)
            rec["direct_ids_active"] = [t["nct_id"] for t in hit if t["status_group"] == "Active"]
            rec["mapping_state"] = "mapped: matching trials identified" if hit else \
                "mapped: no matching trials identified in searched dataset"
            if not hit and (rec["genetic_score"] or 0) >= 0.5:
                rec["flag"] = "review: no trial found via curated assets"
        else:
            rec["mapping_state"] = "mapping not implemented"
        b = BUCKET_FOR_GENE.get(key)
        if b and key in assessed:
            rec["bucket"] = b
            rec["bucket_active"] = mech_active.get(b, 0)
            rec["bucket_all"] = mech_all.get(b, 0)
        targets.append(rec)
    off_list = sorted(display[g] for g in set(assessed) - {r["target"]["approvedSymbol"].upper() for r in ot_rows})

    # --- sponsors & geography --------------------------------------------------------
    spons = collections.Counter(sponsor_group(t["sponsor_class"]) for t in active)
    countries = collections.Counter(c for t in active for c in t["countries"])
    asia_trials = [t for t in active if any(c in ASIA for c in t["countries"])]
    no_sites = sum(1 for t in active if not t["countries"])
    top_sponsors = collections.Counter(t["sponsor"] for t in active if t["sponsor_class"] == "INDUSTRY")

    tdp = [t for t in drug if "TDP-43 biology" in t["classified_mechanisms"]]

    review_queue = collections.Counter()
    review_meta = {}
    for t in drug:
        for r in t["test"]:
            if r["mechanism"] == "Unclassified" or r["modality_tier"] in (TIER_NONE,):
                key = r["name"].strip()
                review_queue[key] += 1
                review_meta.setdefault(key, (r["registry_type"], r["modality"], r["modality_tier"],
                                             r.get("mechanism_candidate") or "", r.get("candidate_rule") or ""))

    return {"analysis_date": analysis_date, "horizon": horizon, "drug": drug, "active": active,
            "active_classified": active_classified, "coverage": coverage,
            "grid_active": grid_active, "grid_all": grid_all, "grid_ids": grid_ids, "grid_fallback": grid_fallback,
            "mech_active": mech_active, "mech_all": mech_all,
            "multi_mech_active": multi_mech_active, "multi_cell_active": multi_cell_active,
            "phase_cols": phase_cols, "phase_prof": phase_prof, "phase_unc": phase_unc,
            "milestones": milestones, "passed": passed, "ms_summary": ms_summary,
            "targets": targets, "assessed_genes": sorted(display.values()), "assessed_off_list": off_list,
            "sponsors": spons, "countries": countries, "asia_trials": asia_trials, "no_sites": no_sites,
            "top_sponsors": top_sponsors, "tdp": tdp,
            "review_queue": review_queue, "review_meta": review_meta}


# --------------------------------------------------------------------------------------
# Software / data checks (NOT accuracy)
# --------------------------------------------------------------------------------------
def integrity_checks(ct, trials, a):
    checks = []

    def add(name, status, detail):
        checks.append({"name": name, "status": status, "detail": detail})

    n_fetched, n_total = len(ct["studies"]), ct.get("totalCount")
    add("Record count reconciles with API total", "PASS" if n_total == n_fetched else "FAIL",
        f"{n_fetched} studies in snapshot; API reported totalCount = {n_total}")
    ids = [t["nct_id"] for t in trials]
    bad = [i for i in ids if not re.fullmatch(r"NCT\d{8}", i)]
    dup = [i for i, c in collections.Counter(ids).items() if c > 1]
    add("Trial IDs well-formed and unique", "PASS" if not bad and not dup else "FAIL",
        f"{len(ids)} interventional trials; {len(bad)} malformed, {len(dup)} duplicated")
    bad_dates = sum(1 for t in a["active"] if t["pcd_raw"] and not t["pcd"])
    add("Registered primary-completion dates parse (with precision kept)", "PASS" if not bad_dates else "WARN",
        f"{len(a['active']) - bad_dates}/{len(a['active'])} active drug-type trials parse; "
        f"{sum(1 for t in a['active'] if t['pcd'] and t['pcd']['precision'] == 'month')} have month-only precision")
    leaked = [t["nct_id"] for t in trials for r in t["test"] if r["role"] not in TEST_ROLES]
    add("Only test-role interventions drive trial labels", "PASS" if not leaked else "FAIL",
        f"{len(leaked)} non-test interventions found in trial labels")
    flagged_unmapped = [t["symbol"] for t in a["targets"] if t["flag"] and t["mapping_state"] == "mapping not implemented"]
    add("Unmapped targets carry no review/opportunity flag", "PASS" if not flagged_unmapped else "FAIL",
        f"{len(flagged_unmapped)} unmapped targets flagged")
    return checks


def regression_answer_checks(trials):
    """Known-answer checks: useful for catching regressions, NOT evidence of overall accuracy."""
    out = []
    by_id = {t["nct_id"]: t for t in trials}
    for nct, mech, mod, label in KNOWN_TRIALS:
        t = by_id.get(nct)
        if not t:
            out.append({"name": f"Known trial: {label}", "status": "WARN", "detail": f"{nct} not in snapshot"})
            continue
        ok = (mech, mod) in t["pairs"]
        out.append({"name": f"Known trial: {label}", "status": "PASS" if ok else "FAIL",
                    "detail": f"{nct}: expected {mech} / {mod}; got "
                              + (", ".join(f"{m} / {d}" for m, d in t["pairs"]) or "none")})
    ivs = [r for t in trials for r in t["test"]]
    for label, rx, mech, mod in KNOWN_ASSETS:
        hits = [r for r in ivs if re.search(rx, r["name"], re.I)]
        if not hits:
            out.append({"name": f"Known asset: {label}", "status": "WARN", "detail": "not in snapshot"})
            continue
        wrong = [r for r in hits if (r["mechanism"], r["modality"]) != (mech, mod)]
        out.append({"name": f"Known asset: {label}", "status": "PASS" if not wrong else "FAIL",
                    "detail": f"{len(hits)} test intervention record(s); expected {mech} / {mod}; {len(wrong)} mismatched"})
    t = by_id.get("NCT04454840")
    if t:
        ok = (t["test"] and all(r["modality"] != "Small molecule" and r["mechanism"] != "Excitotoxicity"
                                for r in t["test"]) and any("Riluzole" in b for b in t["background"]))
        out.append({"name": "Plasma + riluzole trial (NCT04454840) not labelled as riluzole",
                    "status": "PASS" if ok else "FAIL",
                    "detail": "test: " + "; ".join(f"{r['classified_text']} -> {r['mechanism']} / {r['modality']}"
                                                   for r in t["test"])
                              + " | background: " + ("; ".join(t["background"]) or "none")
                              + " | comparator/control: " + ("; ".join(t["comparators"]) or "none")})
    return out


# --------------------------------------------------------------------------------------
# Validation status (human-approved labels only)
# --------------------------------------------------------------------------------------
def fingerprint(trials):
    h = hashlib.sha256()
    for t in trials:
        h.update(f"{t['nct_id']}|{';'.join(t['classified_mechanisms'])}|{';'.join(t['modalities'])}\n".encode())
    return h.hexdigest()[:16]


def _norm_set(s):
    """Label set; 'N/A' (no therapeutic intervention) means an EMPTY expected set."""
    out = frozenset(x.strip().lower() for x in (s or "").split(";") if x.strip())
    return frozenset() if out == {"n/a"} else out


def evaluate_pack(pack_rows, trials, split):
    """Agreement against HUMAN-approved labels only, per axis, with separate denominators.

    A row counts for an axis only if the reviewer approved it (status approved/corrected) AND entered
    a label for that axis. Rows reviewed for one axis are never counted as errors on the other.
    """
    by_id = {t["nct_id"]: t for t in trials}
    res = {}
    for subset in ("heldout", "dev"):
        out = {}
        for axis, col, pred_key in (("mechanism", "human_mechanism", "mechanisms"),
                                    ("modality", "human_modality", "modalities")):
            n = agree = 0
            disagreements = []
            for r in pack_rows:
                if split.get(r["nct_id"]) != subset:
                    continue
                if (r.get("human_review_status") or "").strip().lower() not in ("approved", "corrected"):
                    continue
                if not (r.get(col) or "").strip():
                    continue  # this axis not reviewed for this row: never counted, never an error
                gold = _norm_set(r.get(col))
                t = by_id.get(r["nct_id"])
                pred = _norm_set(";".join(t[pred_key])) if t else frozenset()
                n += 1
                if pred == gold:
                    agree += 1
                else:
                    disagreements.append({"nct_id": r["nct_id"], "predicted": sorted(pred), "human": sorted(gold)})
            out[axis] = {"n": n, "agree": agree, "rate": (agree / n) if n else None, "disagreements": disagreements}
        res[subset] = out
    return res


def load_split(valid_dir):
    path = os.path.join(valid_dir, "eval_split.csv")
    if not os.path.exists(path):
        return {}
    with open(path, newline="") as f:
        return {r["nct_id"]: r["split"] for r in csv.DictReader(f)}


def validation_status(valid_dir, trials):
    """Return a dict describing what has (and has not) been validated, and whether it is current."""
    path = os.path.join(valid_dir, "validation_results.json")
    fp = fingerprint(trials)
    if not os.path.exists(path):
        return {"state": "not measured", "headline": "Classification accuracy not yet measured",
                "detail": "No human-approved reference labels exist yet. The review pack contains AI-proposed "
                          "labels only; they are not a reference standard.", "results": None}
    with open(path) as f:
        res = json.load(f)
    stale = (res.get("classification_version") != CLASSIFICATION_VERSION or res.get("trials_fingerprint") != fp)
    held = res["results"]["heldout"]
    n_any = held["mechanism"]["n"] + held["modality"]["n"]
    if stale:
        return {"state": "stale", "headline": "Classification accuracy not yet measured for this version",
                "detail": f"validation_results.json belongs to classifier {res.get('classification_version')} / "
                          f"fingerprint {res.get('trials_fingerprint')}; current is {CLASSIFICATION_VERSION} / {fp}. "
                          "Re-run --evaluate.", "results": res}
    if n_any == 0:
        return {"state": "not measured", "headline": "Classification accuracy not yet measured",
                "detail": "No human-approved labels in the held-out subset.", "results": res}
    return {"state": "measured", "headline": "Agreement with human-reviewed labels (held-out subset)",
            "detail": "", "results": res}


def cmd_evaluate(data_dir, trials):
    valid_dir = os.path.join(data_dir, "validation")
    pack_path = os.path.join(valid_dir, "label_review_pack.csv")
    if not os.path.exists(pack_path):
        sys.exit("No review pack yet - run the pipeline first.")
    with open(pack_path, newline="") as f:
        rows = list(csv.DictReader(f))
    split = load_split(valid_dir)
    res = evaluate_pack(rows, trials, split)
    out = {"classification_version": CLASSIFICATION_VERSION, "trials_fingerprint": fingerprint(trials),
           "review_pack_sha256": sha256_file(pack_path), "evaluated_at": now_iso(), "results": res,
           "note": "Only rows with human_review_status approved/corrected and a non-empty human label for the axis "
                   "are scored. Held-out rows must not be used to change rules or overrides; dev-subset agreement "
                   "is NOT independent validation."}
    with open(os.path.join(valid_dir, "validation_results.json"), "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
    for subset in ("heldout", "dev"):
        for axis in ("mechanism", "modality"):
            r = res[subset][axis]
            rate = "n/a" if r["rate"] is None else f"{r['rate']:.0%}"
            print(f"  {subset:8s} {axis:9s}: {r['agree']}/{r['n']} agree ({rate})")
    print("Re-run the pipeline to show the result on the report.")


# --------------------------------------------------------------------------------------
# Review pack (AI-proposed labels kept separate from human-approved labels)
# --------------------------------------------------------------------------------------
PACK_COLS = ["nct_id", "split", "registry_url", "title", "status_in_snapshot", "phase",
             "arms_context", "interventions_with_roles", "test_interventions",
             "pred_mechanisms", "pred_modalities", "pred_rules_and_tiers", "pred_flags",
             "proposed_mechanism", "proposed_modality", "proposed_confidence", "proposed_reasoning",
             "proposed_source_urls", "proposed_by", "open_questions", "triage_pred_vs_proposal",
             "human_mechanism", "human_modality", "human_review_status", "human_reviewer", "human_review_date",
             "human_notes"]
HUMAN_COLS = ["human_mechanism", "human_modality", "human_review_status", "human_reviewer", "human_review_date",
              "human_notes"]


def build_review_pack(data_dir, ct, trials):
    valid_dir = os.path.join(data_dir, "validation")
    sample_path = os.path.join(data_dir, "label_sample.csv")
    if not os.path.exists(sample_path):
        return None
    with open(sample_path, newline="") as f:
        sample_ids = [r["nct_id"] for r in csv.DictReader(f)]
    split = load_split(valid_dir)
    proposals = {}
    ppath = os.path.join(valid_dir, "proposed_labels.csv")
    if os.path.exists(ppath):
        with open(ppath, newline="") as f:
            proposals = {r["nct_id"]: r for r in csv.DictReader(f)}
    pack_path = os.path.join(valid_dir, "label_review_pack.csv")
    existing = {}
    if os.path.exists(pack_path):  # never lose human review work
        with open(pack_path, newline="") as f:
            existing = {r["nct_id"]: r for r in csv.DictReader(f)}
    studies = {s["protocolSection"]["identificationModule"]["nctId"]: s for s in ct["studies"]}
    by_id = {t["nct_id"]: t for t in trials}
    rows = []
    for nct in sample_ids:
        t = by_id.get(nct)
        s = studies.get(nct, {})
        arms = (s.get("protocolSection", {}).get("armsInterventionsModule", {}) or {}).get("armGroups") or []
        p = proposals.get(nct, {})
        row = {
            "nct_id": nct, "split": split.get(nct, ""), "registry_url": f"https://clinicaltrials.gov/study/{nct}",
            "title": t["title"] if t else "", "status_in_snapshot": t["status"] if t else "",
            "phase": t["phase"] if t else "",
            "arms_context": " || ".join(f"[{a.get('type')}] {a.get('label')}: {', '.join(a.get('interventionNames') or [])}"
                                        for a in arms),
            "interventions_with_roles": " || ".join(f"{r['name']} [{r['registry_type']}] -> {r['role']}"
                                                    for r in (t["records"] if t else [])),
            "test_interventions": " || ".join(r["name"] for r in (t["test"] if t else [])),
            "pred_mechanisms": "; ".join(t["mechanisms"]) if t else "",
            "pred_modalities": "; ".join(t["modalities"]) if t else "",
            "pred_rules_and_tiers": " || ".join(
                f"{r['name']}: mech={r['mechanism']} ({r['mechanism_tier']}; {r['mechanism_rule'] or '-'}), "
                f"mod={r['modality']} ({r['modality_tier']}; {r['modality_rule']})"
                + (f", candidate={r['mechanism_candidate']} via {r['candidate_rule']}" if r.get("mechanism_candidate") else "")
                for r in (t["test"] if t else [])),
            "pred_flags": " || ".join(f for r in (t["test"] if t else []) for f in r["flags"]),
            "proposed_mechanism": p.get("proposed_mechanism", ""), "proposed_modality": p.get("proposed_modality", ""),
            "proposed_confidence": p.get("proposed_confidence", ""), "proposed_reasoning": p.get("proposed_reasoning", ""),
            "proposed_source_urls": p.get("proposed_source_urls", ""),
            "proposed_by": p.get("proposed_by", ""), "open_questions": p.get("open_questions", ""),
        }
        if p:
            mm = [ax for ax, pk, pc in (("mechanism", "pred_mechanisms", "proposed_mechanism"),
                                         ("modality", "pred_modalities", "proposed_modality"))
                  if _norm_set(row[pk]) != _norm_set(row[pc])]
            row["triage_pred_vs_proposal"] = ("differs on " + " and ".join(mm) if mm else "same") + \
                " (triage aid only; AI proposal is not a reference label)"
        else:
            row["triage_pred_vs_proposal"] = ""
        prev = existing.get(nct, {})
        for c in HUMAN_COLS:
            row[c] = prev.get(c, "")
        rows.append(row)
    os.makedirs(valid_dir, exist_ok=True)
    with open(pack_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=PACK_COLS)
        w.writeheader()
        w.writerows(rows)
    return rows


# --------------------------------------------------------------------------------------
# Observations (deterministic; every statement carries its denominator)
# --------------------------------------------------------------------------------------
def pct(a, b):
    return f"{a / b:.0%}" if b else "n/a"


def observations(a):
    obs = []
    cov = a["coverage"]
    n_act = cov["n_active"]
    n_cls = len(a["active_classified"])
    n_unc = n_act - n_cls
    obs.append(f"<b>Coverage first:</b> a mechanism could be assigned for <b>{n_cls}</b> of <b>{n_act}</b> active "
               f"drug-type trials ({pct(n_cls, n_act)}). The other <b>{n_unc}</b> ({pct(n_unc, n_act)}) have no classified "
               "mechanism, so statements below describe the classified subset only, not the whole field.")
    ma = most_common_sorted(a["mech_active"], 3)
    if ma:
        def assets_in(m):
            return len({r["asset_id"] for t in a["active"] for r in t["test"]
                        if r["mechanism"] == m and r.get("asset_id") and r.get("asset_kind") != "drug_class"})
        obs.append("Most frequent mechanisms within the classified subset (unique active trials; a trial can count in "
                   "more than one): " + ", ".join(f"<b>{esc(m)}</b> ({c} trials; {assets_in(m)} distinct curated asset{'s' if assets_in(m) != 1 else ''} or compound{'s' if assets_in(m) != 1 else ''})"
                                                  for m, c in ma)
                   + ". Trial counts are not programme counts: several individualised (n-of-1) or extension studies "
                     "can share one asset.")
    gd = [t for t in a["active"] if set(t["classified_mechanisms"]) & set(GENE_DIRECTED_BUCKETS)]
    if n_act:
        lo, hi = len(gd), len(gd) + n_unc
        obs.append(f"Gene- or target-directed buckets (SOD1, C9orf72, FUS, ATXN2, TDP-43 biology, other gene-directed) "
                   f"cover <b>{len(gd)}</b> of the {n_cls} classified active trials ({pct(len(gd), n_cls)}). Across all "
                   f"{n_act} active drug-type trials the share is only bounded: between {pct(lo, n_act)} and "
                   f"{pct(hi, n_act)}, depending on what the unclassified trials turn out to be.")
    ivc = cov["iv_active"]
    obs.append(f"Modality is explicit (curated asset or keyword) for <b>{ivc['mod_explicit']}</b> of {ivc['n']} active "
               f"test-intervention records; <b>{ivc['mod_fallback']}</b> rely on the registry intervention type "
               f"(e.g. DRUG presumed small molecule) and <b>{ivc['mod_unresolved']}</b> are unresolved.")
    rev = [t for t in a["targets"] if t["flag"]]
    nm = sum(1 for t in a["targets"] if t["mapping_state"] == "mapping not implemented")
    obs.append(f"Target mapping is implemented for {len(a['targets']) - nm} of {len(a['targets'])} Open Targets genes "
               f"listed; {nm} are <b>not assessed</b> and carry no flag. "
               + (("Mapped genes with no trial found via curated assets: "
                   + ", ".join(f"<b>{esc(t['symbol'])}</b>" for t in rev) + ".") if rev else
                  "Every mapped gene has at least one matching trial (any status) in the snapshot.")
               + ((" Mapped genes whose matching trials are all inactive in the snapshot: "
                   + ", ".join(f"<b>{esc(t['symbol'])}</b>" for t in a["targets"]
                               if t["direct_all"] and not t["direct_active"]) + ".")
                  if any(t["direct_all"] and not t["direct_active"] for t in a["targets"]) else ""))
    if a["milestones"]:
        by_year = collections.Counter(t["pcd"]["start"].year for t in a["milestones"])
        obs.append("Registered primary-completion milestones (active, highest phase ≥ 2) by year: "
                   + ", ".join(f"{y}: <b>{c}</b>" for y, c in sorted(by_year.items()))
                   + ". These dates mark the end of primary-outcome data collection as registered, not result announcements.")
    if n_act:
        obs.append(f"<b>{len(a['asia_trials'])}</b> active drug-type trials ({pct(len(a['asia_trials']), n_act)}) list at "
                   f"least one site in Asia; {a['no_sites']} list no sites at all.")
    return obs


# --------------------------------------------------------------------------------------
# HTML rendering
# --------------------------------------------------------------------------------------
def esc(x):
    return html.escape("" if x is None else str(x))


def nct_link(i):
    return f'<a href="https://clinicaltrials.gov/study/{esc(i)}" target="_blank" rel="noopener">{esc(i)}</a>'


def fmt_num(v, missing="not retrieved"):
    return f"<span class='na'>{esc(missing)}</span>" if v is None else esc(v)


SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
       "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]


def seq_color(v, vmax):
    if not v:
        return None
    i = min(len(SEQ) - 1, int(round((v / vmax) ** 0.6 * (len(SEQ) - 1)))) if vmax else 0
    return SEQ[i]


def note(text):
    return f"<p class='note'>{text}</p>"


def caveat(text):
    return f"<div class='caveat'>{text}</div>"


def render_coverage(a):
    cov = a["coverage"]
    ia, ta = cov["iv_active"], cov["trial_active"]
    iall, tall = cov["iv_all"], cov["trial_all"]
    nt = len(a["targets"])
    nm = sum(1 for t in a["targets"] if t["mapping_state"] != "mapping not implemented")

    def tier_rows(counter, tiers):
        return "".join(f"<tr><td class='ind'>{esc(t)}</td><td class='num'>{counter.get(t, 0)}</td></tr>" for t in tiers)

    mech_tiers = [TIER_CURATED_VERIFIED, TIER_CURATED_REGISTRY, TIER_CURATED_UNVERIFIED, TIER_OVERRIDE, TIER_KW_NAME, TIER_NONE]
    mod_tiers = [TIER_CURATED_VERIFIED, TIER_CURATED_REGISTRY, TIER_CURATED_UNVERIFIED, TIER_OVERRIDE, TIER_KW_NAME,
                 TIER_KW_DESC, TIER_TYPE, TIER_NONE]
    roles = cov["roles"]
    role_order = [ROLE_EXPERIMENTAL, ROLE_EXPERIMENTAL_INFERRED, ROLE_UNRESOLVED, ROLE_BACKGROUND, ROLE_COMPARATOR,
                  ROLE_CONTROL, ROLE_PLACEBO, ROLE_NONTHERAPEUTIC]
    return f"""
<div class="two">
<div>
<table class="data compact cov"><thead><tr><th>Trial-level (unique trials)</th><th class='num'>Active</th><th class='num'>All-time</th></tr></thead><tbody>
<tr><td>Interventional trials in snapshot</td><td class='num'>–</td><td class='num'>{cov['n_interventional']}</td></tr>
<tr><td>…with a drug-type test intervention (drug, biologic, genetic, combination product or dietary supplement)</td><td class='num'>{cov['n_active']}</td><td class='num'>{cov['n_drug']}</td></tr>
<tr><td class='ind'>every test intervention has a mechanism</td><td class='num'>{ta['all']}</td><td class='num'>{tall['all']}</td></tr>
<tr><td class='ind'>some test interventions have a mechanism</td><td class='num'>{ta['partial']}</td><td class='num'>{tall['partial']}</td></tr>
<tr class='hl'><td class='ind'><b>no classified mechanism</b></td><td class='num'><b>{ta['none']}</b> ({pct(ta['none'], ta['n'])})</td><td class='num'>{tall['none']} ({pct(tall['none'], tall['n'])})</td></tr>
</tbody></table>
<table class="data compact cov"><thead><tr><th>Units</th><th class='num'>Active</th><th class='num'>All-time</th></tr></thead><tbody>
<tr><td>Test-intervention records (one per registry intervention entry)</td><td class='num'>{ia['n']}</td><td class='num'>{iall['n']}</td></tr>
<tr><td>Unique curated assets or compounds (drug-class rules excluded; uncurated interventions are not de-duplicated)</td><td class='num'>{len(cov['unique_assets_active'])}</td><td class='num'>{len(cov['unique_assets_all'])}</td></tr>
</tbody></table>
<table class="data compact cov"><thead><tr><th>Target subset</th><th class='num'>Genes</th></tr></thead><tbody>
<tr><td>Open Targets genes listed (top {N_OT_TARGETS} by association score)</td><td class='num'>{nt}</td></tr>
<tr><td class='ind'>mapped to trials through curated assets</td><td class='num'>{nm}</td></tr>
<tr><td class='ind'>mapping not implemented (not assessed)</td><td class='num'>{nt - nm}</td></tr>
</tbody></table>
</div>
<div>
<table class="data compact cov"><thead><tr><th>Mechanism, by evidence tier (active test-intervention records)</th><th class='num'>n</th></tr></thead><tbody>
{tier_rows(ia['mech_tiers'], mech_tiers)}
<tr><td class='ind'>…of the unresolved, with a description-keyword candidate for review</td><td class='num'>{ia['mech_candidates']}</td></tr>
</tbody></table>
<table class="data compact cov"><thead><tr><th>Modality, by evidence tier (active test-intervention records)</th><th class='num'>n</th></tr></thead><tbody>
{tier_rows(ia['mod_tiers'], mod_tiers)}
</tbody></table>
<table class="data compact cov"><thead><tr><th>Role of drug-type intervention records (all trials)</th><th class='num'>n</th></tr></thead><tbody>
{''.join(f"<tr><td class='ind'>{esc(r)}{' (included as test)' if r in TEST_ROLES else ' (excluded)'}</td><td class='num'>{roles.get(r, 0)}</td></tr>" for r in role_order)}
</tbody></table>
</div></div>
{note("Mechanism is assigned only from curated assets, manual overrides or keywords in the intervention <i>name</i>. "
      "Keywords found only in a description are stored as review candidates, not labels, because descriptions often "
      "mention other arms or disease models. 'Registry-type fallback' modalities are presumptions (e.g. DRUG → small "
      "molecule), not classifications. GENETIC and COMBINATION_PRODUCT entries with no other evidence stay Unresolved.")}"""


def render_heatmap(a):
    mechs = [m for m in MECHANISMS if m != "Unclassified" and any(a["grid_all"][(m, d)] for d in MODALITIES)]
    mods = [d for d in MODALITIES if any(a["grid_all"][(m, d)] for m in mechs)]
    vmax = max([a["grid_active"][(m, d)] for m in mechs for d in mods] + [1])
    rows = []
    for m in mechs:
        cells = []
        for d in mods:
            v, tot, fb = a["grid_active"][(m, d)], a["grid_all"][(m, d)], a["grid_fallback"][(m, d)]
            col = seq_color(v, vmax)
            ink = "#fff" if col and SEQ.index(col) >= 6 else "var(--text-primary)"
            tip = f"{m} × {d}: {v} active / {tot} all-time unique trials"
            if fb:
                tip += f"; {fb} active trial(s) here rely on registry-type modality fallback"
            ids = a["grid_ids"][(m, d)]
            if ids:
                tip += " | " + ", ".join(ids[:12]) + (" …" if len(ids) > 12 else "")
            style = f"background:{col};color:{ink}" if col else ""
            label = str(v) if v else ("·" if tot else "")
            mark = "<span class='fb'>*</span>" if fb else ""
            sub = f'<span class="sub">{tot}</span>' if tot else ""
            cells.append(f'<td class="hm" style="{style}" title="{esc(tip)}">{label}{mark}{sub}</td>')
        rows.append(f'<tr><th class="rowh">{esc(m)}</th>{"".join(cells)}<td class="tot">{a["mech_active"].get(m, 0)}</td></tr>')
    head = "".join(f"<th class='colh'>{esc(d)}</th>" for d in mods)
    n_act, n_cls = len(a["active"]), len(a["active_classified"])
    return (f'<div class="unc-banner"><span class="big">{n_act - n_cls}</span> of {n_act} active drug-type trials '
            f'({pct(n_act - n_cls, n_act)}) have <b>no classified mechanism</b> and do not appear in this grid. '
            f'The grid covers {n_cls} classified active trials.</div>'
            f'<div class="scrollx"><table class="heat"><thead><tr><th></th>{head}<th class="colh">Unique active trials</th></tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div>'
            + note("Large number = unique active trials (registry status in snapshot); small number = unique trials, all "
                   "statuses; '·' = all-time records only; <b>*</b> = at least one trial in the cell relies on a "
                   "registry-type modality fallback. Each trial counts once per cell, but a trial whose test "
                   f"interventions span several mechanisms or modalities appears in several cells "
                   f"(active trials in more than one cell: {len(a['multi_cell_active'])}; spanning more than one "
                   f"mechanism: {len(a['multi_mech_active'])}), so row totals need not sum to the classified total.")
            + caveat("An empty cell means no matching <i>classified</i> registry record in this snapshot. It does not "
                     "establish commercial white space: unclassified trials, other registries, preclinical and "
                     "undisclosed programmes are all invisible here, and registry activity says nothing about efficacy."))


def render_phase(a):
    mechs = [m for m in MECHANISMS if a["phase_prof"].get(m)]
    head = "".join(f"<th class='num'>{c}</th>" for c in a["phase_cols"])
    rows = []
    for m in mechs:
        pr = a["phase_prof"][m]
        rows.append(f"<tr><th class='rowh'>{esc(m)}</th>"
                    + "".join(f"<td class='num'>{pr.get(c, 0) or ''}</td>" for c in a["phase_cols"]) + "</tr>")
    pu = a["phase_unc"]
    rows.append("<tr class='hl'><th class='rowh'>No classified mechanism</th>"
                + "".join(f"<td class='num'>{pu.get(c, 0) or ''}</td>" for c in a["phase_cols"]) + "</tr>")
    return (f'<table class="data compact"><thead><tr><th>Mechanism</th>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table>'
            + note("Highest listed phase per active trial: Phase 1/2 counts as Ph 2, Phase 2/3 as Ph 3, Early Phase 1 "
                   "as Ph 1. A trial spanning two mechanisms appears in both rows."))


def iv_summary(t, n=3):
    names = list(dict.fromkeys(r["name"] for r in t["test"]))  # repeated cohorts listed once
    return "; ".join(names[:n]) + (" …" if len(names) > n else "")


def date_cell(t):
    if not t["pcd"]:
        return "<span class='na'>none registered</span>"
    prec = " (month only)" if t["pcd"]["precision"] == "month" else ""
    typ = (t["pcd_type"] or "type not given").lower()
    return f"{esc(t['pcd_raw'])}<div class='sub2'>{esc(typ)}{prec}</div>"


def render_trial_table(trials, kind="milestone", empty="No trials in this snapshot."):
    if not trials:
        return f"<p class='note'>{empty}</p>"
    rows = []
    for t in trials:
        extra = ""
        if kind == "milestone":
            extra = (f"<td class='nowrap'>{esc(t['last_update_posted'])}</td>"
                     f"<td>{'yes' if t['has_results'] else 'no'}</td>")
        elif kind == "tdp":
            extra = f"<td>{esc(', '.join(t['direct_targets']) or '—')}</td>"
        rows.append(
            "<tr>"
            + (f"<td class='nowrap'>{date_cell(t)}</td>" if kind in ("milestone", "tdp") else "")
            + f"<td class='nowrap'>{nct_link(t['nct_id'])}</td>"
            f"<td>{esc(iv_summary(t))}"
            + (f"<div class='sub2'>background: {esc('; '.join(t['background']))}</div>" if t["background"] else "")
            + "</td>"
            f"<td>{esc(', '.join(t['mechanisms']))}</td>"
            f"<td>{esc(', '.join(t['modalities']))}</td>"
            f"<td class='nowrap'>{esc(t['phase'])}</td>"
            f"<td>{esc(t['sponsor'])}</td>"
            f"<td class='num'>{esc(t['enrollment'] if t['enrollment'] is not None else '')}</td>"
            f"<td class='nowrap'>{esc(t['status'].replace('_', ' ').title())}</td>"
            + extra + "</tr>")
    head = (("<th>Registered primary completion</th>" if kind in ("milestone", "tdp") else "")
            + "<th>Trial</th><th>Test intervention(s)</th><th>Mechanism</th><th>Modality</th><th>Phase</th>"
              "<th>Lead sponsor</th><th class='num'>n</th><th>Status in snapshot</th>"
            + ("<th>Last update posted</th><th>Results posted on CT.gov</th>" if kind == "milestone" else "")
            + ("<th>Direct target (curated)</th>" if kind == "tdp" else ""))
    return f'<div class="scroll"><table class="data"><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>'


def bar(v, vmax, w=120):
    px = 0 if not vmax or v is None else max(2, round(v / vmax * w))
    return f'<span class="bar" style="width:{px}px"></span>'


def render_targets(a, ot):
    if not a["targets"]:
        return "<p class='note'><b>Open Targets data unavailable</b> (not retrieved). No target statements are made.</p>"
    tmax = max((t["ot_score"] or 0) for t in a["targets"]) or 1
    rows = []
    for t in a["targets"]:
        gs = t["genetic_score"]
        nm = t["mapping_state"] == "mapping not implemented"
        if nm:
            trial_cell = "<span class='na'>not assessed</span>"
        else:
            ids = t["direct_ids_active"]
            trial_cell = (f"<b>{t['direct_active']}</b> active / {t['direct_all']} all-time"
                          + (f"<div class='sub2'>{', '.join(nct_link(i) for i in ids[:6])}{' …' if len(ids) > 6 else ''}</div>" if ids else ""))
        proxy = (f"{esc(t['bucket'])}: {t['bucket_active']} active / {t['bucket_all']} all-time"
                 f"<div class='sub2'>bucket-level proxy, not gene-specific</div>") if t["bucket"] else ""
        pm_t = fmt_num(t["pubmed_total"], t["pubmed_status"] if t["pubmed_total"] is None else "")
        pm_r = fmt_num(t["pubmed_recent"], t["pubmed_status"] if t["pubmed_recent"] is None else "")
        state_cls = "st-na" if nm else ("st-warn" if "no matching" in t["mapping_state"] else "st-pass")
        rows.append(
            f"<tr><td><b>{esc(t['symbol'])}</b><div class='sub2'>{esc(t['name'][:48])}</div></td>"
            f"<td class='nowrap'>{bar(t['ot_score'], tmax)} {(t['ot_score'] or 0):.2f}</td>"
            f"<td class='num'>{fmt_num(None if gs is None else f'{gs:.2f}', 'not reported')}</td>"
            f"<td class='num'>{pm_t}</td><td class='num'>{pm_r}</td>"
            f"<td><span class='st {state_cls}'>{esc(t['mapping_state'])}</span>"
            + (f"<div class='sub2'>assets: {esc(t['curated_assets'])}</div>" if t["curated_assets"] else "") + "</td>"
            f"<td>{trial_cell}</td><td>{proxy}</td>"
            f"<td>{'<span class=tag>' + esc(t['flag']) + '</span>' if t['flag'] else ''}</td></tr>")
    off = a["assessed_off_list"]
    return ('<div class="scroll"><table class="data"><thead><tr><th>Target</th><th>Open Targets overall</th>'
            '<th class="num">Genetic association</th><th class="num">PubMed co-mentions (all)</th>'
            '<th class="num">PubMed (recent window)</th><th>Mapping state</th><th>Trials via curated direct assets</th>'
            '<th>Bucket-level proxy</th><th></th></tr></thead><tbody>' + "".join(rows) + "</tbody></table></div>"
            + note("Mapping is implemented only for genes that a curated asset targets directly "
                   f"({', '.join(a['assessed_genes'])}). Counts are unique trials with an experimental-role intervention "
                   "whose curated direct target is the gene; unrecognised assets aimed at the same gene would be missed. "
                   "The TDP-43 bucket proxy includes downstream approaches (e.g. STMN2 restoration) and is shown separately "
                   "from direct TARDBP relationships."
                   + (f" Curated direct targets outside this Open Targets list: {', '.join(off)}." if off else ""))
            + caveat("Open Targets association scores are prioritisation evidence, not probabilities of therapeutic "
                     "success. PubMed co-mention counts measure literature volume, not evidence quality or direction. "
                     "'Not assessed' and 'not retrieved' mean no statement is made; they are not zeros. A mapped gene "
                     "with no matching trial does not establish absence of development: preclinical, non-registered "
                     "and other-registry programmes are not visible."))


def render_bars(items):
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


def render_validation(vs, pack_rows):
    n_pack = len(pack_rows or [])
    n_prop = sum(1 for r in (pack_rows or []) if r.get("proposed_mechanism") or r.get("proposed_modality"))
    n_human = sum(1 for r in (pack_rows or []) if (r.get("human_review_status") or "").strip())
    body = f"<div class='acc-box acc-{vs['state'].replace(' ', '-')}'><b>{esc(vs['headline'])}</b>"
    if vs["detail"]:
        body += f"<div>{esc(vs['detail'])}</div>"
    if vs["state"] == "measured":
        for subset, label in (("heldout", "Held-out subset (independent of rule tuning)"),
                              ("dev", "Development subset (NOT independent: may have informed rules)")):
            r = vs["results"]["results"][subset]
            parts = []
            for axis in ("mechanism", "modality"):
                x = r[axis]
                parts.append(f"{axis}: {x['agree']}/{x['n']}" + (f" ({x['rate']:.0%})" if x["rate"] is not None else ""))
            body += f"<div>{label}: {'; '.join(parts)}</div>"
    body += "</div>"
    body += note(f"Review pack: {n_pack} sampled trials; AI-proposed reference labels drafted for {n_prop} "
                 f"(not human-approved); human review fields completed for {n_human}. Agreement is computed only from "
                 "human-approved rows, separately for mechanism and modality, each with its own denominator.")
    return body


CSS = """
.viz-root{color-scheme:light;--surface-0:#f6f5f2;--surface-1:#fcfcfb;--text-primary:#0b0b0b;--text-secondary:#52514e;
--muted:#6b6a65;--rule:#e4e2dc;--accent:#2a78d6;--good:#0b7d0b;--warn:#8a5d00;--crit:#c0392b;--hl:#fff6e0}
@media (prefers-color-scheme:dark){.viz-root{color-scheme:dark;--surface-0:#121211;--surface-1:#1a1a19;
--text-primary:#f0efec;--text-secondary:#c3c2b7;--muted:#a3a299;--rule:#33332f;--accent:#5aa0ef;--hl:#2c2617;--warn:#e0b04a}}
*{box-sizing:border-box}body{margin:0}
.viz-root{background:var(--surface-0);color:var(--text-primary);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",
Inter,Helvetica,Arial,sans-serif;padding:32px 20px 64px}
.wrap{max-width:1180px;margin:0 auto}
h1{font-size:26px;margin:0 0 4px;letter-spacing:-.01em}h2{font-size:18px;margin:0 0 4px}
.lede{color:var(--text-secondary);margin:0 0 14px;max-width:860px}
.meta{color:var(--text-secondary);font-size:12px;margin-bottom:18px;border:1px solid var(--rule);border-radius:8px;padding:10px 14px;background:var(--surface-1)}
.meta table td{padding:1px 10px 1px 0;vertical-align:top}.meta table td:first-child{white-space:nowrap;color:var(--muted)}
.card{background:var(--surface-1);border:1px solid var(--rule);border-radius:10px;padding:20px 22px;margin:0 0 18px}
.card .why{color:var(--text-secondary);margin:0 0 14px;font-size:13px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin-bottom:18px}
.kpi{background:var(--surface-1);border:1px solid var(--rule);border-radius:10px;padding:14px 16px}
.kpi .v{font-size:28px;font-weight:650;line-height:1.1}.kpi .l{color:var(--text-secondary);font-size:12px}
.kpi.warnk{border-color:#d9a400}
ul.obs{margin:0;padding-left:18px}ul.obs li{margin:6px 0}
table{border-collapse:collapse;width:100%}
table.data th,table.data td{border-bottom:1px solid var(--rule);padding:6px 8px;text-align:left;vertical-align:top;font-size:13px}
table.data thead th{color:var(--text-secondary);font-weight:600;font-size:12px;position:sticky;top:0;background:var(--surface-1)}
table.compact{width:auto}table.compact td,table.compact th{padding:4px 12px}
table.cov{margin-bottom:12px;width:100%}td.ind{padding-left:22px!important}tr.hl td,tr.hl th{background:var(--hl)}
.num{text-align:right;font-variant-numeric:tabular-nums}.nowrap{white-space:nowrap}
.scroll{max-height:520px;overflow:auto;border:1px solid var(--rule);border-radius:8px}.scrollx{overflow-x:auto}
table.heat{width:auto;border-collapse:separate;border-spacing:2px}
table.heat td.hm{width:84px;height:40px;text-align:center;font-weight:650;font-size:15px;border-radius:4px;
background:transparent;outline:1px solid var(--rule);outline-offset:-1px;cursor:default;font-variant-numeric:tabular-nums}
table.heat td.hm .sub{display:block;font-size:10px;font-weight:400;opacity:.75}
table.heat td.hm .fb{font-size:11px;vertical-align:super}
table.heat th.colh{font-size:11px;color:var(--text-secondary);font-weight:600;vertical-align:bottom;padding:0 2px 6px;width:84px;max-width:84px;overflow-wrap:anywhere;hyphens:auto}
table.heat th.rowh{text-align:right;font-size:12px;padding-right:10px;font-weight:500;white-space:nowrap}
table.heat td.tot{padding-left:10px;font-weight:650;color:var(--text-secondary);font-variant-numeric:tabular-nums}
.rowh{font-weight:500}
.note{color:var(--muted);font-size:12px;margin:10px 0 0}
.caveat{font-size:12px;margin:10px 0 0;padding:8px 12px;border-left:3px solid #d9a400;background:var(--hl);color:var(--text-primary)}
.sub2{color:var(--muted);font-size:11px}.na{color:var(--muted);font-style:italic}
.bar{display:inline-block;height:10px;background:var(--accent);border-radius:0 4px 4px 0;vertical-align:middle}
table.bars{width:auto}table.bars th{text-align:right;font-weight:500;padding:3px 10px 3px 0;font-size:13px}
table.bars td{padding:3px 6px}
.two{display:grid;grid-template-columns:1fr 1fr;gap:18px}@media(max-width:820px){.two{grid-template-columns:1fr}}
.tag{font-size:11px;border:1px solid var(--warn);color:var(--warn);border-radius:10px;padding:0 7px;white-space:nowrap}
.st{font-size:11px;font-weight:650;white-space:nowrap}.st-pass{color:var(--good)}.st-fail{color:var(--crit)}
.st-warn{color:var(--warn)}.st-todo,.st-na{color:var(--muted)}
.banner{background:#fff3cd;color:#4a3700;border:1px solid #e9c46a;border-radius:8px;padding:10px 14px;margin-bottom:16px}
.unc-banner{background:var(--hl);border:1px solid #d9a400;border-radius:8px;padding:10px 14px;margin:0 0 14px}
.unc-banner .big{font-size:22px;font-weight:700;margin-right:6px}
.acc-box{border:1px solid var(--rule);border-radius:8px;padding:10px 14px;margin:6px 0}
.acc-not-measured,.acc-stale{border-color:#d9a400;background:var(--hl)}
h3{font-size:14px;margin:16px 0 6px}
a{color:var(--accent)}
footer{color:var(--muted);font-size:12px;margin-top:28px}
"""


def render_html(ct, ot, pm, trials, a, checks, regress, obs, vs, pack_rows, prov, fixture):
    cov = a["coverage"]
    ms = a["ms_summary"]
    n_act, n_cls = cov["n_active"], len(a["active_classified"])
    kpis = [
        (cov["n_interventional"], "interventional ALS trials in snapshot (unique trials, all statuses)", ""),
        (cov["n_drug"], "with a drug-type test intervention (drug, biologic, genetic, combination product or dietary supplement)", ""),
        (n_act, "of those, active by registry status in snapshot", ""),
        (f"{n_cls}/{n_act}", f"active trials with a classified mechanism ({pct(n_cls, n_act)})", "warnk"),
        (ms["in_window"], f"registered primary-completion milestones, {a['analysis_date']:%b %Y} – {a['horizon']:%b %Y}", ""),
    ]
    kpi_html = "".join(f"<div class='kpi {c}'><div class='v'>{esc(v)}</div><div class='l'>{esc(l)}</div></div>" for v, l, c in kpis)
    banner = ("<div class='banner'><b>SYNTHETIC TEST FIXTURE - NOT PUBLIC DATA.</b> Built from invented records to "
              "exercise the pipeline. No statement on this page describes the real ALS landscape.</div>") if fixture else ""
    banner += ("<div class='banner'><b>Exploratory prototype.</b> Classification accuracy not yet measured against "
               "human-reviewed labels; rule-based labels are provisional. Not investment advice.</div>")
    meta_rows = "".join(f"<tr><td>{esc(k)}</td><td>{v}</td></tr>" for k, v in prov)
    ms_text = (f"Of {ms['pool']} active trials with highest phase ≥ 2 (Phase 1/2 is included; it counts as Phase 2), "
               f"<b>{ms['in_window']}</b> have a registered primary-completion date inside the window "
               f"({ms['in_window_month_precision']} recorded to month precision only; {ms['in_window_current_month']} "
               f"fall in a period that contains the analysis date). <b>{ms['passed']}</b> have a registered date whose "
               f"whole period ended before the analysis date ({ms['passed_actual']} marked actual, "
               f"{ms['passed_estimated']} marked estimated); {ms['beyond']} fall after the window; {ms['no_date']} have none.")
    passed_list = render_trial_table(a["passed"], kind="milestone", empty="None.")
    sections = f"""
<div class="card"><h2>1 · Coverage and counting units</h2>
<p class="why">Read this first. It shows how much of the trial set the rules could classify, on what evidence, and
what each number counts.</p>
{render_coverage(a)}</div>

<div class="card"><h2>Key observations</h2>
<p class="why">Generated by fixed rules from the tables below (no language model). Each carries its denominator.</p>
<ul class="obs">{''.join(f'<li>{o}</li>' for o in obs)}</ul></div>

<div class="card"><h2>2 · Classified trial activity by mechanism and modality</h2>
<p class="why">Unique active trials per mechanism × modality, restricted to trials with at least one classified test
intervention.</p>
{render_heatmap(a)}</div>

<div class="two">
<div class="card"><h2>3 · Phase profile of active trials</h2>
<p class="why">Highest listed phase per active trial, classified subset, with unclassified trials shown separately.</p>
{render_phase(a)}</div>
<div class="card"><h2>4 · Sponsors and sites (active drug-type trials)</h2>
<p class="why">Lead-sponsor class, top industry lead sponsors, and countries listed among sites. {len(a['asia_trials'])} of
{n_act} list a site in Asia; {a['no_sites']} list no sites.</p>
{render_bars(most_common_sorted(a['sponsors'], 6))}<div style="height:12px"></div>{render_bars(most_common_sorted(a['top_sponsors'], 8))}
<div style="height:12px"></div>{render_bars(most_common_sorted(a['countries'], 10))}
{caveat("Site lists come from ClinicalTrials.gov and may be incomplete or out of date; trials registered only in "
        "jRCT, ChiCTR, CTIS or other registries are absent, which under-represents Asian activity in particular.")}</div>
</div>

<div class="card"><h2>5 · Registered primary-completion milestones, {a['analysis_date']:%d %b %Y} – {a['horizon']:%d %b %Y}</h2>
<p class="why">{ms_text}</p>
{caveat("A primary-completion date is when collection of primary-outcome data is registered to finish. It does not say "
        "when, or whether, results will be announced or published. Estimated dates often move. A month-only date is "
        "treated as covering the whole month, so it is never counted as passed before that month has ended. "
        "'Status in snapshot' is the status recorded on the retrieval date; registry update history was not examined, "
        "so a passed date alone does not establish missing results or a stale record.")}
{render_trial_table(a['milestones'], kind='milestone')}
<h3>Active Phase 2+ trials whose registered primary-completion period has passed (as of the analysis date)</h3>
{passed_list}</div>

<div class="card"><h2>6 · Target evidence and mapped clinical activity</h2>
<p class="why">Top {len(a['targets'])} ALS-associated targets from Open Targets ({esc(ot.get('disease', {}).get('id', 'not retrieved'))}),
literature volume from PubMed, and trial activity only where a curated direct relationship exists.</p>
{render_targets(a, ot)}</div>

<div class="card"><h2>7 · TDP-43 biology: classified trials</h2>
<p class="why">Trials whose test intervention is classified in the TDP-43 biology bucket. The 'direct target' column
separates approaches aimed at TDP-43 itself (TARDBP) from downstream ones (e.g. STMN2 restoration).</p>
{render_trial_table(sorted(a['tdp'], key=lambda t: (t['status_group'] != 'Active', -(t['phase_rank'] or 0), t['nct_id'])), kind='tdp')}
{caveat("Only trials the rules could classify appear here. Programmes acting on TDP-43 indirectly through other "
        "mechanisms, unrecognised asset names and preclinical work are not captured.")}</div>

<div class="card"><h2>8 · Method, validation and reproducibility</h2>
<p class="why">Three different things, reported separately: software/data checks, classification coverage (section 1),
and agreement with human-reviewed labels.</p>
<h3>a · Agreement with human-reviewed labels</h3>
{render_validation(vs, pack_rows)}
<h3>b · Software and data-integrity checks</h3>
{render_checks(checks)}
<h3>c · Known-answer regression checks</h3>
{note("These confirm that specific, previously fixed cases still behave as expected. They are not evidence of overall accuracy.")}
{render_checks(regress)}
<h3>d · Method</h3>
<p class="note">Interventions are first given a role from the registry's arm structure (experimental, background,
comparator, control-arm, placebo/sham, non-therapeutic, or unresolved). Only test-role drug-type interventions
(experimental, experimental-inferred, unresolved-role) drive trial labels; background, comparator and control-arm
treatments are kept in separate columns. Each test intervention is then labelled by: manual override → curated asset table
(data/curation/asset_curation.csv: {len(load_curation())} rules, each with its source status) → name keywords ({len(KEYWORD_RULES)}
rules) → description keywords (modality only; mechanism as candidate) → registry-type fallback (modality only).
The tier and rule behind every label are in data/interventions.csv.</p>
{caveat("ClinicalTrials.gov coverage is not a complete global or preclinical landscape. Registry activity does not "
        "establish clinical efficacy. This prototype has not been scientifically validated.")}</div>
"""
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>ALS Landscape Navigator{' - SYNTHETIC TEST FIXTURE' if fixture else ''}</title><style>{CSS}</style></head>
<body><div class="viz-root"><div class="wrap">
{banner}
<h1>ALS Landscape Navigator{' <span style="color:#c0392b">(synthetic test fixture)</span>' if fixture else ''}</h1>
<p class="lede">An exploratory view of registered ALS interventional trials, joined to Open Targets genetic evidence and
PubMed literature volume. Every label is traceable to a rule or curated source; coverage and limitations are shown
next to each result.</p>
<div class="meta"><table>{meta_rows}</table></div>
<div class="kpis">{kpi_html}</div>
{sections}
<footer>Sources: ClinicalTrials.gov API v2 · Open Targets Platform GraphQL API · NCBI PubMed E-utilities.
Built by als_navigator.py (standard-library Python) from cached raw responses in data/raw/.</footer>
</div></div></body></html>"""


# --------------------------------------------------------------------------------------
# CSV outputs
# --------------------------------------------------------------------------------------
def write_csvs(trials, a, data_dir):
    def w_csv(name, header, rows):
        with open(os.path.join(data_dir, name), "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(header)
            w.writerows(rows)

    w_csv("trials.csv",
          ["nct_id", "title", "status_in_snapshot", "status_group", "phase", "phase_rank", "enrollment", "lead_sponsor",
           "sponsor_class", "start", "primary_completion", "primary_completion_type", "primary_completion_precision",
           "last_update_posted", "results_posted", "countries", "has_drug_type_test_intervention",
           "test_interventions", "mechanisms", "mechanism_status", "modalities", "curated_assets",
           "curated_direct_targets", "background_therapy", "comparator_or_control_treatments", "rules_fired"],
          [[t["nct_id"], t["title"], t["status"], t["status_group"], t["phase"], t["phase_rank"],
            "" if t["enrollment"] is None else t["enrollment"], t["sponsor"], t["sponsor_class"], t["start_raw"],
            t["pcd_raw"], t["pcd_type"], t["pcd"]["precision"] if t["pcd"] else "", t["last_update_posted"],
            t["has_results"], "; ".join(t["countries"]), t["is_drug"],
            "; ".join(r["name"] for r in t["test"]), "; ".join(t["mechanisms"]), t["mech_status"],
            "; ".join(t["modalities"]), "; ".join(t["assets"]), "; ".join(t["direct_targets"]),
            "; ".join(t["background"]), "; ".join(t["comparators"]),
            "; ".join(f"{r['mechanism_rule'] or '-'}/{r['modality_rule']}" for r in t["test"])] for t in trials])

    w_csv("interventions.csv",
          ["nct_id", "iv_index", "intervention_name", "registry_type", "arm_labels", "arm_types", "role", "role_basis",
           "classified_text", "background_components", "mechanism", "mechanism_tier", "mechanism_rule",
           "mechanism_candidate", "candidate_rule", "modality", "modality_tier", "modality_rule", "asset_id",
           "direct_targets", "target_relationship", "source_urls", "flags"],
          [[r["nct_id"], r["iv_index"], r["name"], r["registry_type"], "; ".join(r["arm_labels"]),
            "; ".join(r["arm_types"]), r["role"], r["role_basis"], r.get("classified_text", ""),
            "; ".join(r.get("background_components", [])), r.get("mechanism", ""), r.get("mechanism_tier", ""),
            r.get("mechanism_rule", ""), r.get("mechanism_candidate") or "", r.get("candidate_rule") or "",
            r.get("modality", ""), r.get("modality_tier", ""), r.get("modality_rule", ""), r.get("asset_id") or "",
            "; ".join(r.get("direct_targets", [])), r.get("target_relationship", ""), r.get("source_urls", ""),
            " | ".join(r["flags"])]
           for t in trials for r in t["records"] if r["drug_type"]])

    w_csv("targets.csv",
          ["symbol", "name", "ensembl", "ot_score", "genetic_association_score", "literature_score",
           "pubmed_total", "pubmed_recent", "pubmed_status", "mapping_state", "curated_assets",
           "direct_active_trials", "direct_all_trials", "direct_active_trial_ids",
           "bucket_proxy", "bucket_proxy_active_trials", "bucket_proxy_all_trials", "flag"],
          [[t["symbol"], t["name"], t["ensembl"], t["ot_score"], "" if t["genetic_score"] is None else t["genetic_score"],
            "" if t["literature_score"] is None else t["literature_score"],
            "" if t["pubmed_total"] is None else t["pubmed_total"], "" if t["pubmed_recent"] is None else t["pubmed_recent"],
            t["pubmed_status"], t["mapping_state"], t["curated_assets"],
            "not assessed" if t["direct_active"] is None else t["direct_active"],
            "not assessed" if t["direct_all"] is None else t["direct_all"], "; ".join(t["direct_ids_active"]),
            t["bucket"], "" if t["bucket_active"] is None else t["bucket_active"],
            "" if t["bucket_all"] is None else t["bucket_all"], t["flag"]] for t in a["targets"]])

    w_csv("milestones.csv",
          ["nct_id", "position_vs_analysis_date", "primary_completion", "primary_completion_type", "precision",
           "status_in_snapshot", "last_update_posted", "results_posted", "phase", "test_interventions", "mechanisms",
           "modalities", "lead_sponsor"],
          [[t["nct_id"], t["_ms_pos"], t["pcd_raw"], t["pcd_type"], t["pcd"]["precision"] if t["pcd"] else "",
            t["status"], t["last_update_posted"], t["has_results"], t["phase"], "; ".join(r["name"] for r in t["test"]),
            "; ".join(t["mechanisms"]), "; ".join(t["modalities"]), t["sponsor"]]
           for t in sorted(a["milestones"] + a["passed"], key=lambda t: (t["_ms_pos"], t["pcd"]["start"], t["nct_id"]))])

    w_csv("unclassified_interventions.csv",
          ["intervention_name", "n_trials", "registry_type", "current_modality", "modality_tier",
           "description_candidate_mechanism", "candidate_rule", "mechanism", "modality", "source", "rationale"],
          [[name, n, *a["review_meta"][name], "", "", "", ""] for name, n in most_common_sorted(a["review_queue"])])


# --------------------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------------------
def resolve_analysis_date(ct, explicit):
    if explicit:
        return dt.date.fromisoformat(explicit), "set by --analysis-date"
    fa = (ct.get("fetched_at") or "")[:10]
    if fa:
        return dt.date.fromisoformat(fa), "defaults to the ClinicalTrials.gov retrieval date"
    raise SystemExit("No analysis date: snapshot has no fetched_at; pass --analysis-date.")


def build_timestamp():
    sde = os.environ.get("SOURCE_DATE_EPOCH")
    if sde:
        return dt.datetime.fromtimestamp(int(sde), dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).strftime("%Y-%m-%dT%H:%M:%SZ")


def build(ct, ot, pm, data_dir, out_dir, analysis_date, date_source, fixture=False, raw_paths=None,
          html_name="als_landscape.html", write_pack=True):
    curation = load_curation()
    overrides = load_overrides(data_dir)
    trials = normalise_trials(ct, curation, overrides)
    a = analyse(ot, pm, trials, curation, analysis_date)
    checks = integrity_checks(ct, trials, a)
    regress = regression_answer_checks(trials) if not fixture else []
    obs = observations(a)
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    write_csvs(trials, a, data_dir)
    pack_rows = build_review_pack(data_dir, ct, trials) if write_pack else None
    vs = validation_status(os.path.join(data_dir, "validation"), trials)
    built = build_timestamp()
    hashes = {os.path.basename(p): sha256_file(p) for p in (raw_paths or []) if os.path.exists(p)}
    hashes["asset_curation.csv"] = sha256_file(CURATION_PATH)
    pm_window = (pm or {}).get("window")
    prov = [
        ("Retrieved", f"ClinicalTrials.gov {esc((ct.get('fetched_at') or 'unknown')[:10])} · Open Targets "
                      f"{esc((ot.get('fetched_at') or 'not retrieved')[:10])} · PubMed {esc(((pm or {}).get('fetched_at') or 'not retrieved')[:10])}"),
        ("Queries and limits", "<br>".join(f"{esc(k)}: {esc(v)}" for k, v in QUERY_DESCRIPTION.items())
         + (f"<br>PubMed recent window in this snapshot: {pm_window[0]}–{pm_window[1]}" if pm_window else "")),
        ("Analysis date", f"{analysis_date.isoformat()} ({esc(date_source)}); all date-sensitive results use this date"),
        ("Classification", f"version {CLASSIFICATION_VERSION}; trial-label fingerprint {fingerprint(trials)}"),
        ("Validation status", esc(vs["headline"])),
        ("Build", f"{built} · input SHA-256: " + ", ".join(f"{esc(k)} {v[:12]}" for k, v in sorted(hashes.items()))),
    ]
    if fixture:
        prov[0] = ("Retrieved", "<b>synthetic fixture: no public data retrieved</b>")
        prov[1] = ("Queries and limits", "not applicable (invented records in tests/fixtures.py)")
    page = render_html(ct, ot, pm, trials, a, checks, regress, obs, vs, pack_rows, prov, fixture)
    out_html = os.path.join(out_dir, html_name)
    with open(out_html, "w") as f:
        f.write(page)
    analytic = {}
    for name in ("trials.csv", "interventions.csv", "targets.csv", "milestones.csv", "unclassified_interventions.csv"):
        analytic[name] = sha256_file(os.path.join(data_dir, name))
    manifest = {
        "classification_version": CLASSIFICATION_VERSION, "analysis_date": analysis_date.isoformat(),
        "analysis_date_source": date_source, "build_timestamp_utc": built, "fixture": fixture,
        "source_retrieved": {"ClinicalTrials.gov": ct.get("fetched_at"), "Open Targets": ot.get("fetched_at"),
                             "PubMed": (pm or {}).get("fetched_at")},
        "queries": QUERY_DESCRIPTION, "input_sha256": hashes, "output_sha256": analytic,
        "trials_fingerprint": fingerprint(trials), "validation_state": vs["state"],
        "counts": {"interventional": a["coverage"]["n_interventional"], "drug_type": a["coverage"]["n_drug"],
                   "active": a["coverage"]["n_active"], "active_classified": len(a["active_classified"]),
                   "milestones_in_window": a["ms_summary"]["in_window"]},
    }
    with open(os.path.join(out_dir, "run_manifest.json" if not fixture else "run_manifest_FIXTURE.json"), "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    return {"trials": trials, "analysis": a, "checks": checks, "regress": regress, "html": out_html,
            "manifest": manifest, "validation": vs}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true", help="ignore cache and re-download (creates a new snapshot)")
    ap.add_argument("--analysis-date", help="YYYY-MM-DD; default = ClinicalTrials.gov retrieval date of the snapshot")
    ap.add_argument("--evaluate", action="store_true", help="score HUMAN-approved labels in the review pack")
    args = ap.parse_args(argv)

    ct_path = os.path.join(RAW_DIR, "ctgov_studies.json")
    ot_path = os.path.join(RAW_DIR, "opentargets.json")
    pm_path = os.path.join(RAW_DIR, "pubmed_counts.json")
    try:
        ct = load_or_fetch(ct_path, fetch_ctgov, args.refresh)
    except Exception as e:
        sys.exit(f"Could not load or download ClinicalTrials.gov data: {e}")
    analysis_date, date_source = resolve_analysis_date(ct, args.analysis_date)
    try:
        ot = load_or_fetch(ot_path, fetch_opentargets, args.refresh)
    except Exception as e:
        print(f"[Open Targets] unavailable: {e}")
        ot = {"fetched_at": None, "rows": []}
    symbols = [r["target"]["approvedSymbol"] for r in ot.get("rows", [])][:N_PUBMED_TARGETS]
    symbols += [g for g in REQUIRED_PUBMED_GENES if g not in symbols]
    try:
        pm = load_or_fetch(pm_path, lambda: fetch_pubmed(symbols, analysis_date.year), args.refresh)
    except Exception as e:
        print(f"[PubMed] unavailable: {e}")
        pm = {"fetched_at": None, "counts": None}

    if args.evaluate:
        trials = normalise_trials(ct, load_curation(), load_overrides(DATA_DIR))
        cmd_evaluate(DATA_DIR, trials)
        return

    res = build(ct, ot, pm, DATA_DIR, OUT_DIR, analysis_date, date_source,
                fixture=bool(ct.get("_fixture")), raw_paths=[ct_path, ot_path, pm_path])
    a = res["analysis"]
    print(f"Classification v{CLASSIFICATION_VERSION} | analysis date {analysis_date} ({date_source})")
    print(f"  interventional {a['coverage']['n_interventional']} | drug-type {a['coverage']['n_drug']} | "
          f"active {a['coverage']['n_active']} | active with classified mechanism {len(a['active_classified'])}")
    for c in res["checks"] + res["regress"]:
        print(f"  [{c['status']}] {c['name']}: {c['detail']}")
    print(f"  Validation: {res['validation']['headline']}")
    print(f"\n  Report:   {res['html']}\n  Manifest: {os.path.join(OUT_DIR, 'run_manifest.json')}")


if __name__ == "__main__":
    main()
