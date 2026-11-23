# ALS Landscape Navigator (v2, repaired): exploratory prototype

**[▶ Open the live report](https://sourvegie.github.io/als-landscape-navigator/)** · [What was fixed in v2](CHANGELOG.md) · [Test results](docs/test_results.txt)

**In brief.** A Python pipeline (standard library only) that pulls 724 interventional ALS trials from ClinicalTrials.gov, the top 30 ALS-associated targets from Open Targets and literature counts from PubMed. It classifies each trial's test interventions by mechanism and modality, and builds one self-contained interactive HTML report. Every label records the rule and evidence tier behind it. The build is reproducible from a cached snapshot and covered by 31 regression tests plus a mutation check that confirms the tests catch the v1 defects.

An exploratory scientific-intelligence prototype. It joins registered ALS interventional trials (ClinicalTrials.gov), ALS-associated targets (Open Targets) and literature volume (PubMed) into one self-contained HTML report plus CSV tables. Every label can be traced to a rule or a curated source.

> **Status:** exploratory. **Classification accuracy has not yet been measured** against human-reviewed labels. Programme notes are provisional and pending scientific review. This is not investment advice.

| Source | What it contributes | Snapshot retrieved |
|---|---|---|
| ClinicalTrials.gov API v2 | ALS studies: phase, status, arms, interventions, dates, sites | 2026-11-18 |
| Open Targets Platform (GraphQL) | Top 30 ALS-associated targets and association scores (prioritisation evidence) | 2026-11-18 |
| PubMed E-utilities | Title/abstract co-mention counts per target (literature volume only) | 2026-11-18 |

## Reproduce the report from the cached snapshot

Requires Python 3.8+ and the standard library only.

```bash
python3 als_navigator.py                          # analysis date defaults to the CT.gov retrieval date (2026-11-18)
python3 -m unittest discover -s tests -v          # 31 regression checks
python3 tests/mutation_smoke.py                   # confirms the tests catch re-introduced v1 defects
python3 tools/build_test_preview.py               # SYNTHETIC preview -> test_preview/ (kept apart from output/)
```

The analysis date controls every date-sensitive result. Override it explicitly with `--analysis-date YYYY-MM-DD`. Two runs with the same inputs and analysis date produce byte-identical CSVs. Set `SOURCE_DATE_EPOCH` as well if you want a byte-identical HTML build stamp too.

`--refresh` downloads a **new** snapshot. That changes the inputs, so results will differ and any validation results will be marked stale. Back up `data/raw/` first if you want to keep the 2026-11-18 snapshot.

## Outputs

| File | Contents |
|---|---|
| `output/als_landscape.html` | Self-contained report: coverage, observations, classified activity grid, phase profile, sponsors and sites, registered primary-completion milestones, target evidence with mapping states, TDP-43 trials, method and validation |
| `output/run_manifest.json` | Retrieval dates, queries and limits, analysis date, classification version, validation state, build time, input and output SHA-256 hashes |
| `data/trials.csv` | One row per interventional trial: test interventions, labels, background and comparator treatments, curated assets and targets, date type and precision |
| `data/interventions.csv` | One row per drug-type intervention record: role and its basis, label, **evidence tier**, rule, description-keyword candidate, curated source URLs, flags |
| `data/targets.csv` | Open Targets genes with mapping state, direct curated-asset trial counts ("not assessed" when unmapped) and a labelled bucket proxy |
| `data/milestones.csv` | Active Phase 2+ trials with registered primary-completion dates in the window or already passed |
| `data/unclassified_interventions.csv` | Review queue: interventions with no mechanism or an unresolved modality, with description candidates |
| `data/curation/asset_curation.csv` | Curated asset rules, each with mechanism, modality, direct targets, relationship, source status, source URLs and rationale |
| `data/validation/` | Frozen evaluation split, AI-proposed labels, review pack, and (after human review) `validation_results.json` |
| `docs/programme_assessments.md` | Provisional, source-linked notes on QRL-201, AP-101 and CNM-Au8 |
| `docs/REVIEW_TASKS.md` | Outstanding scientific questions and human-review tasks |
| `docs/test_results.txt` | Output of the regression checks and mutation smoke check for this build |
| `CHANGELOG.md` | Analytical errors corrected since v1 |
| `data/original_outputs/` | v1 code and outputs, unchanged, for comparison |
| `test_preview/` | **Synthetic** fixture preview. Not public data. |

## How classification works

```
ClinicalTrials.gov study (interventional only)
  → role per intervention from arm types: experimental | experimental (inferred) | unresolved
                                          background | comparator | control-arm | placebo/sham | non-therapeutic
  → test-role drug-type interventions only (DRUG, BIOLOGICAL, GENETIC, COMBINATION_PRODUCT, DIETARY_SUPPLEMENT)
      split "X + riluzole" style names: the standard-of-care part becomes background
  → label, by first evidence found:
      manual override (data/overrides.csv, with source)
      curated asset table (verified source | registry text only | v1 rule not yet re-verified)
      keyword in the intervention NAME
      keyword in the description: modality only; a mechanism found here becomes a review CANDIDATE
      registry-type fallback: modality only (DRUG → small molecule presumed); GENETIC/COMBINATION stay Unresolved
  → trial labels = union over test interventions; trials counted once per category
```

**Counting units.** *Unique trials* are used for all trial counts, grid cells and phase rows. *Intervention records* are registry intervention entries, so repeated cohorts count separately here. *Unique curated assets* are de-duplicated only where a curated asset rule matched.

**Target mapping.** A gene is "mapped" only if a curated asset targets it directly. Otherwise it is "not assessed" and gets no flag. The TDP-43 bucket count is a labelled proxy and is never attributed to TARDBP.

## Validation workflow (human review required)

1. Open `data/validation/label_review_pack.csv`. Each of the 50 sampled trials has a registry link, arm context, roles, current predictions with rules and tiers, and an **AI-proposed** label with reasoning and sources (`proposed_*`). The proposals are drafts, not a reference standard. `triage_pred_vs_proposal` is a triage aid only.
2. A scientist fills `human_mechanism`, `human_modality` (use `N/A` when no therapeutic intervention exists), `human_review_status` (`approved` / `corrected` / `unresolved`), `human_reviewer` and `human_review_date`. Re-running the pipeline keeps these fields.
3. `python3 als_navigator.py --evaluate` scores only approved or corrected rows, **separately per axis with separate denominators**, and reports **held-out** and **dev** subsets separately (`data/validation/eval_split.csv`, frozen before any tuning).
4. Rules or overrides may be changed using **dev** rows only. Never tune on held-out rows. Agreement on rows used for tuning is not independent validation.
5. Results are tied to the classifier version and a label fingerprint. Any change marks them **stale** on the report.

The three report panels measure different things: software checks (data integrity), coverage (how much the rules could label), and agreement with human-reviewed labels (the only accuracy-type measure).

## Limitations (also shown next to each result in the report)

- ClinicalTrials.gov is not a global or preclinical census. jRCT, ChiCTR, CTIS-only and undisclosed programmes are missing.
- 70% of active drug-type trials have no classified mechanism. Field-wide distributions are unknown.
- Open Targets scores are prioritisation evidence, not probabilities of success. PubMed counts measure volume, not quality.
- Registry status and dates reflect the snapshot date. A passed primary-completion date does not establish missing results.
- Many curated rules carried over from v1 have not had their sources re-verified. Their tier says so.

## Sharing

Shared as a portfolio prototype. Scientific review of the labels and programme notes is still pending (see `docs/REVIEW_TASKS.md`), so treat figures as exploratory.
