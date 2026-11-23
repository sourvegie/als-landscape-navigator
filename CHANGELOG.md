# Change log: v1 → v2.0.0 (classification version 2.0.0)

Repaired on 2026-11-20 from the cached snapshot retrieved 2026-11-18. The analysis date is fixed at 2026-11-18. The raw source files in `data/raw/` are byte-identical to the originals. v1 code and outputs are kept in `data/original_outputs/` for comparison.

Headline numbers (unique trials, active = registry status in snapshot):

| | v1 code re-run | v2 |
|---|---|---|
| Interventional trials | 724 | 724 |
| Trials in the drug landscape | 465 ("drug" = any drug-type intervention) | 455 (≥1 drug-type **test** intervention) |
| Active | 106 | 105 |
| Active with no classified mechanism | 79 (75%) | 73 (70%) |
| Targets flagged "open" | 26 of 30, of which 24 were unmapped | 0 "open" flags. One cautious "review" flag is possible only for a mapped gene with no trials; none occurs in this snapshot |
| Trials whose labels changed | — | 59 of the v1 drug trials (47 on mechanism) |

## 1. Target mapping and opportunity flags
- **Error.** v1 joined genes to trials only through five mechanism buckets. Any unmapped gene got `active_trials = 0` and could be tagged **open**. In the v1 run, 24 of the 26 "open" tags were on genes that had never been mapped at all.
- **Fix.** Every Open Targets gene now has an explicit state: `mapping not implemented`, `mapped: matching trials identified`, or `mapped: no matching trials identified in searched dataset`. Unmapped genes show "not assessed", never 0, and cannot be flagged. The word "open" is gone. A mapped gene with no trials gets only a "review" flag.
- **Error.** v1 gave every TDP-43-bucket trial to TARDBP and to STMN2. **Fix.** Direct relationships now come only from curated assets (`data/curation/asset_curation.csv`, with source and rationale). The bucket count is shown separately and labelled as a proxy. TARDBP: 3 direct active trials (VTx-002 and individualised TARDBP ASOs) vs 4 in the TDP-43 bucket (adds QRL-201, which targets STMN2).
- v1 counted the SOD1 and FUS buckets from description keywords. A blood-sampling study (NCT07567664) and an IVIG trial whose description mentions FUS (NCT07193953) were counted as SOD1/FUS programmes. Both are now excluded.

## 2. Intervention classification
- **Roles from arm structure.** Each intervention now gets a role: experimental, experimental (inferred), unresolved, background, comparator, control-arm treatment, placebo/sham, or non-therapeutic. Only test roles drive trial labels. Background and comparator treatments are listed in separate columns.
- **NCT04454840 (plasma + riluzole).** v1 labelled the combined plasma intervention Excitotoxicity / Small molecule because its name contains "riluzole". v2 splits off the standard-of-care component as background. The plasma component is Unclassified / Other biologic, marked as a registry-type fallback. The stand-alone riluzole arm is the **comparator**; v1 had treated it as background.
- **GENETIC ≠ AAV.** GENETIC and COMBINATION_PRODUCT entries without other evidence stay **Unresolved** (e.g. INS1202). Generic "gene therapy" or "plasmid" wording maps to *Gene therapy (non-AAV/unspecified)*. Only explicit AAV wording maps to AAV.
- **CNM-Au8.** Its registry composition is elemental gold nanocrystals in buffered water, so v2 adds a new modality, *Inorganic nanocrystal*. Mechanism stays Mitochondria & cell stress, marked as sponsor-proposed (NADH oxidation).
- **Non-therapeutic entries** (blood sampling, PET/SPECT tracers including 89Zr-AP-101, ultrasound microspheres) are excluded from the therapeutic landscape.
- **Placebo detection.** "X or placebo" (blinded naming) is no longer dropped as placebo. "Placebos" (plural) is now recognised.
- **Supportive co-medication** (e.g. prophylactic methylprednisolone in the VTx-002 trial) is background.
- **Description text no longer assigns mechanisms.** Description-only keyword hits are stored as review *candidates*. v1 errors this removes:
  - "calpain" matched "pain", so the calpain-2 ASO AMX0114 became Symptomatic;
  - creatine mentioned in the celecoxib description;
  - SOD1 mice mentioned in the triheptanoin description;
  - IL-2 mentioned in a glucose-water control description.
- **Every label now records its evidence tier and rule**: curated asset (source verified / registry-only / carried-over unverified v1 rule), override, name keyword, description keyword (modality only), registry-type fallback, unresolved.
- **Curated-table changes.** apilimod moved from TDP-43 biology to Proteostasis (source: Hung et al. 2023). The AMX0035 rule was split from TUDCA and phenylbutyrate monotherapy, and troriluzole from riluzole. ARO-SOD1 (siRNA), AMT-162, NTF001, XT-150, AMX0114 and individualised n-Lorem ASOs were added with registry or sponsor sources. The unverifiable `AS-202` pattern was removed.

## 3. Coverage statistics and conclusions
- v2 reports separately: intervention-level mechanism coverage by tier, modality coverage (explicit vs fallback vs unresolved), trial-level coverage (all / partial / none classified), and the size of the mapped target subset.
- **Error.** v1 reported "Genetically defined programmes = 8% … most activity still targets broad downstream pathways" using all active trials as the denominator, 75% of which were unclassified. **Fix.** The statement is now restricted to the classified subset. The whole-field share is given only as a bound (12%–82% in this snapshot).
- The "Gap map" is renamed **Classified trial activity by mechanism and modality**. The unclassified share is shown prominently above it, and an inline caveat says that empty cells are not white space.
- Counting units are explicit: unique trials, intervention records, unique curated assets. Cells count each trial once. The report explains how many trials appear in several cells.
- KPI labels now name supplements and the other registry types included.

## 4. Clinical milestones
- "Upcoming results" is renamed **Registered primary-completion milestones**, with an inline note that a primary-completion date marks the end of primary-outcome data collection and says nothing about when results will be announced.
- Each date keeps its estimated/actual type and its precision. Month-only dates cover the whole month, so they are never "passed" before the month ends (v1 parsed "2026-11" as 1 November).
- The claim that trials "show no update" is removed. Passed dates are listed with last-update-posted date and results-posted flag, and no inference is drawn from them.
- Status is described as "status in snapshot". The report explains why Phase 1/2 counts as Phase 2+.

## 5. Validation
- No manual labels were invented. `data/label_sample.csv` is unchanged. Its `auto_*` columns are from an early v1 run and are stale.
- A held-out/dev split (`data/validation/eval_split.csv`, 24/26, hash-based) was frozen **before** any review-informed change. No rule was changed using the label sample.
- `data/validation/label_review_pack.csv`: registry link, arm context, roles, current predictions with rules and tiers, AI-proposed labels with reasoning and sources (`proposed_*`, kept separate), open questions, and empty `human_*` fields for review.
- `--evaluate` scores only human-approved rows, separately for mechanism and modality, each with its own denominator. A row with only one axis reviewed is never an error on the other axis. Held-out and dev are reported separately. Results carry the classifier version and a label fingerprint and are shown as **stale** if either changes.
- The dashboard shows **"Classification accuracy not yet measured"**. Software checks, coverage and agreement are in separate panels.

## 6. Regression checks
`tests/test_navigator.py` has 31 tests, all passing. They cover: unmapped targets never flagged; the plasma + riluzole case; unknown GENETIC entries staying unresolved; CNM-Au8 source-supported classification; experimental vs background vs comparator roles; unique trial counting across repeated cohorts; month-only dates; separate per-axis denominators and stale-result invalidation; unavailable data shown as unavailable rather than zero; byte-identical outputs across rebuilds; and known-answer regression checks (labelled as not accuracy). `tests/mutation_smoke.py` re-introduces four v1 defects and confirms that each one is caught.

## 7. Reproducibility
- **Error.** The v1 saved HTML/CSV did not match the v1 code: 7 trials carried "Excitotoxicity" from an older rule that matched riluzole in descriptions. **Fix.** All outputs are now regenerated together.
- The analysis date is explicit (default: the ClinicalTrials.gov retrieval date; override with `--analysis-date`). Nothing uses today's date.
- `output/run_manifest.json` records retrieval dates, queries and limits, analysis date, classification version, validation state, build timestamp, and SHA-256 hashes of inputs and analytical outputs. `SOURCE_DATE_EPOCH` makes the build timestamp reproducible too.
- The synthetic preview moved to `test_preview/` and carries a "SYNTHETIC TEST FIXTURE – NOT PUBLIC DATA" banner.

## 8. Interpretation
Caveats now sit beside the results they qualify:
- Open Targets scores are prioritisation evidence;
- PubMed counts measure volume, not quality;
- registry activity does not show efficacy;
- missing mapped trials do not show absence of development;
- ClinicalTrials.gov is not a global or preclinical census.

Provisional, source-linked assessments of QRL-201, AP-101 and CNM-Au8 are in `docs/programme_assessments.md`.
