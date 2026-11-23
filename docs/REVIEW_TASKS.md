# Outstanding human-review tasks and scientific questions

Prepared 2026-11-20. Nothing listed here has been reviewed by a human scientist.

## A. Validation (blocking any accuracy claim)
1. **Review all 50 rows** in `data/validation/label_review_pack.csv`. Fill `human_*` fields; mark each `approved`, `corrected` or `unresolved`. Treat `proposed_*` as drafts: several identities are flagged "from model background knowledge; not source-checked" (e.g. GSK1223249 = ozanezumab, VM202 = HGF plasmid, sNN0029 = VEGF, prosetin target).
2. **Review the held-out subset (24 rows) first and do not tune rules on it.** Rule gaps already noticed in **dev** rows (allowed to fix after review): `MSC_NTF` underscore not matched by the NurOwn rule; `AVP-923` and `dronabinol` not in the symptomatic list; "progenitor cells" and "HYNR-CS" not recognised as cell therapy; antibodies registered as DRUG (ANX005) fall back to Small molecule.
3. Decide the labelling convention for multi-agent or platform trials (e.g. ROAR-DIGAP NCT06429059, HEALEY master protocol NCT04297683): trial-level union, or one row per regimen?
4. Decide whether healthy-volunteer PK studies indexed to ALS (e.g. NCT05797701) and imaging studies belong in the landscape.

## B. Curated asset table (`data/curation/asset_curation.csv`)
5. Re-verify the 19 rules marked `carried_over_unverified` (edaravone, riluzole, AMX0035, masitinib, ibudilast, verdiperstat, zilucoplan, IL-2, Treg, eIF2B activators, troponin activators, proteostasis compounds, NurOwn, methylcobalamin, symptomatic list and others). Add a source URL to each.
6. Confirm three judgement calls:
   - apilimod moved from TDP-43 biology to Proteostasis & autophagy;
   - AMX0114 (calpain-2 ASO) placed in Other neuroprotective;
   - pridopidine → SIGMAR1 as a direct target, given the reported D2/D3 binding.
7. **Zilucoplan modality.** It is a synthetic macrocyclic peptide, currently labelled "Other biologic". Add a *Peptide* category?
8. **Category gaps to decide.** Herbal and Kampo/TCM preparations (TJ-68, C. cicadae), microbiome products (MTP-101C), and inorganic salts (lithium) filed under Small molecule.
9. Identify INS1202 (NCT07290062). It is registered as GENETIC, enrols SOD1 carriers or patients with no known mutation, and is currently Unresolved.

## C. Role resolution
10. Check the five `unresolved role` records and the 60 intervention records flagged `no EXPERIMENTAL arm in registry` in `data/interventions.csv`, especially:
    - NCT00403104 (ONO-2506PO linked only to a placebo arm);
    - NCT00647296 (dexpramipexole linked only to a placebo arm);
    - NCT06391645 (NGF linked to both experimental and comparator arms).

## D. Target mapping
11. Extend direct-target curation beyond the 6 mapped Open Targets genes only where a sourced asset exists (e.g. any TBK1-, OPTN- or UNC13A-directed clinical assets). Until then those genes stay "not assessed".
12. Decide whether a "review" flag should also apply to mapped genes whose only trials are inactive (currently C9orf72 and ATXN2). The current rule flags only genes with no trials at all.

## E. Programme assessments (`docs/programme_assessments.md`)
13. **QRL-201.** How was tissue target engagement measured, and what do the prespecified intention-to-treat analyses of the dose-range-finding phase show?
14. **AP-101.** Read the results posted on ClinicalTrials.gov for NCT05039099. Check how the 12-month survival/NIV analysis handled crossover after week 24.
15. **CNM-Au8.** Is there a prespecified, randomised NfL effect in either RCT (JAMA HEALEY paper; posted results for NCT04098406 and NCT04414345)?
16. A scientist should sign off on, or rewrite, each "key uncertainty" before any external use.

## F. Data and scope
17. The snapshot's registry statuses are as of 2026-11-18. Some are known to be out of date relative to later news (e.g. VHB937 was reportedly discontinued in Sept 2026 while its record still shows active, not recruiting). Decide whether to note such cases, or to rely on registry status only.
18. Visual check: the report was rendered and inspected in headless Chrome in the build environment. A human should also check it in their own browser, including dark mode.
