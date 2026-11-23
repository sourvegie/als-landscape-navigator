"""Regression checks for the ALS Landscape Navigator.

Run:  python3 -m unittest discover -s tests -v

These tests guard specific analytical behaviours. Known-answer checks are regression checks only;
passing them is NOT evidence of overall classification accuracy.
"""
import copy
import csv
import datetime as dt
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import als_navigator as N  # noqa: E402
from fixtures import synthetic_dataset, study  # noqa: E402

SNAPSHOT = os.path.join(ROOT, "data", "raw", "ctgov_studies.json")
AD = dt.date(2026, 11, 18)


def run_synthetic(ct=None, ot=None, pm=None, analysis_date=AD):
    ct0, ot0, pm0 = synthetic_dataset()
    ct = ct if ct is not None else ct0
    ot = ot if ot is not None else ot0
    pm = pm0 if pm is None else pm
    cur = N.load_curation()
    trials = N.normalise_trials(ct, cur, {})
    a = N.analyse(ot, pm, trials, cur, analysis_date)
    return trials, a


def by_id(trials):
    return {t["nct_id"]: t for t in trials}


class SnapshotMixin:
    @classmethod
    def setUpClass(cls):
        with open(SNAPSHOT) as f:
            cls.ct = json.load(f)
        cls.cur = N.load_curation()
        cls.trials = N.normalise_trials(cls.ct, cls.cur, {})
        cls.T = by_id(cls.trials)


# 1 -------------------------------------------------------------------------------------------
class TestTargetMappingStates(unittest.TestCase):
    def test_unmapped_targets_never_flagged(self):
        _, a = run_synthetic()
        t = {x["symbol"]: x for x in a["targets"]}
        self.assertEqual(t["UNMAPPEDGENE"]["mapping_state"], "mapping not implemented")
        self.assertEqual(t["UNMAPPEDGENE"]["flag"], "")
        self.assertIsNone(t["UNMAPPEDGENE"]["direct_active"])  # not assessed, NOT zero
        for x in a["targets"]:
            if x["mapping_state"] == "mapping not implemented":
                self.assertEqual(x["flag"], "", x["symbol"])

    def test_three_explicit_states(self):
        _, a = run_synthetic()
        t = {x["symbol"]: x for x in a["targets"]}
        self.assertEqual(t["TARDBP"]["mapping_state"], "mapped: matching trials identified")
        self.assertEqual(t["ATXN2"]["mapping_state"], "mapped: no matching trials identified in searched dataset")
        self.assertTrue(t["ATXN2"]["flag"].startswith("review"))  # cautious wording, never "open"
        self.assertNotIn("open", t["ATXN2"]["flag"])

    def test_bucket_proxy_is_separate_from_direct(self):
        ct, ot, pm = synthetic_dataset()
        # a TDP-43-bucket trial that is NOT a direct TARDBP relationship (STMN2 restoration)
        ct["studies"].append(study("NCT99000010", [("A", "EXPERIMENTAL")], [("DRUG", "QRL-201", ["A"], "")]))
        _, a = run_synthetic(ct, ot, pm)
        t = {x["symbol"]: x for x in a["targets"]}["TARDBP"]
        self.assertEqual(t["direct_active"], 1)          # VTx-002 only
        self.assertEqual(t["bucket"], "TDP-43 biology")
        self.assertEqual(t["bucket_active"], 2)          # VTx-002 + QRL-201, labelled as proxy

    def test_snapshot_unmapped_not_flagged(self):
        with open(SNAPSHOT) as f:
            ct = json.load(f)
        with open(os.path.join(ROOT, "data", "raw", "opentargets.json")) as f:
            ot = json.load(f)
        with open(os.path.join(ROOT, "data", "raw", "pubmed_counts.json")) as f:
            pm = json.load(f)
        cur = N.load_curation()
        a = N.analyse(ot, pm, N.normalise_trials(ct, cur, {}), cur, AD)
        unm = [x for x in a["targets"] if x["mapping_state"] == "mapping not implemented"]
        self.assertGreater(len(unm), 0)
        self.assertTrue(all(x["flag"] == "" and x["direct_all"] is None for x in unm))


# 2 -------------------------------------------------------------------------------------------
class TestPlasmaPlusRiluzole(SnapshotMixin, unittest.TestCase):
    def test_snapshot_NCT04454840(self):
        t = self.T["NCT04454840"]
        self.assertEqual(len(t["test"]), 1)
        r = t["test"][0]
        self.assertNotEqual(r["mechanism"], "Excitotoxicity")
        self.assertNotEqual(r["modality"], "Small molecule")
        self.assertEqual(r["classified_text"], "Plasma from healthy young people treatment")
        self.assertIn("Riluzole", r["background_components"])
        roles = {x["name"]: x["role"] for x in t["records"]}
        self.assertEqual(roles["Riluzole"], N.ROLE_COMPARATOR)  # stand-alone riluzole arm is the comparator

    def test_synthetic_combination(self):
        trials, _ = run_synthetic()
        t = by_id(trials)["NCT99000006"]
        self.assertEqual([r["modality"] for r in t["test"]], ["Other biologic"])
        self.assertEqual(t["test"][0]["modality_tier"], N.TIER_TYPE)  # labelled as fallback, not explicit
        self.assertNotIn("Excitotoxicity", t["mechanisms"])


# 3 -------------------------------------------------------------------------------------------
class TestGeneticNotAAV(unittest.TestCase):
    def test_unknown_genetic_unresolved(self):
        trials, _ = run_synthetic()
        r = by_id(trials)["NCT99000003"]["test"][0]
        self.assertEqual(r["modality"], "Unresolved")
        self.assertEqual(r["modality_tier"], N.TIER_NONE)

    def test_generic_gene_therapy_keyword_not_aav(self):
        res = N.classify_intervention({"name": "Gene therapy product GT-1", "type": "GENETIC"}, N.load_curation(), {})
        self.assertEqual(res["modality"], "Gene therapy (non-AAV/unspecified)")

    def test_sampling_genetic_excluded(self):
        trials, _ = run_synthetic()
        t = by_id(trials)["NCT99000005"]
        self.assertFalse(t["is_drug"])
        self.assertEqual(t["records"][0]["role"], N.ROLE_NONTHERAPEUTIC)


class TestGeneticSnapshot(SnapshotMixin, unittest.TestCase):
    def test_ins1202_unresolved(self):
        r = self.T["NCT07290062"]["test"][0]
        self.assertEqual(r["modality"], "Unresolved")


# 4 -------------------------------------------------------------------------------------------
class TestCNMAu8(SnapshotMixin, unittest.TestCase):
    def test_all_records(self):
        recs = [r for t in self.trials for r in t["test"] if r.get("asset_id") == "CNM-Au8"]
        self.assertGreaterEqual(len(recs), 4)
        for r in recs:
            self.assertEqual(r["modality"], "Inorganic nanocrystal")
            self.assertNotEqual(r["modality"], "Other biologic")
            self.assertEqual(r["modality_tier"], N.TIER_CURATED_VERIFIED)
            self.assertIn("clinicaltrials.gov", r["source_urls"])


# 5 -------------------------------------------------------------------------------------------
class TestRoles(unittest.TestCase):
    def test_background_vs_experimental(self):
        trials, _ = run_synthetic()
        t = by_id(trials)["NCT99000001"]
        roles = {x["name"]: x["role"] for x in t["records"]}
        self.assertEqual(roles["Masitinib"], N.ROLE_EXPERIMENTAL)
        self.assertEqual(roles["Riluzole"], N.ROLE_BACKGROUND)
        self.assertEqual(roles["Placebo"], N.ROLE_PLACEBO)
        self.assertEqual(t["mechanisms"], ["Neuroinflammation & immune"])
        self.assertIn("Riluzole", t["background"])

    def test_supportive_premedication_is_background(self):
        trials, _ = run_synthetic()
        t = by_id(trials)["NCT99000009"]
        self.assertEqual([r["name"] for r in t["test"]], ["VTx-002"])

    def test_blinded_or_placebo_name_is_test(self):
        trials, _ = run_synthetic()
        t = by_id(trials)["NCT99000008"]
        self.assertTrue(t["is_drug"])
        self.assertEqual(t["mechanisms"], ["Neuroinflammation & immune"])

    def test_placebo_to_match_is_placebo(self):
        self.assertTrue(N.is_placebo_name("Placebo to match riluzole"))
        self.assertTrue(N.is_placebo_name("Matching placebo 25 mg and Matching placebo 50 mg"))
        self.assertFalse(N.is_placebo_name("PLL001 or placebo daily subcutaneous injections"))

    def test_description_keyword_is_candidate_only(self):
        trials, _ = run_synthetic()
        r = by_id(trials)["NCT99000004"]["test"][0]
        self.assertEqual(r["mechanism"], "Unclassified")
        self.assertEqual(r["mechanism_candidate"], "SOD1")


# 6 -------------------------------------------------------------------------------------------
class TestUniqueTrialCounting(unittest.TestCase):
    def test_repeated_cohorts_count_once(self):
        trials, a = run_synthetic()
        t = by_id(trials)["NCT99000002"]
        self.assertEqual(len(t["test"]), 3)                    # three intervention records
        self.assertEqual(a["grid_active"][("SOD1", "ASO")], 1)  # one unique trial
        self.assertEqual(a["mech_active"]["SOD1"], 1)
        self.assertEqual(a["coverage"]["unique_assets_active"].count("tofersen"), 1)


# 7 -------------------------------------------------------------------------------------------
class TestMonthOnlyDates(unittest.TestCase):
    def test_month_precision_not_overdue_before_month_end(self):
        d = N.parse_registry_date("2026-11")
        self.assertEqual(d["precision"], "month")
        horizon = dt.date(2028, 12, 31)
        self.assertEqual(N.milestone_position(d, dt.date(2026, 11, 18), horizon), "in window")
        self.assertEqual(N.milestone_position(d, dt.date(2026, 11, 30), horizon), "in window")
        self.assertEqual(N.milestone_position(d, dt.date(2026, 12, 1), horizon), "passed")

    def test_day_precision(self):
        d = N.parse_registry_date("2026-11-17")
        self.assertEqual(N.milestone_position(d, dt.date(2026, 11, 18), dt.date(2028, 12, 31)), "passed")

    def test_pipeline_uses_analysis_date(self):
        _, a = run_synthetic()
        ids = [t["nct_id"] for t in a["milestones"]]
        self.assertIn("NCT99000002", ids)          # 2026-11 is in window on 2026-11-18
        self.assertNotIn("NCT99000008", [t["nct_id"] for t in a["milestones"]])  # 2026-10 passed
        self.assertIn("NCT99000008", [t["nct_id"] for t in a["passed"]])


# 8 -------------------------------------------------------------------------------------------
class TestSeparateDenominators(unittest.TestCase):
    def test_partial_rows_not_errors(self):
        trials, _ = run_synthetic()
        split = {"NCT99000001": "heldout", "NCT99000002": "heldout", "NCT99000003": "heldout",
                 "NCT99000005": "heldout", "NCT99000007": "dev"}
        rows = [
            # mechanism only -> counts for mechanism, ignored for modality
            {"nct_id": "NCT99000001", "human_mechanism": "Neuroinflammation & immune", "human_modality": "",
             "human_review_status": "approved"},
            # modality only -> counts for modality, ignored for mechanism
            {"nct_id": "NCT99000002", "human_mechanism": "", "human_modality": "ASO", "human_review_status": "approved"},
            # both, but NOT approved -> ignored entirely
            {"nct_id": "NCT99000003", "human_mechanism": "SOD1", "human_modality": "AAV gene therapy",
             "human_review_status": ""},
            # N/A means "no therapeutic intervention": matches an empty prediction
            {"nct_id": "NCT99000005", "human_mechanism": "N/A", "human_modality": "N/A",
             "human_review_status": "corrected"},
            {"nct_id": "NCT99000007", "human_mechanism": "Mitochondria & cell stress",
             "human_modality": "Inorganic nanocrystal", "human_review_status": "approved"},
        ]
        res = N.evaluate_pack(rows, trials, split)
        h = res["heldout"]
        self.assertEqual(h["mechanism"]["n"], 2)
        self.assertEqual(h["modality"]["n"], 2)
        self.assertEqual(h["mechanism"]["agree"], 2)
        self.assertEqual(h["modality"]["agree"], 2)
        self.assertEqual(res["dev"]["mechanism"]["n"], 1)   # dev reported separately, never pooled

    def test_status_not_measured_without_results(self):
        trials, _ = run_synthetic()
        with tempfile.TemporaryDirectory() as d:
            vs = N.validation_status(d, trials)
            self.assertEqual(vs["headline"], "Classification accuracy not yet measured")

    def test_stale_results_invalidated(self):
        trials, _ = run_synthetic()
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "validation_results.json"), "w") as f:
                json.dump({"classification_version": "1.0.0", "trials_fingerprint": "old",
                           "results": {"heldout": {"mechanism": {"n": 10}, "modality": {"n": 10}}}}, f)
            vs = N.validation_status(d, trials)
            self.assertEqual(vs["state"], "stale")
            self.assertIn("not yet measured", vs["headline"])


# 9 -------------------------------------------------------------------------------------------
class TestUnavailableIsNotZero(unittest.TestCase):
    def _build(self, ot=None, pm=None):
        ct, ot0, pm0 = synthetic_dataset()
        d = tempfile.mkdtemp()
        res = N.build(ct, ot0 if ot is None else ot, pm0 if pm is None else pm, os.path.join(d, "data"),
                      os.path.join(d, "out"), AD, "test", fixture=True, write_pack=False)
        with open(res["html"]) as f:
            return res, f.read(), d

    def test_pubmed_missing_symbol_is_not_retrieved(self):
        res, page, _ = self._build()
        t = {x["symbol"]: x for x in res["analysis"]["targets"]}
        self.assertIsNone(t["ATXN2"]["pubmed_total"])
        self.assertEqual(t["ATXN2"]["pubmed_status"], "not retrieved")
        self.assertEqual(t["TARDBP"]["pubmed_total"], 0)       # a real retrieved zero stays zero
        self.assertIn("not retrieved", page)

    def test_pubmed_unavailable(self):
        res, page, _ = self._build(pm={"fetched_at": None, "counts": None})
        self.assertTrue(all(x["pubmed_status"] == "unavailable" for x in res["analysis"]["targets"]))
        self.assertIn("unavailable", page)

    def test_open_targets_unavailable(self):
        res, page, _ = self._build(ot={"fetched_at": None, "rows": []})
        self.assertIn("Open Targets data unavailable", page)

    def test_missing_genetic_score_not_zero(self):
        res, page, d = self._build()
        with open(os.path.join(d, "data", "targets.csv")) as f:
            rows = {r["symbol"]: r for r in csv.DictReader(f)}
        self.assertEqual(rows["NOGENETIC"]["genetic_association_score"], "")
        self.assertEqual(rows["UNMAPPEDGENE"]["direct_active_trials"], "not assessed")
        self.assertIn("not reported", page)


# 10 ------------------------------------------------------------------------------------------
class TestReproducibility(unittest.TestCase):
    def test_identical_inputs_identical_outputs(self):
        ct, ot, pm = synthetic_dataset()
        outs = []
        for epoch in ("1795132800", "1795219200"):  # different build times
            os.environ["SOURCE_DATE_EPOCH"] = epoch
            d = tempfile.mkdtemp()
            N.build(copy.deepcopy(ct), ot, pm, os.path.join(d, "data"), os.path.join(d, "out"), AD, "test",
                    fixture=True, write_pack=False)
            files = {}
            for name in ("trials.csv", "interventions.csv", "targets.csv", "milestones.csv",
                         "unclassified_interventions.csv"):
                with open(os.path.join(d, "data", name), "rb") as f:
                    files[name] = f.read()
            outs.append(files)
        os.environ.pop("SOURCE_DATE_EPOCH", None)
        self.assertEqual(outs[0], outs[1])

    def test_analysis_date_controls_results(self):
        _, a1 = run_synthetic(analysis_date=dt.date(2026, 11, 18))
        _, a2 = run_synthetic(analysis_date=dt.date(2026, 12, 1))
        self.assertNotEqual(a1["ms_summary"], a2["ms_summary"])


# Known-answer regression checks (NOT accuracy) ---------------------------------------------------
class TestKnownAnswerRegression(SnapshotMixin, unittest.TestCase):
    def test_known_answers(self):
        res = N.regression_answer_checks(self.trials)
        failed = [c for c in res if c["status"] == "FAIL"]
        self.assertEqual(failed, [], failed)

    def test_curation_file_well_formed(self):
        rules = N.load_curation()
        self.assertTrue(all(r["source_status"] for r in rules))
        verified = [r for r in rules if r["tier"] == N.TIER_CURATED_VERIFIED]
        self.assertTrue(all(r["source_urls"].startswith("http") for r in verified))


if __name__ == "__main__":
    unittest.main(verbosity=2)
