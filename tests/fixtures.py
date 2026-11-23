"""Synthetic ClinicalTrials.gov-shaped records for regression tests and the TEST-FIXTURE preview.

Every record here is invented. Nothing in this file describes a real trial; NCT numbers use the
NCT99xxxxxx range so they cannot be mistaken for registry IDs used in the public snapshot.
"""


def study(nct, arms, interventions, status="RECRUITING", phases=("PHASE2",), pcd=None, pcd_type="ESTIMATED",
          sponsor="Synthetic Sponsor", sponsor_class="INDUSTRY", countries=("United States",), title=None):
    """arms: list of (label, type); interventions: list of (type, name, [arm labels], description)."""
    return {
        "protocolSection": {
            "identificationModule": {"nctId": nct, "briefTitle": title or f"Synthetic study {nct}"},
            "statusModule": {"overallStatus": status,
                             "primaryCompletionDateStruct": ({"date": pcd, "type": pcd_type} if pcd else {}),
                             "lastUpdatePostDateStruct": {"date": "2026-01-01", "type": "ACTUAL"}},
            "designModule": {"studyType": "INTERVENTIONAL", "phases": list(phases),
                             "enrollmentInfo": {"count": 40}},
            "sponsorCollaboratorsModule": {"leadSponsor": {"name": sponsor, "class": sponsor_class}},
            "armsInterventionsModule": {
                "armGroups": [{"label": lb, "type": ty} for lb, ty in arms],
                "interventions": [{"type": ty, "name": nm, "armGroupLabels": list(al), "description": ds}
                                  for ty, nm, al, ds in interventions],
            },
            "contactsLocationsModule": {"locations": [{"country": c} for c in countries]},
        },
        "hasResults": False,
    }


def ot_row(symbol, score=0.8, genetic=0.9):
    dts = [] if genetic is None else [{"id": "genetic_association", "score": genetic}]
    return {"score": score, "target": {"id": f"ENSG_SYN_{symbol}", "approvedSymbol": symbol,
                                       "approvedName": f"synthetic {symbol}"}, "datatypeScores": dts}


def synthetic_dataset():
    S = []
    # Add-on design: riluzole in both arms is background; X is experimental.
    S.append(study("NCT99000001", [("Drug X", "EXPERIMENTAL"), ("Placebo", "PLACEBO_COMPARATOR")],
                   [("DRUG", "Masitinib", ["Drug X"], "oral"),
                    ("DRUG", "Riluzole", ["Drug X", "Placebo"], "background therapy in both arms"),
                    ("DRUG", "Placebo", ["Placebo"], "")], pcd="2027-03"))
    # Repeated cohorts of the same asset: one trial.
    S.append(study("NCT99000002", [("Cohort A", "EXPERIMENTAL"), ("Cohort B", "EXPERIMENTAL"), ("OLE", "EXPERIMENTAL")],
                   [("DRUG", "Tofersen cohort A", ["Cohort A"], ""), ("DRUG", "Tofersen cohort B", ["Cohort B"], ""),
                    ("DRUG", "Tofersen", ["OLE"], "")], phases=("PHASE1", "PHASE2"), pcd="2026-11"))
    # Unknown GENETIC product.
    S.append(study("NCT99000003", [("Gene", "EXPERIMENTAL")],
                   [("GENETIC", "ZZX-901", ["Gene"], "Suspension for injection.")], phases=("PHASE1",)))
    # Description mentions SOD1 mouse model: must NOT become SOD1.
    S.append(study("NCT99000004", [("T", "EXPERIMENTAL"), ("P", "PLACEBO_COMPARATOR")],
                   [("DRUG", "Compound Q-77", ["T"], "Improved survival in the SOD1 G93A mouse; anti-inflammatory"),
                    ("DRUG", "Matching placebo", ["P"], "")], status="COMPLETED", pcd="2020-05-01", pcd_type="ACTUAL"))
    # Sampling-only study registered as GENETIC.
    S.append(study("NCT99000005", [("All", "OTHER")],
                   [("GENETIC", "Blood sample for genetic analysis", ["All"], "SOD1, C9orf72, FUS panel")],
                   phases=()))
    # Combined name with standard of care, plus a riluzole comparator arm.
    S.append(study("NCT99000006", [("Combo", "EXPERIMENTAL"), ("SoC", "ACTIVE_COMPARATOR")],
                   [("BIOLOGICAL", "Donor plasma infusion + Riluzole", ["Combo"], "plasma every 2 weeks"),
                    ("DRUG", "Riluzole", ["SoC"], "")], phases=("EARLY_PHASE1",)))
    # Gold nanocrystal suspension.
    S.append(study("NCT99000007", [("Au", "EXPERIMENTAL"), ("P", "PLACEBO_COMPARATOR")],
                   [("DRUG", "CNM-Au8", ["Au"], "aqueous suspension of gold nanocrystals"),
                    ("DRUG", "Placebo", ["P"], "")], phases=("PHASE3",), pcd="2028-06"))
    # Blinded naming: "X or placebo".
    S.append(study("NCT99000008", [("Blinded", "EXPERIMENTAL")],
                   [("DRUG", "Ibudilast or placebo capsules", ["Blinded"], "")], pcd="2026-10"))
    # Supportive co-medication alongside a vectorised antibody.
    S.append(study("NCT99000009", [("Low", "EXPERIMENTAL"), ("High", "EXPERIMENTAL")],
                   [("GENETIC", "VTx-002", ["Low", "High"], "gene therapy"),
                    ("DRUG", "Prophylactic corticosteroid premedication", ["Low", "High"], "")],
                   phases=("PHASE1", "PHASE2"), pcd="2027-10-15"))
    ct = {"fetched_at": "2026-11-18T00:00:00", "totalCount": len(S), "studies": S, "_fixture": True}
    ot = {"fetched_at": "2026-11-18T00:00:00", "disease": {"id": "MONDO_SYNTHETIC", "name": "synthetic"},
          "rows": [ot_row("SOD1"), ot_row("TARDBP"), ot_row("ATXN2"), ot_row("UNMAPPEDGENE", genetic=0.99),
                   ot_row("NOGENETIC", genetic=None)]}
    pm = {"fetched_at": "2026-11-18T00:00:00", "window": [2022, 2026],
          "counts": {"SOD1": {"total": 10, "recent": 3}, "TARDBP": {"total": 0, "recent": 0}}}
    return ct, ot, pm
