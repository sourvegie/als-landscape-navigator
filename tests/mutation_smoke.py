"""Mutation smoke check: re-introduce four v1 defects in memory and confirm the matching regression
test fails for each. Run: python3 tests/mutation_smoke.py  (expects every count > 0)."""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import als_navigator as N  # noqa: E402
import test_navigator as T  # noqa: E402


def failures(names):
    suite = unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromName(n, T) for n in names])
    with open(os.devnull, "w") as sink:
        r = unittest.TextTestRunner(stream=sink).run(suite)
    return len(r.failures) + len(r.errors)


results = []
N.TYPE_FALLBACK_MODALITY["GENETIC"] = "AAV gene therapy"
results.append(("GENETIC defaults to AAV", failures(["TestGeneticNotAAV.test_unknown_genetic_unresolved"])))
del N.TYPE_FALLBACK_MODALITY["GENETIC"]

orig = N.milestone_position
N.milestone_position = lambda d, a, h: "passed" if d and d["start"] < a else orig(d, a, h)
results.append(("month-only date overdue from day 1", failures(
    ["TestMonthOnlyDates.test_month_precision_not_overdue_before_month_end"])))
N.milestone_position = orig

orig_split = N.split_background_components
N.split_background_components = lambda n: (n, [])
results.append(("plasma+riluzole not split", failures(["TestPlasmaPlusRiluzole.test_snapshot_NCT04454840"])))
N.split_background_components = orig_split

orig_roles = N.resolve_roles


def all_experimental(ivs, arms):
    out = orig_roles(ivs, arms)
    for o in out:
        if o["role"] != N.ROLE_PLACEBO:
            o["role"] = N.ROLE_EXPERIMENTAL
    return out


N.resolve_roles = all_experimental
results.append(("background treated as experimental", failures(["TestRoles.test_background_vs_experimental"])))
N.resolve_roles = orig_roles

ok = all(n > 0 for _, n in results)
for name, n in results:
    print(f"  {'caught' if n else 'MISSED'}: {name} ({n} failing test(s))")
sys.exit(0 if ok else 1)
