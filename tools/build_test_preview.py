#!/usr/bin/env python3
"""Build the SYNTHETIC test-fixture preview, kept apart from the public-data report.

Output goes to test_preview/ only. Nothing here reads or writes data/ or output/.
"""
import datetime as dt
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))
import als_navigator as N  # noqa: E402
from fixtures import synthetic_dataset  # noqa: E402

if __name__ == "__main__":
    ct, ot, pm = synthetic_dataset()
    base = os.path.join(ROOT, "test_preview")
    res = N.build(ct, ot, pm, os.path.join(base, "data_SYNTHETIC"), base, dt.date(2026, 11, 18),
                  "fixed for the synthetic fixture", fixture=True, write_pack=False,
                  html_name="preview_SYNTHETIC_TEST_FIXTURE.html")
    print("Synthetic preview:", res["html"])
