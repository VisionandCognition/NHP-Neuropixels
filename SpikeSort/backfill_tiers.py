"""Add usable / su_like / tier columns to units.csv of every finished session (idempotent).

    python backfill_tiers.py            # all sessions under config.OUTPUT_ROOT
Writes SpikeSorting/_tiers_summary.csv (units per tier per session).
"""
import csv

import pandas as pd

import config
import curation


def main():
    root = config.OUTPUT_ROOT
    rows = []
    for f in sorted(root.glob("*/*/*/units.csv")) + sorted(root.glob("*/*/*/*/units.csv")):
        q = pd.read_csv(f)
        if "good" not in q:
            continue
        q = q.drop(columns=[c for c in ("usable", "su_like", "tier") if c in q])
        q["peak_channel_id"] = q["peak_channel_id"].astype(str)
        curation.add_tiers(q)
        q.to_csv(f, index=False)
        rows.append(dict(session=str(f.parent.relative_to(root)), n_units=len(q), good=int(q["good"].sum()),
                         usable=int(q["usable"].sum()), su_like=int(q["su_like"].sum())))
        print(rows[-1])
    with open(root / "_tiers_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["session", "n_units", "good", "usable", "su_like"])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
