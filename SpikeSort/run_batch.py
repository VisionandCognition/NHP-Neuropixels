"""Sort every sortable (monkey, day, node, experiment) found under Data_collection.

    python run_batch.py --list                    # show what would be sorted, with sizes
    python run_batch.py                           # sort everything not yet DONE (resumable)
    python run_batch.py --monkey m032 --day 20260312 --max-sessions 3
    python run_batch.py --force                   # redo already finished sessions

Safe to interrupt and restart: a session counts as finished only if its DONE.json exists.
Each session logs to <session out dir>/sort.log; a one-line-per-attempt history is appended
to _manifest.csv. Failures are recorded and the batch continues.
"""
from __future__ import annotations

import argparse
import csv
import logging
import traceback
from datetime import datetime

import config
from find_recordings import exp_number, scan_all
import sort_session as ss

MANIFEST_FIELDS = ["time", "monkey", "day", "node", "experiment", "status", "n_units", "n_good", "error"]


def sort_units(rows):
    """Group sortable recording rows into sort units (monkey, day, node, experiment)."""
    units = {}
    for r in rows:
        if int(r["sortable"] or 0):
            key = (r["monkey"], r["day"], r["node"], exp_number(r))
            units.setdefault(key, []).append(r)
    return dict(sorted(units.items()))


def log_manifest(**kw):
    config.MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    new = not config.MANIFEST.exists()
    with open(str(config.MANIFEST), "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        if new:
            w.writeheader()
        w.writerow(dict({k: "" for k in MANIFEST_FIELDS}, time=datetime.now().isoformat(timespec="seconds"), **kw))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--monkey")
    ap.add_argument("--day")
    ap.add_argument("--max-sessions", type=int)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    units = sort_units(scan_all(config.DATA_COLLECTION))
    todo = []
    for (monkey, day, node, exp_n), recs in units.items():
        if (args.monkey and monkey != args.monkey) or (args.day and day != args.day):
            continue
        done = (ss.out_dir_for(monkey, day, exp_n, node) / "DONE.json").exists()
        gb = sum(int(r["ap_bytes"]) for r in recs) / 1e9
        if args.list:
            print("%-5s %s %-12s %6.1f GB %s" % (monkey, day, ss.session_name(exp_n, node), gb,
                                                 "DONE" if done else "todo"))
        if args.force or not done:
            todo.append((monkey, day, node, exp_n))
    if args.list:
        print("%d of %d sort units to do" % (len(todo), len(units)))
        return
    if args.max_sessions:
        todo = todo[:args.max_sessions]

    for monkey, day, node, exp_n in todo:
        out = ss.out_dir_for(monkey, day, exp_n, node)
        out.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(str(out / "sort.log"))
        fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        logging.getLogger().addHandler(fh)
        try:
            res = ss.sort_session(monkey, day, exp_n, node)
            log_manifest(monkey=monkey, day=day, node=node, experiment=exp_n, status="ok",
                         n_units=res["n_units"], n_good=res["n_good"])
        except Exception as e:
            logging.getLogger().error("FAILED %s %s exp%d: %s\n%s", monkey, day, exp_n, e, traceback.format_exc())
            log_manifest(monkey=monkey, day=day, node=node, experiment=exp_n, status="failed",
                         error=str(e)[:300])
        finally:
            logging.getLogger().removeHandler(fh)
            fh.close()


if __name__ == "__main__":
    main()
