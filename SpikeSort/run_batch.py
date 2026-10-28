"""Sort every sortable (monkey, day, node, experiment) found under Data_collection.

    python run_batch.py --list                    # show what would be sorted, with sizes
    python run_batch.py                           # sort everything not yet DONE (resumable)
    python run_batch.py --monkey m032 --day 20260312 --max-sessions 3
    python run_batch.py --force                   # redo already finished sessions
    python run_batch.py --until 07:30             # overnight: only start sessions expected to finish by 07:30
    touch SpikeSort/STOP                          # ask a running batch to stop after the current session

Sessions are staged (copied from the server to local disk) one ahead: while session N sorts on the GPU,
session N+1 is already being copied. With --until, a session is started only if its estimated runtime
(EST_MIN_PER_GB x size) ends before the deadline; otherwise the next one that fits is taken.

Safe to interrupt and restart: a session counts as finished only if its DONE.json exists.
Each session logs to <session out dir>/sort.log; a one-line-per-attempt history is appended
to _manifest.csv. Failures are recorded and the batch continues.
"""
from __future__ import annotations

import argparse
import csv
import logging
import threading
import time
import traceback
from datetime import datetime, timedelta
from pathlib import Path

import config
from find_recordings import exp_number, scan_all
import sort_session as ss
import stage_raw

STOP_FILE = Path(__file__).parent / "STOP"
EST_MIN_PER_GB = 1.2     # observed 0.8 (14.8 GB) to 1.1 (45 GB, incl. Phy) min per GB of raw AP data, GPU-bound
EST_FIXED_MIN = 5.0


def est_minutes(gb):
    return EST_FIXED_MIN + EST_MIN_PER_GB * gb


def parse_deadline(hhmm):
    """Next occurrence of HH:MM (today if still ahead, else tomorrow)."""
    h, m = [int(x) for x in hhmm.split(":")]
    now = datetime.now()
    d = now.replace(hour=h, minute=m, second=0, microsecond=0)
    return d if d > now else d + timedelta(days=1)


def prefetch(monkey, day, node, exp_n):
    """Stage a session in a background thread; sort_session's own staging then finds it complete."""
    def run():
        try:
            stage_raw.stage_rows(ss.select_recordings(monkey, day, exp_n, node))
            logging.getLogger().info("prefetch done: %s %s %s", monkey, day, ss.session_name(exp_n, node))
        except Exception as e:     # not fatal: sort_session stages it itself
            logging.getLogger().warning("prefetch failed (%s %s exp%d): %s", monkey, day, exp_n, e)
    t = threading.Thread(target=run, daemon=True)
    t.start()
    return t


MANIFEST_FIELDS = ["time", "monkey", "day", "node", "experiment", "status", "n_units", "n_good", "profile", "flags", "error"]


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
    ap.add_argument("--only", action="append", metavar="DAY:EXP", help="restrict to these sessions, e.g. 20260303:1 (repeatable)")
    ap.add_argument("--until", metavar="HH:MM", help="do not start a session that would end after this time")
    ap.add_argument("--no-prefetch", action="store_true", help="do not stage the next session during a sort")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    units = sort_units(scan_all(config.DATA_COLLECTION))
    todo = []
    for (monkey, day, node, exp_n), recs in units.items():
        if (args.monkey and monkey != args.monkey) or (args.day and day != args.day):
            continue
        if args.only and "%s:%d" % (day, exp_n) not in args.only:
            continue
        done = (ss.out_dir_for(monkey, day, exp_n, node) / "DONE.json").exists()
        gb = sum(int(r["ap_bytes"]) for r in recs) / 1e9
        if args.list:
            print("%-5s %s %-12s %6.1f GB %s" % (monkey, day, ss.session_name(exp_n, node), gb,
                                                 "DONE" if done else "todo"))
        if args.force or not done:
            todo.append((monkey, day, node, exp_n, gb))
    if args.list:
        print("%d of %d sort units to do" % (len(todo), len(units)))
        return
    if args.max_sessions:
        todo = todo[:args.max_sessions]

    deadline = parse_deadline(args.until) if args.until else None
    log = logging.getLogger()
    if deadline:
        log.info("batch window ends %s", deadline.strftime("%Y-%m-%d %H:%M"))
    pending = list(todo)
    prefetched = None            # (key, thread) of the session staged ahead
    while pending:
        if STOP_FILE.exists():
            log.info("STOP file found, ending batch (remove %s to run again)", STOP_FILE)
            break
        pick = None
        for k, item in enumerate(pending):      # chronological, first one that fits the window
            if deadline is None or datetime.now() + timedelta(minutes=est_minutes(item[4])) <= deadline:
                pick = k
                break
        if pick is None:
            log.info("no remaining session fits before %s, ending batch", deadline.strftime("%H:%M"))
            break
        monkey, day, node, exp_n, gb = pending.pop(pick)
        key = (monkey, day, node, exp_n)
        if prefetched:
            prefetched[1].join()                  # never copy two sessions at once
        prefetched = None
        out = ss.out_dir_for(monkey, day, exp_n, node)
        out.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(str(out / "sort.log"))
        fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        log.addHandler(fh)
        try:
            if config.STAGE:
                # stage this session first (full link speed), then copy the next one while it sorts
                stage_raw.stage_rows(ss.select_recordings(monkey, day, exp_n, node))
                if not args.no_prefetch and pending:
                    nxt = pending[0]
                    if deadline is None or datetime.now() + timedelta(minutes=est_minutes(gb) + est_minutes(nxt[4])) <= deadline:
                        prefetched = ((nxt[0], nxt[1], nxt[2], nxt[3]), prefetch(*nxt[:4]))
            res = ss.sort_session(monkey, day, exp_n, node)
            log_manifest(monkey=monkey, day=day, node=node, experiment=exp_n, status="ok",
                         n_units=res["n_units"], n_good=res["n_good"], profile=res["profile"], flags=res["flags"])
        except Exception as e:
            log.error("FAILED %s %s exp%d: %s\n%s", monkey, day, exp_n, e, traceback.format_exc())
            log_manifest(monkey=monkey, day=day, node=node, experiment=exp_n, status="failed",
                         error=str(e)[:300])
        finally:
            log.removeHandler(fh)
            fh.close()
    if prefetched:
        prefetched[1].join()
        log.info("note: %s %s %s exp%d is staged in %s but not sorted (next run reuses it)", *prefetched[0], config.STAGE_ROOT)


if __name__ == "__main__":
    main()
