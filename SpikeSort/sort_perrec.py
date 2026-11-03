"""Sort every recording of one experiment separately, then link units across recordings.

    python sort_perrec.py m032 20260312 1                 # all recordings, then matching
    python sort_perrec.py m032 20260312 1 --match-only    # redo only the matching

Why: Kilosort's drift correction does not bridge abrupt changes between recordings, so a session sorted as one
concatenated block splits the same neuron into several time-limited units (20260312 exp1: presence_ratio ~0.5).
Output: <session>__rec<N>/ per recording (normal sort_session layout) and <session>__matched/ with the cross-
recording links (see match_units.py). The ordinary concatenated result in <session>/ is left untouched.
"""
from __future__ import annotations

import argparse
import logging
import time
import traceback

import config
import match_units
import sort_session as ss
import stage_raw

log = logging.getLogger()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("monkey")
    ap.add_argument("day")
    ap.add_argument("exp", type=int)
    ap.add_argument("--node", default=".")
    ap.add_argument("--match-only", action="store_true")
    ap.add_argument("--force", action="store_true", help="redo recordings that already have DONE.json")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    rows = ss.select_recordings(args.monkey, args.day, args.exp, args.node)
    rows = [r for r in rows if float(r["duration_s"]) >= config.MIN_RECORDING_S]
    base = ss.out_dir_for(args.monkey, args.day, args.exp, args.node)
    failed = []
    if not args.match_only:
        t0 = time.time()
        staged = stage_raw.stage_rows(ss.select_recordings(args.monkey, args.day, args.exp, args.node))
        for r in sorted(rows, key=lambda r: float(r["duration_s"])):   # shortest first: quick smoke test of the pipeline
            tag = "rec%d" % int("".join(c for c in r["recording"] if c.isdigit()))
            out = base.with_name(base.name + "__" + tag)
            if (out / "DONE.json").exists() and not args.force:
                log.info("%s already done, skipping", tag)
                continue
            try:
                ss.sort_session(args.monkey, args.day, args.exp, args.node, rec_subset=[r["recording"]],
                                tag=tag, keep_stage=True)
            except Exception as e:
                log.error("FAILED %s: %s\n%s", tag, e, traceback.format_exc())
                failed.append(tag)
        log.info("per-recording sorting finished in %.1f min, failed: %s", (time.time() - t0) / 60, failed or "none")
        if not failed:
            stage_raw.remove_staged(staged)
    match_units.match_session(base)


if __name__ == "__main__":
    main()
