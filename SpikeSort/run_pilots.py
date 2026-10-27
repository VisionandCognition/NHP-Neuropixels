"""The queued geometry pilots (rationale and results: PILOT_LOG.md). Sequential, GPU-bound.

    python run_pilots.py            # run the whole queue; failures are logged, the queue continues
"""
import logging

import config
import pilot_slice as ps

QUEUE = [
    # (monkey, day, exp, variants, reuse_pp)
    # T1 single column: does dminx=103 (lateral contact distance) change yield vs the 32 um default?
    ("m032", "20260323", 1, ["default", "dminx103"], config.SCRATCH / "smoke" / "pp"),
    # T2 four_bank_even: regular 40 um line over 15 mm
    ("m032", "20260819", 1, ["default", "dminx103"], None),
    # T3 custom patches: strays dropped? sort patches separately?
    ("m032", "20260317", 1, ["default", "dminx103", "drop_strays", "drop_strays_dminx103", "by_patch"], None),
    # T4 very fragmented (38-51 segments, mixed pitch)
    ("m032", "20260325", 1, ["default", "drop_strays_dminx103", "by_patch"], None),
]

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    for monkey, day, exp, variants, reuse in QUEUE:
        try:
            ps.run_pilot(monkey, day, exp, variants, reuse_pp=reuse)
        except Exception:
            logging.getLogger("pilot").exception("pilot %s %s exp%d FAILED", monkey, day, exp)
