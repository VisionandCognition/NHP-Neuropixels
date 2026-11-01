"""Unit tiers written to units.csv: `good` (strict, whole-session) > `usable` > other; `artifact` =
peak on a PREPROC.always_bad channel (overrides `good`, which is left untouched as the raw flag).

`good` is computed in sort_session.unit_table. `usable` / `su_like` / `tier` are derived from the
metric columns alone, so add_tiers() can also be applied to finished sessions (backfill_tiers.py).
"""
import numpy as np

import config


def add_tiers(q):
    u = config.USABLE
    amp_ok = ~(q["amplitude_cutoff"] >= u["amplitude_cutoff_max"])   # NaN passes
    art = q["peak_channel_id"].astype(str).isin(config.PREPROC.get("always_bad", ()))
    ok = (amp_ok
          & (q["isi_violations_ratio"] < u["isi_violations_ratio_max"])
          & (q["snr"] > u["snr_min"])
          & (q["num_spikes"] >= u["min_spikes"])
          & ~art)
    q["usable"] = ok.astype(int)
    q["su_like"] = (ok & (q["isi_violations_ratio"] < config.SU_ISI_MAX)).astype(int)
    q["tier"] = np.where(art, "artifact", np.where(q["good"] == 1, "good", np.where(ok, "usable", "other")))
    return q
