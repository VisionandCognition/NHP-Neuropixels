"""Link units across separately sorted recordings of one experiment (see sort_perrec.py).

For each recording k (in order) every unit is compared with the most recent template of each linked unit seen
in the previous MAX_GAP recordings. A pair is a candidate when the units sit within MAX_DIST_UM of each other,
their peak amplitudes agree within AMP_RATIO and the templates' cosine similarity is >= MIN_COS; candidates are
assigned one-to-one (Hungarian). Unmatched units start a new linked unit. Only units with >= MIN_SPIKES spikes
take part; the rest keep a link id of their own and are flagged `matched_eligible = 0`.

Output in <session>__matched/: unit_matches.csv (one row per recording-unit), linked_units.csv (one row per
linked unit), spikes_matched.npz (spike times per linked unit, seconds on the OpenEphys timestamp clock).
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

log = logging.getLogger(__name__)

MIN_COS = 0.85
MAX_DIST_UM = 40.0
AMP_RATIO = (0.6, 1 / 0.6)
MAX_GAP = 3
MIN_SPIKES = 50


def _rec_dirs(base):
    base = Path(base)
    pat = re.compile(re.escape(base.name) + r"__rec(\d+)$")
    found = []
    for d in base.parent.iterdir():
        m = pat.fullmatch(d.name)
        if m and (d / "DONE.json").exists() and (d / "units.csv").exists():
            found.append((int(m.group(1)), d))
    return sorted(found)


def _load(d):
    import spikeinterface.full as si
    units = pd.read_csv(str(d / "units.csv"), index_col=0)
    an = si.load_sorting_analyzer(str(d / "analyzer"))
    tmpl = np.asarray(an.get_extension("templates").get_data(operator="average"), dtype=np.float32)
    chans = [str(c) for c in an.channel_ids]
    ids = [int(u) for u in an.unit_ids]
    return units, tmpl, chans, ids


def match_session(base, min_cos=MIN_COS):
    base = Path(base)
    recs = _rec_dirs(base)
    if len(recs) < 2:
        raise RuntimeError("need >= 2 finished per-recording sorts next to %s, found %d" % (base, len(recs)))
    out = base.with_name(base.name + "__matched")
    out.mkdir(parents=True, exist_ok=True)

    data = {n: _load(d) for n, d in recs}
    all_ch = sorted({c for _, (_, _, chans, _) in data.items() for c in chans},
                    key=lambda c: int(re.sub(r"\D", "", c) or 0))
    cidx = {c: i for i, c in enumerate(all_ch)}
    nsamp = {t.shape[1] for _, t, _, _ in data.values()}
    if len(nsamp) != 1:
        raise RuntimeError("template lengths differ across recordings: %s" % nsamp)

    def full(tmpl, chans):
        a = np.zeros((tmpl.shape[0], tmpl.shape[1], len(all_ch)), dtype=np.float32)
        a[:, :, [cidx[c] for c in chans]] = tmpl
        return a.reshape(tmpl.shape[0], -1)

    linked = []          # dict(id, tmpl (flat), x, y, amp, last_k)
    rows = []
    for k, (n, _) in enumerate(recs):
        units, tmpl, chans, ids = data[n]
        flat = full(tmpl, chans)
        amp = -tmpl.min(axis=(1, 2))
        uu = units.reindex(ids)
        elig = (uu["num_spikes"].values >= MIN_SPIKES)
        live = [g for g in linked if k - g["last_k"] <= MAX_GAP]
        cost = np.full((len(ids), len(live)), 10.0)
        sims = np.zeros_like(cost)
        for i in range(len(ids)):
            if not elig[i]:
                continue
            for j, g in enumerate(live):
                dist = np.hypot(uu["x_um"].values[i] - g["x"], uu["y_um"].values[i] - g["y"])
                r = amp[i] / max(g["amp"], 1e-6)
                if dist > MAX_DIST_UM or not (AMP_RATIO[0] <= r <= AMP_RATIO[1]):
                    continue
                c = float(flat[i] @ g["tmpl"] / (np.linalg.norm(flat[i]) * np.linalg.norm(g["tmpl"]) + 1e-12))
                sims[i, j] = c
                if c >= min_cos:
                    cost[i, j] = 1 - c
        assign = {}
        if len(live) and len(ids):
            ri, ci = linear_sum_assignment(cost)
            assign = {i: j for i, j in zip(ri, ci) if cost[i, j] < 10.0}
        for i, uid in enumerate(ids):
            if i in assign:
                g = live[assign[i]]
                sim = sims[i, assign[i]]
            else:
                g = dict(id=len(linked), n_rec=0)
                linked.append(g)
                sim = np.nan
            g.update(tmpl=flat[i], x=float(uu["x_um"].values[i]), y=float(uu["y_um"].values[i]),
                     amp=float(amp[i]), last_k=k)
            g["n_rec"] = g.get("n_rec", 0) + 1
            rows.append(dict(linked_id=g["id"], recording_num=n, unit_id=uid, sim_to_prev=sim,
                             matched_eligible=int(elig[i]), num_spikes=int(uu["num_spikes"].values[i]),
                             snr=float(uu["snr"].values[i]), presence_ratio=float(uu["presence_ratio"].values[i]),
                             peak_channel_id=str(uu["peak_channel_id"].values[i]), x_um=g["x"], y_um=g["y"],
                             tier=str(uu["tier"].values[i]) if "tier" in uu else ""))
    m = pd.DataFrame(rows)
    m.to_csv(str(out / "unit_matches.csv"), index=False)

    agg = m.groupby("linked_id").agg(n_recordings=("recording_num", "nunique"), first_rec=("recording_num", "min"),
                                     last_rec=("recording_num", "max"), n_spikes=("num_spikes", "sum"),
                                     best_snr=("snr", "max"), any_usable=("tier", lambda t: int(t.isin(["good", "usable"]).any())),
                                     eligible=("matched_eligible", "max"))
    agg["recordings"] = m.groupby("linked_id")["recording_num"].apply(lambda s: ",".join(str(x) for x in sorted(s)))
    agg.to_csv(str(out / "linked_units.csv"))

    # merged spike trains
    lid, tt, rr = [], [], []
    for n, d in recs:
        z = np.load(str(d / "spikes.npz"))
        sub = m[m["recording_num"] == n].set_index("unit_id")["linked_id"]
        sel = np.isin(z["spike_unit"], sub.index.values)
        lid.append(sub.reindex(z["spike_unit"][sel]).values.astype(np.int64))
        tt.append(z["spike_time_s"][sel])
        rr.append(np.full(sel.sum(), n, dtype=np.int16))
    lid, tt, rr = np.concatenate(lid), np.concatenate(tt), np.concatenate(rr)
    o = np.lexsort((tt, lid))
    np.savez(str(out / "spikes_matched.npz"), linked_ids=agg.index.values.astype(np.int64), spike_linked=lid[o],
             spike_time_s=tt[o], spike_recording=rr[o], params=json.dumps(
                 dict(min_cos=min_cos, max_dist_um=MAX_DIST_UM, amp_ratio=AMP_RATIO, max_gap=MAX_GAP, min_spikes=MIN_SPIKES)))

    n_rec = len(recs)
    el = agg[agg["eligible"] == 1]
    log.info("matched %d recordings: %d unit rows -> %d linked units (%d eligible); eligible linked units present in "
             ">=50%% of recordings: %d (concatenated sort had 13%% of units in all recordings)",
             n_rec, len(m), len(agg), len(el), int((el["n_recordings"] >= n_rec / 2).sum()))
    pairs = m["sim_to_prev"].dropna()
    log.info("accepted-link cosine similarity: n=%d, median %.3f, 10th pct %.3f", len(pairs), pairs.median(),
             pairs.quantile(0.1) if len(pairs) else float("nan"))
    return out
