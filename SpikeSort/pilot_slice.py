"""Short-slice pilots: compare Kilosort variants on a few minutes of one session.

Why: electrode configurations differ a lot between sessions (see geometry.py and PILOT_LOG.md),
and Kilosort's defaults (e.g. dminx=32 um) are tuned for Neuropixels 1/2. A 3-minute slice is
cheap (~4 GB, ~5 min per variant) and shows how parameters / handling choices change the unit
yield before a multi-hour sort is spent on a whole day.

The slice is preprocessed once (cached in config.SCRATCH/pilots/<key>/pp) and every variant is
run on that same data. One row per variant is appended to pilots/results.csv.

Variants (see VARIANTS): Kilosort parameter overrides, dropping isolated stray channels
(patches with fewer than N channels), and sorting each contiguous patch separately.
"""
from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd

import config
import geometry
import sort_session as ss

log = logging.getLogger("pilot")
PILOT_ROOT = config.SCRATCH / "pilots"
RESULTS = Path(__file__).parent / "pilots" / "results.csv"

# Base: nblocks=1 because a 3-min slice cannot show session-scale drift; this isolates the
# effect of the geometry-related choices. Production keeps nblocks=5 (config.KS4_PARAMS).
BASE_KS = dict(config.KS4_PARAMS, nblocks=1)

VARIANTS = {
    "default": dict(),
    "dminx103": dict(ks=dict(dminx=103)),
    "drop_strays": dict(drop_strays=8),
    "drop_strays_dminx103": dict(drop_strays=8, ks=dict(dminx=103)),
    "by_patch": dict(drop_strays=8, by_patch=True, ks=dict(dminx=103)),
}


def prepare_slice(monkey, day, exp, dur_s=180.0, start_s=120.0, reuse_pp=None):
    """Return (si, preprocessed binary recording, info) for a slice of the longest recording."""
    import spikeinterface.full as si
    key = "%s_%s_exp%d_%ds" % (monkey, day, exp, dur_s)
    d = PILOT_ROOT / key
    if reuse_pp is not None and (Path(reuse_pp) / "binary.json").exists():
        return key, si, si.load(str(reuse_pp)), dict(reused=str(reuse_pp))
    if (d / "pp" / "binary.json").exists():
        return key, si, si.load(str(d / "pp")), dict(cached=True)
    rows = ss.select_recordings(monkey, day, exp)
    _, rec = ss.load_recording(rows, exp)
    si.set_global_job_kwargs(n_jobs=config.N_JOBS, chunk_duration="1s", progress_bar=False)
    fs = rec.get_sampling_frequency()
    lens = [rec.get_num_samples(segment_index=i) for i in range(rec.get_num_segments())]
    i = int(np.argmax(lens))
    a = int(min(start_s * fs, max(0, lens[i] - dur_s * fs)))
    sl = rec.select_segments([i]).frame_slice(a, int(a + dur_s * fs))
    pp, info = ss.preprocess(si, sl)
    d.mkdir(parents=True, exist_ok=True)
    saved = pp.save(folder=str(d / "pp"), format="binary", overwrite=True)
    info.update(recording=rows[i]["recording"], start_s=a / fs)
    return key, si, saved, info


def build_variant_recording(si, rec, spec):
    """Apply channel dropping; return (recording, list of (label, channel_ids) groups)."""
    loc = rec.get_channel_locations()
    labels = geometry.patch_labels(loc)
    ids = np.asarray(rec.channel_ids)
    n_drop = spec.get("drop_strays", 0)
    keep = np.ones(len(ids), dtype=bool)
    if n_drop:
        counts = np.bincount(labels)
        keep = counts[labels] >= n_drop
        rec = rec.select_channels(ids[keep])
        labels = labels[keep]
        ids = ids[keep]
    groups = [(int(g), ids[labels == g]) for g in np.unique(labels)]
    return rec, groups, dict(n_patches=len(groups), n_channels=int(len(ids)),
                             patch_sizes=sorted([len(g[1]) for g in groups], reverse=True)[:12])


def run_ks(si, rec, folder, ks):
    return si.run_sorter("kilosort4", rec, folder=str(folder), remove_existing_folder=True,
                         verbose=False, **ks)


def evaluate(si, sorting, rec, secs):
    sorting = sorting.remove_empty_units()
    n_units = sorting.get_num_units()
    if n_units == 0:
        return dict(n_units=0, n_good=0, n_spikes=0, median_fr_hz=0.0), pd.DataFrame()
    an = si.create_sorting_analyzer(sorting, rec, sparse=True, format="memory")
    an.compute("random_spikes", method="uniform", max_spikes_per_unit=300)
    an.compute("waveforms", ms_before=1.0, ms_after=2.0)
    an.compute("templates", operators=["average"])
    an.compute("noise_levels")
    an.compute("spike_amplitudes")
    an.compute("quality_metrics", metric_names=["num_spikes", "firing_rate", "presence_ratio", "snr",
                                                "isi_violation", "amplitude_cutoff"])
    q = an.get_extension("quality_metrics").get_data()
    c = config.CURATION
    good = ((q["amplitude_cutoff"] <= c["amplitude_cutoff_max"]) & (q["isi_violations_ratio"] <= c["isi_violations_ratio_max"])
            & (q["presence_ratio"] >= c["presence_ratio_min"]) & (q["snr"] >= c["snr_min"]))
    q["good"] = good.fillna(False).astype(int)
    summ = dict(n_units=int(len(q)), n_good=int(q["good"].sum()), n_spikes=int(q["num_spikes"].sum()),
                median_fr_hz=float(np.median(q["firing_rate"])), median_snr=float(np.nanmedian(q["snr"])),
                frac_isi_ok=float((q["isi_violations_ratio"] <= c["isi_violations_ratio_max"]).mean()))
    return summ, q


def ks_settings_used(folder):
    out = {}
    for p in Path(folder).rglob("ops.npy"):
        try:
            ops = np.load(str(p), allow_pickle=True).item()
            for k in ("dmin", "dminx", "nearest_chans", "nblocks", "n_chan_bin"):
                if k in ops:
                    out[k] = float(ops[k]) if ops[k] is not None else None
            if "settings" in ops:
                for k in ("dmin", "dminx", "nearest_chans", "nblocks"):
                    out.setdefault(k, ops["settings"].get(k))
        except Exception:
            pass
        break
    return out


def run_variant(si, rec_full, key, name, secs, check_only=False):
    spec = VARIANTS[name]
    ks = dict(BASE_KS, **spec.get("ks", {}))
    rec, groups, gsum = build_variant_recording(si, rec_full, spec)
    log.info("[%s/%s] channels %d -> %d, patches %d, ks=%s", key, name, rec_full.get_num_channels(),
             gsum["n_channels"], gsum["n_patches"], ks)
    if check_only:
        return None
    folder = PILOT_ROOT / key / ("ks_" + name)
    t = time.time()
    failed = []
    if spec.get("by_patch"):
        sortings, failed = [], []
        for g, gids in groups:
            if len(gids) < 20:          # too small for template matching
                continue
            sub = rec.select_channels(gids)
            kg = dict(ks, nblocks=1 if len(gids) < 150 else ks["nblocks"])
            try:
                sortings.append(run_ks(si, sub, folder / ("patch%d" % g), kg))
            except Exception as e:      # e.g. a quiet patch: too few spikes to initialise templates
                log.warning("patch %d (%d ch) failed, skipped: %s", g, len(gids), str(e).strip().splitlines()[-1][:120])
                failed.append("patch%d(%dch)" % (g, len(gids)))
        sorting = si.aggregate_units(sortings) if len(sortings) > 1 else sortings[0]
        used = ks_settings_used(folder)
    else:
        sorting = run_ks(si, rec, folder, ks)
        used = ks_settings_used(folder)
    ks_s = time.time() - t
    summ, q = evaluate(si, sorting, rec, secs)
    if len(q):
        q.to_csv(str(PILOT_ROOT / key / ("units_" + name + ".csv")))
    row = dict(time=pd.Timestamp.now().isoformat(timespec="seconds"), key=key, variant=name,
               ks_runtime_s=round(ks_s, 1), n_channels=gsum["n_channels"], n_patches=gsum["n_patches"],
               ks_params=json.dumps(ks, default=str), ks_used=json.dumps(used),
               note=("failed: " + ",".join(failed)) if failed else "", **summ)
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(str(RESULTS), mode="a", header=not RESULTS.exists(), index=False)
    log.info("[%s/%s] %s", key, name, {k: row[k] for k in ("n_units", "n_good", "n_spikes", "ks_runtime_s")})
    return row


def run_pilot(monkey, day, exp, variants, dur_s=180.0, reuse_pp=None, check_only=False):
    key, si, rec, info = prepare_slice(monkey, day, exp, dur_s, reuse_pp=reuse_pp)
    log.info("slice %s ready: %d channels, %.0f s, %s", key, rec.get_num_channels(), rec.get_total_duration(), info)
    for v in variants:
        try:
            run_variant(si, rec, key, v, rec.get_total_duration(), check_only)
        except Exception:
            log.exception("variant %s/%s FAILED", key, v)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("monkey")
    ap.add_argument("day")
    ap.add_argument("experiment", type=int)
    ap.add_argument("variants", nargs="+", choices=sorted(VARIANTS))
    ap.add_argument("--dur", type=float, default=180.0)
    ap.add_argument("--check-only", action="store_true", help="build variants, no Kilosort")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_pilot(args.monkey, args.day, args.experiment, args.variants, args.dur, check_only=args.check_only)


if __name__ == "__main__":
    main()
