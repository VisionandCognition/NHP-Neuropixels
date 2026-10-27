"""Sort one (monkey, day, node, experiment): all its recordings concatenated.

    python sort_session.py m032 20260312 1                 # flat layout, experiment1
    python sort_session.py m032 20260224 1 --node "m032-2026-02-24_12-59-28/Record Node 110"
    python sort_session.py m032 20260312 1 --ks nblocks=1 --dry-run

Stages: select recordings -> load (SpikeInterface/OpenEphys) -> preprocess -> save binary to
scratch -> Kilosort4 -> SortingAnalyzer + quality metrics -> curation flags -> export.
See README.md for the output layout. Written against SpikeInterface 0.104 / Kilosort 4.1;
NOT yet run against real data -- see README "Validation status".
"""
from __future__ import annotations

import argparse
import json
import logging
import platform
import shutil
import time
from datetime import datetime
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

import config
import export_spikes
from find_recordings import exp_number, rec_number, scan_day

log = logging.getLogger("sort_session")


# --------------------------------------------------------------------------------------
# paths / selection
# --------------------------------------------------------------------------------------
def node_slug(node):
    return "" if node in (".", "") else node.replace("\\", "/").replace("/", "__").replace(" ", "_")


def session_name(exp_n, node="."):
    slug = node_slug(node)
    return ("%s__exp%d" % (slug, exp_n)) if slug else ("exp%d" % exp_n)


def out_dir_for(monkey, day, exp_n, node="."):
    return config.OUTPUT_ROOT / monkey / day / session_name(exp_n, node)


def select_recordings(monkey, day, exp_n, node="."):
    """Rows (sorted by recording number) of all recordings in this experiment/node."""
    day_dir = config.DATA_COLLECTION / monkey / day
    rows = [r for r in scan_day(day_dir, monkey)
            if exp_number(r) == exp_n and r["node"] == node]
    if not rows:
        raise RuntimeError("no recordings for %s %s experiment%d node=%r" % (monkey, day, exp_n, node))
    rows = sorted(rows, key=rec_number)
    missing = [r["recording"] for r in rows if not int(r["has_ap_dat"] or 0)]
    if missing:
        raise RuntimeError("recordings without AP continuous.dat: %s" % missing)
    if not rows[0]["settings_xml"]:
        raise RuntimeError("no settings.xml next to experiment%d (geometry unknown)" % exp_n)
    if len({r["ap_num_channels"] for r in rows}) > 1:
        raise RuntimeError("recordings differ in AP channel count: %s" % [r["ap_num_channels"] for r in rows])
    return rows


# --------------------------------------------------------------------------------------
# timing info straight from the OpenEphys files (same source the MUA pipeline aligns with)
# --------------------------------------------------------------------------------------
def _contiguous(sample_numbers, chunk=10_000_000):
    for s in range(0, len(sample_numbers), chunk):
        a = np.asarray(sample_numbers[s:s + chunk + 1])
        if len(a) > 1 and np.any(np.diff(a) != 1):
            return False
    return True


def segment_info(row):
    d = Path(row["path"]) / "continuous" / row["ap_folder"]
    n_samples = int(round(int(row["ap_bytes"]) / 2 / int(row["ap_num_channels"])))
    info = dict(recording=row["recording"], n_samples=n_samples,
                first_timestamp_s=float("nan"), first_sample_number=-1, contiguous=True)
    try:
        ts = np.load(str(d / "timestamps.npy"), mmap_mode="r")
        sn = np.load(str(d / "sample_numbers.npy"), mmap_mode="r")
        info["first_timestamp_s"] = float(ts[0])
        info["first_sample_number"] = int(sn[0])
        info["contiguous"] = bool(_contiguous(sn))
    except Exception as e:
        log.warning("%s: no usable timestamps/sample_numbers (%s)", row["recording"], e)
    if not info["contiguous"]:
        log.warning("%s: sample_numbers have gaps; spike times assume contiguous samples!", row["recording"])
    return info


# --------------------------------------------------------------------------------------
# load + preprocess
# --------------------------------------------------------------------------------------
def load_recording(rows, exp_n):
    import spikeinterface.full as si
    root = Path(rows[0]["settings_xml"]).parent       # Record-Node folder (holds experimentN/)
    stream_names, _ = si.get_neo_streams("openephysbinary", str(root))
    ap_streams = [s for s in stream_names if config.AP_STREAM_KEY in s]
    if len(ap_streams) != 1:
        raise RuntimeError("expected one %s stream in %s, got %s" % (config.AP_STREAM_KEY, root, stream_names))
    # neo exposes each experiment as a block, in sorted order of the experiment folders present
    exps = sorted(int(p.name.replace("experiment", "")) for p in root.glob("experiment*") if p.is_dir())
    rec = si.read_openephys(str(root), stream_name=ap_streams[0], block_index=exps.index(exp_n))

    # segments must correspond 1:1, in order, to our recording rows
    if rec.get_num_segments() != len(rows):
        raise RuntimeError("SpikeInterface sees %d segments, found %d recording folders"
                           % (rec.get_num_segments(), len(rows)))
    for i, r in enumerate(rows):
        expected = int(round(int(r["ap_bytes"]) / 2 / int(r["ap_num_channels"])))
        if rec.get_num_samples(segment_index=i) != expected:
            raise RuntimeError("segment %d (%s) has %d samples, file implies %d -- segment order "
                               "differs from recording-number order" %
                               (i, r["recording"], rec.get_num_samples(segment_index=i), expected))
    if not rec.has_probe():
        raise RuntimeError("no probe geometry attached by the OpenEphys reader; attach from "
                           "settings.xml via probeinterface before sorting (pilot TODO)")
    return si, rec


def preprocess(si, rec):
    """phase_shift -> highpass -> drop bad channels -> reference. Returns (rec_pp, info)."""
    pp = config.PREPROC
    # NP1 ADC multiplexing shifts. SpikeInterface sets this property for some readers;
    # if the OpenEphys reader did not, compute the standard NP1 (12 ch/ADC, 13 cycles) table.
    if "inter_sample_shift" not in rec.get_property_keys():
        from probeinterface.neuropixels_tools import get_neuropixels_sample_shifts
        shifts = get_neuropixels_sample_shifts(rec.get_num_channels(), 12)
        rec.set_property("inter_sample_shift", shifts)
        log.info("set inter_sample_shift from NP1 ADC table")
    out = si.phase_shift(rec)
    out = si.highpass_filter(out, freq_min=pp["highpass_hz"])
    bad_ids, labels = si.detect_bad_channels(out, method=pp["bad_channel_method"])
    log.info("bad channels removed: %d (%s)", len(bad_ids),
             {k: int((labels == k).sum()) for k in np.unique(labels)})
    if len(bad_ids):
        out = out.remove_channels(bad_ids)
    if pp["reference"] == "destripe":
        out = si.highpass_spatial_filter(out)
    else:
        out = si.common_reference(out, operator="median", reference="global")
    out = si.astype(out, "int16")   # keep scratch the same size as raw
    return out, dict(bad_channel_ids=[str(b) for b in bad_ids],
                     n_channels_kept=int(out.get_num_channels()), preproc=pp)


# --------------------------------------------------------------------------------------
# sorting
# --------------------------------------------------------------------------------------
def run_kilosort(si, rec_saved, sorter_dir, ks_params):
    attempts = [dict(ks_params)] + [dict(ks_params, **fb) for fb in config.KS4_FALLBACKS]
    for i, params in enumerate(attempts):
        try:
            log.info("Kilosort4 attempt %d params=%s", i + 1, params)
            sorting = si.run_sorter("kilosort4", rec_saved, folder=str(sorter_dir),
                                    remove_existing_folder=True, verbose=True, **params)
            return sorting, params
        except Exception as e:
            log.error("Kilosort4 attempt %d failed: %s", i + 1, e)
            if i == len(attempts) - 1:
                raise
            try:
                import torch
                torch.cuda.empty_cache()
            except Exception:
                pass


# --------------------------------------------------------------------------------------
# analyzer, metrics, curation
# --------------------------------------------------------------------------------------
def build_analyzer(si, sorting, rec_saved, analyzer_dir):
    if analyzer_dir.exists():
        shutil.rmtree(str(analyzer_dir))
    an = si.create_sorting_analyzer(sorting, rec_saved, sparse=True,
                                    format="binary_folder", folder=str(analyzer_dir))
    an.compute("random_spikes", method="uniform", max_spikes_per_unit=500)
    an.compute("waveforms", ms_before=1.0, ms_after=2.0)
    an.compute("templates", operators=["average", "std"])
    an.compute("noise_levels")
    an.compute("spike_amplitudes")
    an.compute("correlograms", window_ms=50.0, bin_ms=1.0)
    an.compute("unit_locations", method="monopolar_triangulation")
    an.compute("template_similarity")
    metric_names = ["num_spikes", "firing_rate", "presence_ratio", "snr",
                    "isi_violation", "amplitude_cutoff"]
    if config.UNITREFINE:
        an.compute("spike_locations")
        an.compute("template_metrics", include_multi_channel_metrics=True)
        metric_names += ["rp_violation", "sliding_rp_violation", "amplitude_median",
                         "amplitude_cv", "firing_range", "drift", "sd_ratio"]
    an.compute("quality_metrics", metric_names=metric_names)
    return an


def unit_table(si, an, rec_saved, sorter_out):
    q = an.get_extension("quality_metrics").get_data().copy()
    q.index.name = "unit_id"
    locs = an.get_extension("unit_locations").get_data()
    q["x_um"] = locs[:, 0]
    q["y_um"] = locs[:, 1]          # along the probe, same convention as NP_PROBE.ELECTRODE_YPOS
    ext = si.get_template_extremum_channel(an, peak_sign="neg", outputs="id")
    ch_loc = dict(zip(rec_saved.channel_ids, rec_saved.get_channel_locations()))
    q["peak_channel_id"] = [str(ext[u]) for u in q.index]
    q["peak_channel_y_um"] = [float(ch_loc[ext[u]][1]) for u in q.index]

    c = config.CURATION
    q["good"] = ((q["amplitude_cutoff"] < c["amplitude_cutoff_max"])
                 & (q["isi_violations_ratio"] < c["isi_violations_ratio_max"])
                 & (q["presence_ratio"] > c["presence_ratio_min"])
                 & (q["snr"] > c["snr_min"])).astype(int)

    q["ks_label"] = ""
    tsv = Path(sorter_out) / "cluster_KSLabel.tsv"
    if tsv.is_file():
        ks = pd.read_csv(str(tsv), sep="\t").set_index("cluster_id")["KSLabel"]
        q["ks_label"] = [ks.get(int(u), "") for u in q.index]

    q["unitrefine_noise"] = ""
    q["unitrefine_sua"] = ""
    if config.UNITREFINE:
        try:
            import spikeinterface.curation as sc
            nl = sc.model_based_label_units(sorting_analyzer=an, trust_model=True,
                                            repo_id="SpikeInterface/UnitRefine_noise_neural_classifier")
            q["unitrefine_noise"] = nl["prediction"].reindex(q.index).values
            neural = an.select_units(nl.index[nl["prediction"] != "noise"])
            sl = sc.model_based_label_units(sorting_analyzer=neural, trust_model=True,
                                            repo_id="SpikeInterface/UnitRefine_sua_mua_classifier")
            q["unitrefine_sua"] = sl["prediction"].reindex(q.index).values
        except Exception as e:  # advisory only
            log.error("UnitRefine failed: %s", e)
    return q


def qc_figures(units, spike_unit, spike_time_s, fig_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig_dir.mkdir(parents=True, exist_ok=True)
    y = units["peak_channel_y_um"].to_dict()
    good = units["good"].to_dict()

    fig, ax = plt.subplots(figsize=(11, 5))
    sub = np.random.default_rng(0).choice(len(spike_unit), size=min(len(spike_unit), 300_000), replace=False) \
        if len(spike_unit) else np.array([], dtype=int)
    yy = np.array([y.get(u, np.nan) for u in spike_unit[sub]])
    gg = np.array([good.get(u, 0) for u in spike_unit[sub]], dtype=bool)
    ax.scatter(spike_time_s[sub][~gg], yy[~gg], s=0.3, c="0.7", linewidths=0, label="non-good")
    ax.scatter(spike_time_s[sub][gg], yy[gg], s=0.3, c="C3", linewidths=0, label="good")
    ax.set_xlabel("time on OE timeline (s)")
    ax.set_ylabel("peak channel y (um from tip)")
    ax.set_title("spike raster (drift check): %d units, %d good" % (len(units), int(units["good"].sum())))
    ax.legend(markerscale=20)
    fig.savefig(str(fig_dir / "raster.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)

    fig, axs = plt.subplots(1, 4, figsize=(15, 3.5))
    for ax, col in zip(axs, ["firing_rate", "snr", "isi_violations_ratio", "amplitude_cutoff"]):
        v = units[col].replace([np.inf, -np.inf], np.nan).dropna()
        ax.hist(v, bins=40)
        ax.set_title(col)
    fig.savefig(str(fig_dir / "metrics.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------------------
def package_versions():
    v = {}
    for p in ("spikeinterface", "kilosort", "probeinterface", "neo", "torch", "numpy"):
        try:
            v[p] = metadata.version(p)
        except Exception:
            v[p] = None
    return v


def sort_session(monkey, day, exp_n, node=".", ks_overrides=None, keep_scratch=False, dry_run=False):
    t0 = time.time()
    out = out_dir_for(monkey, day, exp_n, node)
    rows = select_recordings(monkey, day, exp_n, node)
    keep_idx = [i for i, r in enumerate(rows) if float(r["duration_s"]) >= config.MIN_RECORDING_S]
    if not keep_idx:
        raise RuntimeError("all recordings shorter than %.0f s" % config.MIN_RECORDING_S)
    dropped = [rows[i]["recording"] for i in range(len(rows)) if i not in keep_idx]
    total_gb = sum(int(rows[i]["ap_bytes"]) for i in keep_idx) / 1e9
    log.info("%s %s %s: %d recordings kept (dropped %s), %.1f GB, %.2f h", monkey, day,
             session_name(exp_n, node), len(keep_idx), dropped or "none", total_gb,
             sum(float(rows[i]["duration_s"]) for i in keep_idx) / 3600)
    if dry_run:
        return None

    scratch = config.SCRATCH / ("%s_%s_%s" % (monkey, day, session_name(exp_n, node)))
    scratch.mkdir(parents=True, exist_ok=True)
    free_gb = shutil.disk_usage(str(scratch)).free / 1e9
    if free_gb < 1.3 * total_gb:
        raise RuntimeError("scratch has %.0f GB free, need ~%.0f GB (%s)" % (free_gb, 1.3 * total_gb, scratch))
    out.mkdir(parents=True, exist_ok=True)
    if (out / "DONE.json").exists():
        (out / "DONE.json").unlink()   # mark as incomplete while re-running

    segs = [segment_info(rows[i]) for i in keep_idx]
    si, rec = load_recording(rows, exp_n)
    si.set_global_job_kwargs(n_jobs=config.N_JOBS, chunk_duration="1s", progress_bar=True)
    fs = rec.get_sampling_frequency()

    # drop short recordings, then join the rest into one continuous segment
    parts = [rec.select_segments([i]) for i in keep_idx]
    rec = parts[0] if len(parts) == 1 else si.concatenate_recordings(parts)
    seg_starts = np.concatenate([[0], np.cumsum([s["n_samples"] for s in segs])[:-1]])

    rec_pp, pp_info = preprocess(si, rec)
    log.info("saving preprocessed binary to %s", scratch / "preprocessed")
    rec_saved = rec_pp.save(folder=str(scratch / "preprocessed"), format="binary", overwrite=True)

    ks_params = dict(config.KS4_PARAMS, **(ks_overrides or {}))
    sorting, ks_used = run_kilosort(si, rec_saved, scratch / "ks4", ks_params)
    sorting = sorting.remove_empty_units()
    log.info("Kilosort4 found %d units", sorting.get_num_units())

    # --- outputs --------------------------------------------------------------------
    ks_out = scratch / "ks4" / "sorter_output"
    ks_dst = out / "ks4_raw"
    if ks_dst.exists():
        shutil.rmtree(str(ks_dst))
    if ks_out.exists():
        shutil.copytree(str(ks_out), str(ks_dst), ignore=shutil.ignore_patterns("*.dat", "*.bin"))

    if sorting.get_num_units() == 0:
        units = pd.DataFrame(columns=["good"])
        n_units = n_good = 0
    else:
        an = build_analyzer(si, sorting, rec_saved, out / "analyzer")
        units = unit_table(si, an, rec_saved, ks_out)
        n_units, n_good = len(units), int(units["good"].sum())
        units.to_csv(str(out / "units.csv"))

        sv = sorting.to_spike_vector()
        spike_unit = np.asarray(sorting.unit_ids)[sv["unit_index"]].astype(np.int64)
        spike_sample = sv["sample_index"].astype(np.int64)
        first_ts = [s["first_timestamp_s"] for s in segs]
        spike_seg, spike_t = export_spikes.spikes_to_times(spike_sample, seg_starts, first_ts, fs)
        unit_ids = [int(u) for u in units.index]
        export_spikes.write_npz(str(out / "spikes.npz"), unit_ids, spike_unit, spike_sample,
                                spike_seg, spike_t, segs, fs)
        export_spikes.write_mat(str(out / "spikes.mat"), units, unit_ids, spike_unit, spike_t, segs, fs)
        qc_figures(units, spike_unit, spike_t, out / "figures")

    provenance = dict(
        monkey=monkey, day=day, node=node, experiment=exp_n,
        finished=datetime.now().isoformat(timespec="seconds"), runtime_min=round((time.time() - t0) / 60, 1),
        host=platform.node(), versions=package_versions(),
        recordings=[rows[i]["path"] for i in keep_idx], dropped_short=dropped,
        segments=segs, sampling_rate=fs, preprocessing=pp_info,
        kilosort_params=ks_used, curation=config.CURATION, unitrefine=config.UNITREFINE,
        channel_locations_um=rec_saved.get_channel_locations().tolist(),
        channel_ids=[str(c) for c in rec_saved.channel_ids],
        n_units=n_units, n_good=n_good,
    )
    (out / "sorting_params.json").write_text(json.dumps(provenance, indent=1))
    (out / "DONE.json").write_text(json.dumps(dict(n_units=n_units, n_good=n_good,
                                                   finished=provenance["finished"])))
    if not keep_scratch:
        shutil.rmtree(str(scratch), ignore_errors=True)
    log.info("done: %d units (%d good) in %.1f min -> %s", n_units, n_good, (time.time() - t0) / 60, out)
    return dict(n_units=n_units, n_good=n_good, out=str(out))


def parse_kv(items):
    out = {}
    for it in items or []:
        k, v = it.split("=", 1)
        try:
            out[k] = json.loads(v)
        except ValueError:
            out[k] = v
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("monkey")
    ap.add_argument("day")
    ap.add_argument("experiment", type=int)
    ap.add_argument("--node", default=".", help="Record-Node folder relative to the day folder ('.' = flat)")
    ap.add_argument("--ks", action="append", metavar="KEY=VAL", help="Kilosort4 parameter override")
    ap.add_argument("--keep-scratch", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="only select recordings and report size")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    sort_session(args.monkey, args.day, args.experiment, args.node, parse_kv(args.ks),
                 args.keep_scratch, args.dry_run)


if __name__ == "__main__":
    main()
