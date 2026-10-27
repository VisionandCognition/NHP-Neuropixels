"""Writers for sorted-unit output: Python-native (npz + csv) and MATLAB (.mat)."""
from __future__ import annotations

import numpy as np


def spikes_to_times(samples, seg_starts, seg_first_ts, fs):
    """Concatenated-recording sample index -> (segment index, seconds on the OpenEphys timeline).

    seg_starts   : start sample (in the concatenated recording) of each kept segment
    seg_first_ts : first `timestamps.npy` value (s) of each kept segment's AP stream, i.e. the
                   same timeline as the NI-DAQ TTL timestamps used by the MUA pipeline.
    """
    seg_starts = np.asarray(seg_starts, dtype=np.int64)
    seg = np.searchsorted(seg_starts, samples, side="right") - 1
    local = samples - seg_starts[seg]
    t = np.asarray(seg_first_ts, dtype=np.float64)[seg] + local / float(fs)
    return seg.astype(np.int16), t


def write_npz(path, unit_ids, spike_unit, spike_sample, spike_seg, spike_time_s, segments, fs):
    """Flat arrays, one entry per spike (sorted by sample). `spike_unit` holds unit ids."""
    np.savez(path,
             unit_ids=np.asarray(unit_ids, dtype=np.int64),
             spike_unit=np.asarray(spike_unit, dtype=np.int64),
             spike_sample=np.asarray(spike_sample, dtype=np.int64),
             spike_segment=np.asarray(spike_seg, dtype=np.int16),
             spike_time_s=np.asarray(spike_time_s, dtype=np.float64),
             fs=float(fs),
             seg_n_samples=np.asarray([s["n_samples"] for s in segments], dtype=np.int64),
             seg_first_timestamp_s=np.asarray([s["first_timestamp_s"] for s in segments], dtype=np.float64),
             seg_first_sample_number=np.asarray([s["first_sample_number"] for s in segments], dtype=np.int64),
             seg_recording=np.asarray([s["recording"] for s in segments]))


def write_mat(path, units_df, unit_ids, spike_unit, spike_time_s, segments, fs):
    """MATLAB v5 .mat (readable everywhere, 2 GB limit). Spike times per unit as a cell array."""
    from scipy.io import savemat
    spike_unit = np.asarray(spike_unit)
    n_spikes = len(spike_unit)
    if n_spikes * 8 > 1.5e9:
        print("WARNING: %d spikes; spikes.mat may exceed the 2 GB v5 limit" % n_spikes)
    times = np.empty((len(unit_ids), 1), dtype=object)
    for i, u in enumerate(unit_ids):
        times[i, 0] = spike_time_s[spike_unit == u].reshape(-1, 1)
    seg_struct = {
        "recording": np.asarray([s["recording"] for s in segments], dtype=object),
        "n_samples": np.asarray([s["n_samples"] for s in segments], dtype=np.float64),
        "first_timestamp_s": np.asarray([s["first_timestamp_s"] for s in segments], dtype=np.float64),
        "first_sample_number": np.asarray([s["first_sample_number"] for s in segments], dtype=np.float64),
        "contiguous": np.asarray([int(s["contiguous"]) for s in segments], dtype=np.float64),
    }
    mat = {
        "unit_id": np.asarray(unit_ids, dtype=np.float64).reshape(-1, 1),
        "times": times,                       # {nUnits x 1} cell, seconds, OE timeline
        "fs": float(fs),
        "segments": seg_struct,
        "good": units_df["good"].to_numpy(dtype=np.float64).reshape(-1, 1),
        "y_um": units_df["y_um"].to_numpy(dtype=np.float64).reshape(-1, 1),
        "x_um": units_df["x_um"].to_numpy(dtype=np.float64).reshape(-1, 1),
    }
    savemat(path, mat, do_compression=True)
