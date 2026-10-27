"""Central paths and parameters for the SpikeSort pipeline.

Runs on the Linux desktop (RTX 3060) that mounts the server at /media/NETDISKS.
Everything tunable lives here; per-session overrides go through CLI flags.
"""
import os
from pathlib import Path

SERVER_ROOT = Path("/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP")
DATA_COLLECTION = SERVER_ROOT / "Data_collection"             # raw: <monkey>/<YYYYMMDD>/...
OUTPUT_ROOT = SERVER_ROOT / "Data_analysis" / "SpikeSorting"  # spike-sorting results only
MANIFEST = OUTPUT_ROOT / "_manifest.csv"
INVENTORY = OUTPUT_ROOT / "_inventory.csv"

# Local scratch for the preprocessed binary + sorter working files (deleted after a
# successful sort). Needs ~ the size of the raw AP data of the largest experiment (333 GB).
# Home NVMe has only ~260 GB free; /media/DATA1 is a 7.3 TB local HDD (slower, but fits all).
SCRATCH = Path("/media/DATA1/spikesort_scratch")

# The server link is slow (~12 MB/s, see PLAN.md), so each experiment's raw AP data is copied
# once to local disk (stage_raw.py) and sorted from there. Deleted after a successful sort.
STAGE = True
STAGE_ROOT = Path("/media/DATA1/spikesort_raw")
STAGE_STREAMS = 2   # parallel rsync streams; one stream already saturates gigabit (~92 MB/s)
NET_IFACE = "eno1"  # checked before staging: the link has fallen back to 100 Mb/s after flaps
                    # (fix: `sudo ethtool -r eno1`)

AP_STREAM_KEY = "ProbeA-AP"
NIDAQ_STREAM_KEY = "NI-DAQmx"

# Recordings shorter than this (aborted tests) are left out of the sort.
MIN_RECORDING_S = 30.0

N_JOBS = max(1, (os.cpu_count() or 2) // 2)  # SpikeInterface chunked-processing workers

# --- preprocessing (SpikeInterface) -----------------------------------------------------
# phase_shift (NP1 ADC multiplexing) -> highpass -> remove bad channels -> median CMR.
# CMR is done here, so Kilosort's own CAR is switched off in KS4_PARAMS (avoid double CAR).
PREPROC = dict(
    highpass_hz=300.0,
    bad_channel_method="mad",   # "coherence+psd" (IBL, mouse-tuned) flagged all 384 ch on noisy NHP days (20260323)
    reference="median",   # "median" | "destripe" (spatial high-pass, as in XinyuScripts/getmua.m)
)

# --- Kilosort4 --------------------------------------------------------------------------
# nblocks>1 = non-rigid drift correction (awake NHP drifts). If a session has too few
# channels / too sparse geometry the retry ladder in sort_session.py falls back to nblocks=1.
KS4_PARAMS = dict(
    nblocks=5,
    do_CAR=False,
    batch_size=60000,
    Th_universal=9,
    Th_learned=8,
    torch_device="auto",
)
# Tried in order if the sort raises (e.g. CUDA OOM, too few channels for nonrigid).
KS4_FALLBACKS = [dict(nblocks=1), dict(nblocks=1, batch_size=30000)]

# --- curation ---------------------------------------------------------------------------
# A unit is flagged `good` (candidate single unit) if it passes ALL of these.
# Allen-style starting point; to be validated on hand-labelled NHP units (PLAN.md phase 4).
CURATION = dict(
    amplitude_cutoff_max=0.1,
    isi_violations_ratio_max=0.5,
    presence_ratio_min=0.9,
    snr_min=3.0,
)
# Also label units with the SpikeInterface UnitRefine models (needs internet + skops, and
# extra metrics -> slower). Models were trained on mouse; treat as advisory. UNTESTED.
UNITREFINE = False
