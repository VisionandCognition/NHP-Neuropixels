# SpikeSort — generic spike-sorting pipeline plan

Goal: spike-sort **all** NHP Neuropixels recordings (so far only reduced to MUA envelopes, e.g. `SaccadeDirection/`), and save cell-level spike data on the server in a clearly labelled folder next to the existing MUA outputs.

Status/progress log lives in the Claude memory file `spikesorting-effort` (update it every session).

## 1. What the existing pipeline tells us (facts from `SaccadeDirection/`)

- Server root: `/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/` (Linux mount; **not mounted on this Mac**, so nothing below is verified against real files yet).
  - Raw: `Data_collection/<monkey>/<YYYYMMDD>/` (monkey e.g. `m032`)
  - MUA/derived: `Data_analysis/SaccadeDirection/data/` (`saccadeDataDir.m`)
- OpenEphys **binary** format. Recording layout is NOT uniform: usually `<day>/experimentN/recordingM/`, but e.g. 20260224 is nested under `<day>/m032-2026-02-24_12-59-28/Record Node 110/experiment1/recordingM`. `matchRunRecording.m` already handles this with a recursive `recording*` search — reuse that logic (port to Python).
- Streams (from `structure.oebin`): `Neuropix-PXI-100.ProbeA-AP` (30 kHz, 384 ch) and `NI-DAQmx-114.PXIe-6341` (TTL + eye analog). `settings.xml` (`NP_PROBE.ELECTRODE_YPOS`) gives channel→electrode geometry.
- **Online channel selection varies per day** ("quarter density", "4bankeven", "single column", "Bank A", ...). So geometry differs by day; sparse selections (single column, quarter) are poor inputs for Kilosort's template matching → needs pilot testing.
- One OE recording can span several behavioral runs; a day has several recordings (e.g. recording1..4). Raw-data layout is shared with older 2024/2025 sessions (other tasks, e.g. curve-tracing) — all in scope for this pipeline.
- NP1.0-style probe (12-way ADC multiplexing → `phase_shift` needed).
- MUA pipeline timing: events are aligned via NIDAQ TTL timestamps (`timestamps.npy`) vs AP `sample_numbers`. Spike output must keep both so MATLAB code can align exactly like `extractSaccadeSession.m`.

## 2. Tool choices (modern stack, as of 2026)

| Stage | Tool | Why |
|---|---|---|
| Orchestration / IO / preprocessing | **SpikeInterface** (`si.read_openephys`, `run_sorter`, `SortingAnalyzer`) | Standard, reads OE binary, sorter-agnostic, built-in QC/curation/export |
| Sorter | **Kilosort4** (primary) | Current standard for NP, built-in nonrigid drift correction, GPU-fast. Fallback: SpykingCircus2 / Mountainsort5 (CPU) for sessions with sparse geometry where KS4 fails |
| Preproc | `phase_shift` → `highpass_filter(300–400 Hz)` → `detect_bad_channels` (remove) → `common_reference(median)` or `highpass_spatial_filter` (destripe; analog of what `getmua.m` already does) | IBL-style, matches the MUA chain so MUA and spikes are comparable |
| Drift | KS4 internal (`nblocks`) by default; SI `correct_motion` (DREDge) only if drift raster shows large drift | Awake NHP = more drift than mouse; check raster per session and store it |
| QC metrics | SI quality metrics (SNR, ISI-violation ratio, amplitude cutoff, presence ratio, firing rate) | Reproducible thresholds |
| Curation | Stage A: SI **UnitRefine** model-based labels (noise vs neural, SUA vs MUA) + threshold rule (Allen-style). Stage B (optional): Bombcell cross-check; Phy for spot-check of key sessions | Models were trained on mouse; **validate on manual NHP subset** before trusting |
| Cross-session tracking | later: **UnitMatch** (only relevant for same-site repeat penetrations; not v1) | |
| Env | Python ≥3.10 (Mac has 3.9.6, no SI/torch), conda/uv env on the GPU server; pin versions (`spikeinterface`, `kilosort`, `probeinterface`, `neo`) | Reproducibility |

Hardware (confirmed 2026-10-27): runs on the user's Linux Mint desktop (RTX 3060, 12 GB; same machine as the saccade analysis, mounts `/media/NETDISKS`). Mac is for editing/dry runs only. 12 GB VRAM is fine for KS4 on 384 ch, but lower `batch_size` if OOM on long concatenated sessions; run one sort at a time.

## 3. Sorting unit of work

**One sort per (monkey, day, experiment)**, concatenating all recordings of that experiment that share the same probe config (same `settings.xml`/channel map), *not* per behavioral run. Rationale: units are stable within a day; sorting per run would fragment units and the task-agnostic goal ("all data") implies not tying sorts to one paradigm. Segment boundaries are stored so any run/recording can be sliced back out. Recordings with a different channel selection or a probe re-seat → separate sort.

## 4. Output location & format

```
Data_analysis/SpikeSorting/<monkey>/<YYYYMMDD>/exp<N>/
  sorting_params.json        # versions, all parameters, input recording paths, git-free provenance
  ks4_raw/                   # untouched Kilosort4 output (phy-compatible)
  analyzer/                  # SpikeInterface SortingAnalyzer (templates, metrics) 
  units.csv                  # per-unit: id, ks label, UnitRefine label, metrics, depth (um), channel, area (from anatomy table where available), keep flag
  spikes.mat                 # MATLAB-friendly v7 (NOT 7.3 if size allows): per-unit spike times (s, OE timeline) + sample indices, segment boundaries, first timestamp, fs, unit table
  figures/                   # drift raster, unit summary sheets, QC overview
Data_analysis/SpikeSorting/_manifest.csv   # one row per (day, exp): status, sorter, date, #units, #good, errors
```
Assumption: `SpikeSorting/` sits beside `SaccadeDirection/` under `Data_analysis/`. (Rationale: generic, not task-specific.) Easy to change in one place (`config.py`).

MATLAB side: small `loadSortedUnits.m` (+ `spikesortDataDir.m` mirroring `saccadeDataDir.m`) so saccade code can later swap MUA → spikes.

## 5. Code layout (to build in `SpikeSort/`)

```
SpikeSort/
  PLAN.md                  # this file
  config.py                # server roots, sorter params, curation thresholds
  find_recordings.py       # port of matchRunRecording logic: inventory all days/experiments/recordings -> inventory.csv
  sort_session.py          # one (day, exp): load -> preprocess -> sort -> analyzer -> metrics -> curate -> export
  run_batch.py             # loop over inventory, skip-if-done via manifest, resumable, per-session log
  export_matlab.py         # spikes.mat writer
  matlab/loadSortedUnits.m, matlab/spikesortDataDir.m
  environment.yml
```

## 6. Phases

0. **Inventory** (needs server): list every monkey/day/experiment/recording, probe config (channel count, selection scheme from `settings.xml`/oebin), duration, size. Output `inventory.csv`. Decide which sessions are in scope.
1. **Env**: Python ≥3.10 env on the GPU machine; confirm CUDA, KS4 import, server read access.
2. **Pilot** on 2–3 contrasting sessions (full-density chamber-3 day, a "single column" day, a chamber-2 subcortical day): check geometry loads right, drift raster, KS4 runs, unit yield vs MUA.
3. **Validate**: compare against existing MUA results (e.g. do sorted units on a known tuned channel show the same direction tuning? depth match?). Tune KS4 params (`Th_universal/learned`, `nblocks`, `batch_size`), decide preprocessing for sparse-channel days.
4. **Curation**: define thresholds + UnitRefine; hand-label ~100 NHP units to check agreement; freeze rules, store in `config.py`.
5. **Batch**: run all in-scope sessions via `run_batch.py` (resumable, manifest). Expect long GPU runtimes — run in background / tmux.
6. **MATLAB bridge + integration**: `loadSortedUnits.m`; map unit depth→area using `getChannelArea.m` / anatomy CSV; first re-run of a saccade analysis with SUA.

## 7. Decisions from user (2026-10-27)

- Machine: Linux Mint desktop, RTX 3060, mounts server (see §2).
- Output location `Data_analysis/SpikeSorting/` OK; code lives in `./SpikeSort`.
- Scope: **every recording with data**, all sessions/tasks (no task filter). Currently only m032 exists, but the pipeline must discover all monkeys under `Data_collection/` (config-free: any `m*` subject folder).
- Formats: both MATLAB and Python. => `spikes.mat` for MATLAB, plus Python-native `spikes.npz`/`units.csv` and the SI analyzer folder (zarr/binary). NWB not requested; skip unless wanted later.

## 8. Still open

- Is raw data ever only on a local disk (`getmua.m` has an `E:\OpenEphys` local option)? Inventory will report recordings it cannot find a `continuous.dat` for.
- Disk budget on server for `ks4_raw` + preprocessed binary (preprocessed copy ~ same size as raw; write to local scratch, delete after sort).
