# SpikeSort — spike sorting for all NHP Neuropixels recordings

Sorts every OpenEphys AP recording under `Data_collection/` with SpikeInterface + Kilosort4 and
writes cell-level spike data to `Data_analysis/SpikeSorting/` on the server (separate from the MUA
outputs in `Data_analysis/SaccadeDirection/`). Design/rationale: `PLAN.md`.

Runs on the Linux desktop (RTX 3060) that mounts `/media/NETDISKS`. Code is edited anywhere and synced by git.

## Setup (once, on the desktop)

Use the conda env that already has Kilosort4 (Python >= 3.10, CUDA torch). Add the missing pieces:

```bash
conda activate <your-kilosort-env>
pip install "spikeinterface[full]" probeinterface neo scipy pandas matplotlib h5py
# optional, only for UnitRefine model labels (config.UNITREFINE = True):
pip install skops huggingface_hub
python -c "import torch, kilosort, spikeinterface; print(torch.cuda.is_available(), kilosort.__version__, spikeinterface.__version__)"
```
Do **not** let pip replace your CUDA torch (check the check line prints `True`). `environment.yml` is only
for building a fresh env from scratch. Record the working versions under "Validation status" below.

Paths and parameters: `config.py` (server root, scratch dir `/media/DATA1/spikesort_scratch` (local 7.3 TB HDD; home NVMe too small for the largest days), Kilosort/curation settings).
Scratch needs roughly the size of the largest experiment's raw AP data free.

## Usage

```bash
cd SpikeSort
python find_recordings.py                 # 1. inventory -> Data_analysis/SpikeSorting/_inventory.csv
python run_batch.py --list                # 2. what would be sorted, sizes, DONE/todo
python sort_session.py m032 20260312 1 --dry-run   # 3. check selection for one session
python sort_session.py m032 20260312 1    # 4. pilot sort of one session
python run_batch.py                       # 5. everything not yet DONE (resumable)
```
- Non-flat layout (e.g. 20260224): `--node "m032-2026-02-24_12-59-28/Record Node 110"` (the inventory `node` column).
- Kilosort overrides: `--ks nblocks=1 --ks Th_universal=8`.
- Long runs: `nohup python run_batch.py > batch.log 2>&1 &` or inside `tmux`.
- A session is finished only when its `DONE.json` exists; failures are logged in `_manifest.csv` and `sort.log`, and the batch continues.

## Output layout

`Data_analysis/SpikeSorting/<monkey>/<day>/exp<N>/` (non-flat: `<node-slug>__exp<N>/`)

| File | Content |
|---|---|
| `spikes.mat` | MATLAB: `unit_id`, `times` {nUnits x 1} (s), `fs`, `segments`, `good`, `x_um`, `y_um` |
| `spikes.npz` | Python: flat arrays `spike_unit`, `spike_sample`, `spike_segment`, `spike_time_s` + segment table |
| `units.csv` | one row per unit: quality metrics, `x_um`/`y_um` (unit location), `peak_channel_*`, `good`, `ks_label`, optional UnitRefine labels |
| `analyzer/` | SpikeInterface SortingAnalyzer (templates, waveforms, metrics, correlograms) |
| `ks4_raw/` | Kilosort4 output (large `.dat` files removed) |
| `figures/` | `raster.png` (spike raster by depth — drift check), `metrics.png` |
| `sorting_params.json` | versions, all parameters, input recordings, segment table, bad channels, channel geometry |

**Sort unit:** all recordings of one experiment concatenated (recordings < 30 s dropped). Per-recording
boundaries are in `segments` so any run can be sliced back out.

**Timing:** spike times are seconds on the OpenEphys timeline = AP `timestamps.npy`, the same clock as the
NI-DAQ TTL `timestamps.npy` used in `SaccadeDirection/code/extractSaccadeSession.m`. Trial events and spikes can be
compared directly.

**`good` flag** = passes all thresholds in `config.CURATION` (amplitude cutoff, ISI violations, presence ratio, SNR). It marks
*candidate* single units; thresholds are Allen-style defaults, not yet validated on NHP data. Everything not `good` is
multi-unit/noise and still saved.

**Depth:** `y_um` is along the probe, same convention as `NP_PROBE.ELECTRODE_YPOS` used by the MUA pipeline, so units can
be mapped to the MUA channel index / anatomy table by nearest depth.

### Reading the results

```matlab
addpath(genpath('SpikeSort/matlab'));
U = loadSortedUnits('m032', '20260312', 1, 'goodOnly', true);
U.units(1:5, :)         % metrics + depth
spk = U.times{1};       % seconds, OE timeline
```
```python
import numpy as np, pandas as pd
d = "/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_analysis/SpikeSorting/m032/20260312/exp1/"
S, units = np.load(d + "spikes.npz"), pd.read_csv(d + "units.csv", index_col=0)
t = S["spike_time_s"][S["spike_unit"] == units.index[0]]
```

## Validation status (update after each step)

**Written offline, never run against real data.** The first run will likely need small fixes. Check in this order:
1. `find_recordings.py` runs and the inventory looks right (incl. 20260224 nested layout).
2. `sort_session.py ... --dry-run` selects the right recordings.
3. Pilot sort: geometry loads (`rec.has_probe()`), segment order matches recording order (script asserts this), `inter_sample_shift` ok, Kilosort runs, `raster.png` looks sane.
4. Test sparse-channel days ("single column", "quarter"): Kilosort may need `nblocks=1` or other parameters.
5. Compare sorted units with known MUA tuning (e.g. 20260323).
6. Validate `config.CURATION` / UnitRefine on ~100 hand-labelled units.

Working versions (desktop `kilosort` env, 2026-10-27): torch 2.5.1 (CUDA ok), kilosort 4.0.24, spikeinterface 0.103.2, probeinterface 0.3.1, neo 0.14.3.

## Known limitations
- Cross-session unit tracking not implemented (UnitMatch later).
- `ks4_raw/params.py` points at deleted scratch data; for Phy use `spikeinterface.exporters.export_to_phy` on the analyzer, or re-run with `--keep-scratch`.
- Sessions with gaps in `sample_numbers` are flagged in `sorting_params.json` (`contiguous: false`); times after a gap would be wrong.
