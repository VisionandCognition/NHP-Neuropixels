# SpikeSort pilots — what we test and why

Running log (newest last). Results table: `pilots/results.csv`; per-variant unit tables:
`<scratch>/pilots/<key>/units_<variant>.csv`. Code: `pilot_slice.py`, queue: `run_pilots.py`.

## Context (2026-10-27)

- Spike sorting is a **bonus**; the main analysis uses MUA envelopes. Take compute time for quality,
  but don't over-invest in the hardest sessions.
- Pipeline ran end-to-end on `m032 20251010 exp1` (14.8 GB, Single Column): 402 units, 46 good, 38.7 min.
  Raster stable; no bad channels. Kilosort ~1x real time, of which drift correction 267 of 648 s.
  **Drift correction stays on (`nblocks=5`) by default**: drift over a whole session can be substantial
  even if single runs look stable (user). A 10-min raster says little about a 2.7 h session.
- Network link: negotiated 100 Mb/s after link flaps (12 MB/s). `sudo ethtool -r eno1` restored 1000 Mb/s
  (92 MB/s). `stage_raw.check_link()` warns if it drops again.
- Electrode configuration varies a lot per session (coverage prioritised over density on the long NHP probe).
  `geometry.py` classifies: `dense_bankA` (20260915, 20260918), `single_column` (10 days),
  `four_bank_even` (20251022, 0310, 0326, 0819, 0904), `custom_*` (8 days).

## Structure of the sparse layouts (from settings.xml)

| Type | Days | Structure |
|---|---|---|
| four_bank_even | 20251022, 0310, 0326, 0819, 0904 | one continuous line, uniform 40 um pitch, 15.3 mm |
| custom, big patches | 20260224, 0226, 0303, 0317, 0318, 0320 | 2-3 patches of ~50-220 ch (40 um pitch) separated by 1-11 mm, plus a few 1-2 ch strays |
| custom, fragmented | 20260324, 0325 | 38 / 51 patches, mixed 20-80 um pitch (0325: 276 of 365 ch in 4 patches >= 20 ch; 97 ch in 46 tiny patches) |

## Open question 1: Kilosort `dminx`

Installed Kilosort 4.1.7 (`kilosort/parameters.py`): `dmin` is determined automatically; `dminx` defaults
to 32 um ("works well for Neuropixels 1 and 2"; for other probe geometries "try setting this to the median
lateral distance between contacts"). On the NHP probe the two columns are at x = 11 and 114 um (~103 um
apart). Hypothesis: 32 um places redundant template centres in x and may contribute to over-splitting
(402 units on 384 channels in the exp1 run). **Not verified**; this could also affect already-sorted sessions.

## Other design choices under test

- **Drop strays**: patches with < 8 channels cannot support template matching; drop them for sorting
  only (they remain in MUA analyses).
- **Sort per patch**: each contiguous patch is its own "probe" (own drift estimate); patches < 20 ch skipped,
  patches < 150 ch use `nblocks=1`. Costs more compute and bookkeeping.
- Other sorters (SpyKING CIRCUS 2, Tridesclous2, MountainSort5) adapt to irregular geometry but no
  benchmarks on irregular NHP layouts are known to us; kept as comparator only if Kilosort looks poor.

## Queue (3-min slices, base params `nblocks=1` to isolate geometry effects; production keeps 5)

| # | Day / profile | Variants | Question |
|---|---|---|---|
| T1 | 20260323 exp1, single_column (reuses smoke-test slice) | default, dminx103 | does `dminx` change yield? |
| T2 | 20260819 exp1, four_bank_even | default, dminx103 | do regular sparse lines behave on defaults? |
| T3 | 20260317 exp1, custom patches (221/84/56 ch) | default, dminx103, drop_strays, drop_strays_dminx103, by_patch | strays / per-patch handling |
| T4 | 20260325 exp1, fragmented | default, drop_strays_dminx103, by_patch | is anything usable there? |

How to judge: unit count, good fraction, median SNR, ISI-violation fraction, spikes per channel,
runtime. More units is not better by itself (over-splitting); prefer fewer, cleaner units with the same
or higher good count. Final judgement by eye (raster / Phy) on the most promising variant.

## Log

- 2026-10-27: pilots queued behind `20251010 exp2 --phy` (Phy-check run; GPU busy until it finishes).
- 2026-10-27 18:42: `m032 20251010 exp2 --phy` finished: 526 units (62 good), 4.2 M spikes, 62.6 min for 32.6 min of
  data (copy 8.2 min at gigabit, Phy export ~8 min). Phy folder (43 GB incl. recording.dat) is on the server and
  loads headlessly (phylib: 526 clusters, waveform read OK). No PC features in the export, so Phy's feature view is
  empty. Pilot queue started automatically right after.
  Open in Phy on the desktop (DWService): `cd <out>/phy && /home/chris/anaconda3/envs/phy2/bin/phy template-gui params.py`

## Pilot results (2026-10-27, 3-min slices, nblocks=1; `pilots/results.csv`)

| Slice | Variant | Units | Good | Spikes | KS runtime |
|---|---|---|---|---|---|
| 20260323 single column | default (dminx 32) | 416 | 98 | 824 k | 238 s |
| | dminx103 | 328 | 93 | 774 k | 135 s |
| 20260819 four_bank_even | default | 399 | 76 | 404 k | 192 s |
| | dminx103 | 335 | 75 | 385 k | 118 s |
| 20260317 custom patches | default | 385 | 52 | 248 k | 251 s |
| | dminx103 | 368 | 51 | 241 k | 149 s |
| | drop_strays (2 ch only) | 381 | 52 | 246 k | 253 s |
| | drop_strays + dminx103 | 357 | 49 | 240 k | 149 s |
| | by_patch | FAILED first time (quiet 56-ch patch: "n_samples=5 should be >= n_clusters=6"); code now skips failed patches; rerun pending |
| 20260325 fragmented | default (384 ch, 51 patches) | 375 | 82 | 373 k | 227 s |
| | drop_strays + dminx103 (287 ch, 5 patches) | 305 | 74 | 351 k | 105 s |
| | by_patch (287 ch) | 302 | 78 | 350 k | 102 s |

Kilosort's auto `dmin` = 20 um on the single column (20 um pitch) and 40 um on the 40-um-pitch layouts; `dminx` as set.

**Agreement check default vs dminx103 (`compare_two_sorters`, delta 0.4 ms, match >= 0.5):**
- 90 / 92 / 98 % of the default `good` units are re-found with dminx103 (T1 / T2 / T3), and 86 / 65 / 48 of those matches are also `good`.
- 19-31 % of default units have no counterpart in dminx103 (131/416, 103/399, 75/385), but only 10 / 6 / 1 of them were `good`.
- 43 / 40 / 61 dminx103 units have no counterpart in default; only 5 / 6 / 0 are `good`.

**Reading:** dminx=103 keeps nearly all well-isolated units, drops ~15-20 % of units that are almost all
low-quality, and cuts Kilosort time ~40 %. Consistent with the over-splitting hypothesis (32 um places
redundant template centres on this probe), but NOT proof: a lost non-good unit could be a real low-amplitude
neuron. Needs a look in Phy (compare matched vs unmatched units) before changing the production default.

**Caveats:** `good` uses weak thresholds (SNR >= 3 etc., unvalidated on NHP); median SNR is only ~2.6-3.5 and
only 12-25 % of units pass, so most output is MUA-like either way. 3-min slices say nothing about drift.
Dropping the 2 strays on 20260317 had no effect (only 2 of 384 channels). T4 conflates stray-dropping
with dminx (no drop-only variant); the 25 % channel loss on 20260325 (97 ch in 46 tiny patches) cost only
8 of 82 good units and halved runtime.
