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
  Open in Phy on the desktop (DWService): `./phy.sh [phy folder]`. Plain `phy template-gui` segfaulted (QtWebEngine GPU
  thread, `Chrome_InProcGp` in libQt5WebEngineCore, 2026-10-27 20:31); `phy.sh` sets `QTWEBENGINE_CHROMIUM_FLAGS="--disable-gpu --no-sandbox"`.

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
| | by_patch (rerun, failed patches skipped) | 269 | 39 | 198 k | 130 s | (84-ch patch failed: Kilosort could not initialise templates, "n_samples=5 should be >= n_clusters=6"; 22 % of channels lost) |
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

**Per-patch sorting: not worth it.** On 20260317 it lost a whole 84-channel patch (Kilosort crashed on it)
and gave 39 good units vs 49-52 for whole-probe sorting; on 20260325 it matched drop_strays+dminx103 (78 vs 74
good). Whole-probe Kilosort handles the mm-wide gaps fine, so no extra machinery. Failed-patch handling kept in
pilot code only (not in production `sort_session.py`).

**Proposed production changes (pending the Phy look):** `KS4_PROFILES` dminx=103 for all profiles; drop
patches < 8 channels for sorting (matters on fragmented days like 20260324/0325); keep `nblocks=5`.

## Phy-output inspection of 20251010 exp2 (2026-10-28) — noisy channel bug

Headless check of `exp2/phy` (script logic: units.csv + phy template files + recording.dat):
- **349 of 526 units (66 %), 2.35 M of 4.2 M spikes, 0 good, median SNR 0.5, all peak on ONE channel (AP191, y=5740 um, x=103).**
  That channel has rms 1325 vs probe median 151 (9x, highest of 384). The other 177 units sit on 115 channels (max 7/channel),
  62 good, median SNR 3.2. So the "526 units / 62 good" headline is mostly Kilosort templating noise on one channel.
- Cause: `preprocess` ran `detect_bad_channels(mad)` before the median reference; shared noise dominates there
  (ch191 only 2x the median std) so nothing was flagged ("bad_channel_ids": []). After referencing it is flagged
  (also at lower thresholds ch300). **Fix applied in `sort_session.preprocess`: detect on the median-referenced signal.**
  Not yet re-run. exp1 and the 3-min pilots used the old preprocessing, so they may carry the same problem: the
  default-vs-dminx103 unit counts (and over-splitting reading) above are suspect until redone with the fix.
- Real over-splitting among the 177 clean-channel units looks modest: 15 candidate merge pairs (similar template,
  CCG dip, merged ISI violation < 1 %) involving 28 units, 1 of them good.
- Of the 62 good units, 14 have ISI violation ratio > 0.2 (thresholds are loose).

### Re-sort of 20251010 exp2 with fixed bad-channel detection (2026-10-28)
- Preprocessing now: detect on median-referenced signal, `std_mad_threshold=3` (SI default 5 missed AP300 at 3.4x
  median std), plus drop channels with std > 3x probe median (`config.PREPROC`). Removed AP191 and AP300.
  Cutoff 3 flags 0 channels on the 4 other cached pilot days (incl. noisy 20260323); the 3x-median rule would flag 1 on 20260819.
- Result: 518 units, 76 good (was 526 / 62), 50.9 min. 425 units on 223 channels (max 8 per channel, 74 good);
  93 units at the probe ends (y <= 20 or >= 7600 um; 57 at y=7600) hold only 2 good -> edge artefacts, low priority.
- Total spike count stayed ~4.2 M in all three sorts (old, AP191 removed only, both removed): Kilosort's whitening
  rescales each channel, so a noisy channel still yields about the same number of detections; the damage was
  unit labels/clutter, not lost real units.
- sort_session: analyzer now built on local scratch then copied (memmap flush on SMB intermittently gave
  PermissionError); `--resume` reuses scratch preprocessed binary + Kilosort output; an open Phy window holding
  the share's phy/ folder blocks overwriting it (close Phy before re-running with --phy).
- Pilot slices (3 min) and exp1 used the old preprocessing -> the dminx comparison must be redone. 

### Pilot rerun with fixed preprocessing: 20260323 single column (2026-10-28)
New preprocessing removed 0 channels on this slice, so the old pilot was not affected by the noisy-channel bug.
| Variant | Units | Good | Spikes | KS runtime |
|---|---|---|---|---|
| default (dminx 32) | 405 | 95 | 790 k | 220 s |
| dminx103 | 322 | 95 | 746 k | 135 s |
(old: 416/98 and 328/93.) Same picture: -20 % units, same good count, -40 % runtime.
Other cached pilot days: top channel std/median was 1.6-1.9 on 0317, 0325 and smoke(0323); 0819 has one 3.1x channel
(would now be dropped). So only the 0819 pilot is expected to change; 0317/0325 results stand.

### Production defaults applied + exp1 re-sort (2026-10-28)
Applied: `KS4_PROFILES` dminx=103 for all profiles (dense_bankA unpiloted), `PREPROC.min_patch_channels=8` (stray drop, logs
`stray_channel_ids`), bad-channel rules above, nblocks=5 kept.
`20251010 exp1` re-sort: 343 units, 47 good, 1.26 M spikes, 12.1 min total (staging 2.7 min; old run 402 / 46 good, 38.7 min
incl. 21.7 min copy at 100 Mb/s). Bad channels: AP191 and AP300 again (same two as exp2 -> probably permanently noisy on this
probe/day; worth checking other days). Stray-drop path not exercised yet (single column). 25 units still peak at y=7600 (probe tip).

### run_batch test (2026-10-28): --until, prefetch, STOP, stray drop
Ran `--only 20260303:1 --only 20260326:1 --until <+3h>`: both OK. Prefetch of 0326 (8 min copy) overlapped the 0303 sort; its
own staging then took 0.0 min. Stray drop worked (0303: 7 channels in patches < 8 ch). Staged raw removed after each session.
- 20260303 exp1 (custom_sparse, 34 min): 429 units, only 3 good (flag `low_good_frac`), median SNR 0.9, 27 units with SNR >= 3 -> low-yield day, flagged as intended.
- 20260326 exp1 (four_bank_even, 38 GB): 165 units, 52 good, 21.5 min total (~0.57 min/GB; estimate of 1.2 min/GB is conservative).
Full chronological batch started afterwards (no --until; stop with `SpikeSort/STOP`).

### Batch status 2026-10-29 07:55 and staging collision bug
Done so far: 20251010 exp3 (721 units/87 good, 179 min), 20251016 exp1 (582/106), 20251022 exp1 (749/138, four_bank_even),
20260224 Record Node 110 exp1 (434/11 good, custom_sparse, 59 min), plus earlier tests 0303 exp1, 0326 exp1, 20251010 exp1/exp2.
`many_units` flag (> 1.5 units/channel) on the three big single-column/four-bank sessions; counts are 1.5-2 units/ch with
87-138 good units, so the threshold is probably too tight, not a failure (unreviewed).
FAILED: 20260224 exp1 (flat) and 20260226 exp1: staged_root put all experiments of a day (and Record-Node subfolders) in one
folder, so the prefetch of the next session (same day) was being copied while neo scanned that folder for the current one
("mmap length is greater than file size" / missing timestamps.npy). Fix: `stage_raw.staged_root` = day/exp<N>[__<node>], one folder
per session (siblings, no nesting). The running batch (old code) was asked to stop after its current session (STOP file) to be
restarted with the fix; failed sessions have no DONE.json so the restart picks them up again.

### Batch review 2026-10-30: AP191 and 20260312
- **AP191 is a bad ADC channel on most days**, independent of electrode config: units peaking there per session: 1010 exp2 349 (pre-fix),
  0310 106, 0317 91, 0303 exp3 70, 0326 57, 0226 exp2 47, 0303 exp2 43, 0226 exp1 25, 0224 17, 0312 17, 0303 exp1 10 (clean: 1010 exp3, 1022).
  On 0310 these units have 6x the amplitude of others, 0.04 Hz, ISI-violation ratio ~270, i.e. rare large transients; the channel's
  median-based noise is normal (0.7x), so the noise rules do not catch it. `PREPROC.always_bad=("AP191",)` now always drops it.
  AP300 was noisy only on 20251010. Sessions already done keep their AP191 clutter units (flaggable post hoc via peak_channel_id).
- `max_units_per_channel` 1.5 -> 2.5 (flag fired on nearly every session).
- **20260312 exp1** (1539 units, 26 good): not noise. 159 min in 15 recordings; spike rate stable but the unit population turns over:
  only 13 % of units active (>= 20 spikes) in every recording (other days 32-66 %), 30 % in <= 3 recordings. Kilosort's drift
  correction does not bridge abrupt changes between recordings, so units are split by time -> low presence_ratio (median 0.53).
  Cure would be per-recording sorting + unit matching; not done. Treat this day's units with care; MUA unaffected.
