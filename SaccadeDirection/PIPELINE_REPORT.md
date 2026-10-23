# Saccade Direction Pipeline — Report

Subject: m032 ("Kwibus"). Two recording chambers: **chamber 2** (12 sessions, 2026-02-24 to 2026-03-26; posterior burr holes, trajectories mostly subcortical) and **chamber 3** (10 sessions, 2026-08-19 to 2026-09-18; prefrontal burr holes, cortical trajectories). Task: 8-target visually-guided saccade task, run under the shared curve-tracing ("TDCT"/"ct") behavioral codebase with `Stm.difficulty = 1` (single target, no curve, no distractor) and stimulus file `stim_ct_classic_8targ.m`. Data source: `/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_collection/m032`.

All code lives in `SaccadeDirection/code/`. All extracted data and figures live in `SaccadeDirection/data/` and `SaccadeDirection/data/figures/`.

---

## 1. Objective

Extract raw Neuropixels + eye-tracking data for every saccade-direction session, reduce the probe data to a multi-unit activity (MUA) envelope per channel per trial, screen every channel for task responsiveness, and compute saccade-direction tuning for responsive channels — with the eye-tracking data used both as a behavioral readout in its own right (fixation stability, saccade kinematics, pupil) and as the ground truth for calibrating and QC-ing the neural alignment.

## 2. Identifying the task and the sessions

A run is a saccade-direction run if and only if its folder contains `stim_ct_classic_8targ.m`. Scanning all of `m032`'s day folders for this marker found **22 runs**, one per day, in two chambers: 12 in chamber 2 (2026-02-24 to 2026-03-26) and 10 in chamber 3 (2026-08-19 to 2026-09-18, all `run-002`; their stimulus file is byte-identical to the March one).

| Day | Run | Chamber | Day | Run | Chamber |
|---|---|---|---|---|---|
| 20260224 | run-001 | 2 | 20260326 | run-003 | 2 |
| 20260226 | run-003 | 2 | 20260819 | run-002 | 3 |
| 20260303 | run-003 | 2 | 20260821 | run-002 | 3 |
| 20260310 | run-003 | 2 | 20260825 | run-002 | 3 |
| 20260312 | run-003 | 2 | 20260903 | run-002 | 3 |
| 20260317 | run-003 | 2 | 20260904 | run-002 | 3 |
| 20260318 | run-003 | 2 | 20260909 | run-002 | 3 |
| 20260320 | run-003 | 2 | 20260910 | run-002 | 3 |
| 20260323 | run-003 | 2 | 20260915 | run-002 | 3 |
| 20260324 | run-003 | 2 | 20260916 | run-002 | 3 |
| 20260325 | run-003 | 2 | 20260918 | run-002 | 3 |

The stimulus file confirms this is a genuine 8-direction task: `Stm.Targ(1..8).xy_dva` places targets at 0°, 45°, 90°, ..., 315° at 12 dva eccentricity, with `Stm.difficulty = 1` disabling the curve/distractor machinery inherited from the curve-tracing paradigm, leaving a plain visually-guided saccade: fixate → brief target flash (100ms) → 400ms memory delay → go-cue → saccade.

Older 2024/2025 sessions in the same folder tree use the same task family (`task-ct`) but only the curve-tracing variants (`fgdots`/`classic` without `_8targ`), and are excluded.

### 2.1 Run → OpenEphys recording → behavioral log mapping

Each behavioral run folder is matched to its OpenEphys `experimentN/recordingM` folder by `matchRunRecording.m` (closest absolute start time between the run's folder timestamp and the recording's `sync_messages.txt` "Software Time," within a 90s tolerance) and to its own behavioral log `.mat` by `findSaccadeRuns.m` (the `sub-*.mat` file inside the run folder). Recording folders are searched recursively under the day folder, since they are not always in the default `experimentN/recordingM` location directly under it — on 20260224, for instance, the recordings are nested inside an additional `m032-2026-02-24_12-59-28/Record Node 110/` folder. The verified pairing for all 22 sessions:

| Day | Behavioral run folder | Behavioral log (`.mat`) | OE recording folder | Trials |
|---|---|---|---|---|
| 20260224 | `run-001_T-134541` | `sub-m032_ses-20260224_task-ct_stim-classic_run-001.mat` | `m032-2026-02-24_12-59-28/Record Node 110/experiment1/recording4` | 85 |
| 20260226 | `run-003_T-132430` | `sub-m032_ses-20260226_task-ct_stim-classic_run-003.mat` | `experiment2/recording3` | 145 |
| 20260303 | `run-003_T-140041` | `sub-m032_ses-20260303_task-ct_stim-classic_run-003.mat` | `experiment2/recording3` | 114 |
| 20260310 | `run-003_T-143306` | `sub-m032_ses-20260310_task-ct_stim-classic_run-003.mat` | `experiment1/recording4` | 117 |
| 20260312 | `run-003_T-132023` | `sub-m032_ses-20260312_task-ct_stim-classic_run-003.mat` | `experiment1/recording4` | 85 |
| 20260317 | `run-003_T-125341` | `sub-m032_ses-20260317_task-ct_stim-classic_run-003.mat` | `experiment1/recording4` | 88 |
| 20260318 | `run-003_T-134537` | `sub-m032_ses-20260318_task-ct_stim-classic_run-003.mat` | `experiment1/recording4` | 85 |
| 20260320 | `run-003_T-121255` | `sub-m032_ses-20260320_task-ct_stim-classic_run-003.mat` | `experiment1/recording4` | 100 |
| 20260323 | `run-003_T-112901` | `sub-m032_ses-20260323_task-ct_stim-classic_run-003.mat` | `experiment1/recording4` | 95 |
| 20260324 | `run-003_T-121719` | `sub-m032_ses-20260324_task-ct_stim-classic_run-003.mat` | `experiment1/recording4` | 99 |
| 20260325 | `run-003_T-122305` | `sub-m032_ses-20260325_task-ct_stim-classic_run-003.mat` | `experiment1/recording4` | 99 |
| 20260326 | `run-003_T-125258` | `sub-m032_ses-20260326_task-ct_stim-classic_run-003.mat` | `experiment2/recording3` | 101 |
| 20260819 | `run-002_T-141404` | `sub-m032_ses-20260819_task-ct_stim-classic_run-002.mat` | `experiment1/recording3` | 102 |
| 20260821 | `run-002_T-134629` | `sub-m032_ses-20260821_task-ct_stim-classic_run-002.mat` | `experiment1/recording3` | 102 |
| 20260825 | `run-002_T-121203` | `sub-m032_ses-20260825_task-ct_stim-classic_run-002.mat` | `experiment1/recording3` | 115 |
| 20260903 | `run-002_T-135647` | `sub-m032_ses-20260903_task-ct_stim-classic_run-002.mat` | `experiment1/recording3` | 129 |
| 20260904 | `run-002_T-125528` | `sub-m032_ses-20260904_task-ct_stim-classic_run-002.mat` | `experiment1/recording3` | 128 |
| 20260909 | `run-002_T-124444` | `sub-m032_ses-20260909_task-ct_stim-classic_run-002.mat` | `experiment1/recording3` | 110 |
| 20260910 | `run-002_T-141659` | `sub-m032_ses-20260910_task-ct_stim-classic_run-002.mat` | `experiment1/recording3` | 107 |
| 20260915 | `run-002_T-125530` | `sub-m032_ses-20260915_task-ct_stim-classic_run-002.mat` | `experiment1/recording3` | 106 |
| 20260916 | `run-002_T-125529` | `sub-m032_ses-20260916_task-ct_stim-classic_run-002.mat` | `experiment1/recording3` | 117 |
| 20260918 | `run-002_T-131444` | `sub-m032_ses-20260918_task-ct_stim-classic_run-002.mat` | `experiment1/recording3` | 101 |

All 10 chamber-3 runs map to `experiment1/recording3`, which the elab log independently lists as the recording number of that day's 8targ run. Each `_extracted.mat` stores this verified path as `S.LogFile`. **Do not use `S.Par.LogFn`/`S.Par.LogPath` for this purpose** — those fields come from inside the log file's own saved `Par` struct and are stale in every session checked (they point at an unrelated earlier `task-rfmap`/`stim-sparsenoise` run in the same stim-program session, not the 8-target run the data was actually extracted from).

### 2.2 Burr hole, channel selection, and planned target structure per session

**Chamber 2 (12 sessions).**

`Chamber2Holes.docx` (MRI-based trajectory planning) documents 5 burr holes (0, 1, 2, 4, 5), each with a planned insertion depth and a list of structures its trajectory passes through, superficial to deep. `Kwibus_elab.pdf` (the day-by-day lab notebook) gives the burr hole and online channel-selection scheme actually used on each recording day (e.g. "quarter density," "4bankeven," "single column," "sections2/3" — the operator manually restricting acquisition to part of the probe based on what looked good online; this is a per-bank/section choice, not per-electrode). Cross-referencing the two, all 12 chamber-2 sessions:

| Day | Burr hole | Channel selection (online) | Planned structures (superficial → deep) |
|---|---|---|---|
| 20260224 | 5? *(uncertain in source log)* | Ch2H5quarter | S1 BA1/2, BA3, PCC BA23, retrosplenial BA29 → Habenula, Posterior commissure or CM |
| 20260226 | 4 | Ch2H4quarter | S1 BA1/2, BA3, (M1), PCC BA23 → Thalamus MD, CL, CM |
| 20260303 | 5 | Ch2h5quarterD | S1 BA1/2, BA3, PCC BA23, retrosplenial BA29 → Habenula, Posterior commissure or CM |
| 20260310 | 5 | 4bankeven | S1 BA1/2, BA3, PCC BA23, retrosplenial BA29 → Habenula, Posterior commissure or CM |
| 20260312 | 1 | Single column | S1 BA1/2, PCC BA23 → SC |
| 20260317 | 2 | Ch2H2sections3quarterD | S1 BA1/2 or PE BA5 → VIP → Inferior/Anterior pulvinar (Suprageniculate, MGN) |
| 20260318 | 0 | Ch2H0sections2 | S1 BA1/2, BA3 → VLP, VPI, ZI |
| 20260320 | 4 | Ch2H4quarterD | S1 BA1/2, BA3, (M1), PCC BA23 → Thalamus MD, CL, CM |
| 20260323 | 5 | Single column | S1 BA1/2, BA3, PCC BA23, retrosplenial BA29 → Habenula, Posterior commissure or CM |
| 20260324 | 2 | Ch2H2@53mm | S1 BA1/2 or PE BA5 → VIP → Inferior/Anterior pulvinar (Suprageniculate, MGN) |
| 20260325 | 1 | Ch2H1@53 | S1 BA1/2, PCC BA23 → SC |
| 20260326 | 2 | 4bankeven | S1 BA1/2 or PE BA5 → VIP → Inferior/Anterior pulvinar (Suprageniculate, MGN) |

**At planned full depth, every one of these chamber-2 trajectories is mostly or entirely subcortical** (thalamus, pulvinar, SC, habenula, ZI), passing through only a thin S1/PCC cortical rind near the top — relevant context for §6.1.

**Limitation (chamber 2): absolute per-channel depth-in-brain could not be reliably reconstructed.** The elab log has fields for "guide tube final depth," "probe tip enter guide tube depth," and "signal onset depth" (the last being the actual dura/brain-surface crossing — the reference point needed to convert electrode-Y-position-in-µm to true mm-below-cortex). "Signal onset depth" is blank for nearly every one of these 12 days, and only 2 (20260324, 20260325) have an explicit final recording depth at all (encoded directly in the channel-selection string, e.g. "Ch2H2@53mm"). The burr-hole/structure list above is therefore qualitative context (which region a session was *probably* in), not a per-channel anatomical registration — `sessionAnatomyInfo.m` documents this limitation in its header.

**Chamber 3 (10 sessions, prefrontal).** `Chamber3Holes.docx` (MRI-based planning for the second recording chamber) documents burr holes 1, 2, 4 and 5 — the numbers re-use chamber 2's, so a hole is always identified by chamber + number. Planned brain depth and cortical trajectory per hole: hole 1, 16 mm (OFC, VLPFC, DLPFC; lateral BA13, ventral BA46, dorsal BA46); hole 2, 19 mm (OFC, VLPFC, DLPFC, dorsal FEF; medial BA13, BA46, dorsal 8A); hole 4, 13 mm (OFC, VLPFC, DLPFC; medial BA12, ventral BA46, dorsal BA46); hole 5, 10.8 mm (?, VLPFC; medial BA12, (BA45a), ventral BA46). The current `Kwibus_elab.pdf` export (all sessions) gives, per day, the burr hole, online channel selection, probe final depth and signal-onset depth:

| Day | Burr hole | Channel selection (online) | Probe final depth (mm) | Signal onset depth (mm) | In-tissue depth used (mm) | Planned brain depth, hole (mm) |
|---|---|---|---|---|---|---|
| 20260819 | 1 | 4 bank | 46 | 32? | 14 | 16 |
| 20260821 | 5 | single column | 46 | 37–38 | 8.5 | 10.8 |
| 20260825 | 4 | single column | 46.5 | 38? | 8.5 | 13 |
| 20260903 | 2 | single column | 38 | 31? | 7 | 19 |
| 20260904 | 2 | 4bankeven | 46 | blank | unknown | 19 |
| 20260909 | 1 | single column | 48 | 32? | 16 | 16 |
| 20260910 | 4 | single column | 43 | 35 | 8 | 13 |
| 20260915 | 2 | Bank A | 35(.5?) | blank | unknown | 19 |
| 20260916 | 4 | single column | 44 | "don't know" | unknown | 13 |
| 20260918 | 2 | Bank A | 38.5 | 34.5? | 4 | 19 |

(`?` = the log entry itself carries a question mark. 20260916's entry reads "no obvious signal at any depth, not sure if it is the electrode or burr hole, stop early".) **In-tissue depth = probe final depth − signal onset depth**, which is what `sessionAnatomyInfo.m` uses as the estimated insertion depth for chamber 3 (§4.13). It agrees with the planned depth for 20260909 (16 vs 16) and roughly for 20260819 (14 vs 16), but is much shallower than planned in several sessions (20260903: 7 vs 19; 20260918: 4 vs 19), so — unlike chamber 2 — the planned depth is *not* used as a fallback here; where signal onset is blank (20260904, 20260915, 20260916), no absolute depth is assigned.

## 3. Pipeline architecture

**Figure organization**: `data/figures/` is organized **by figure type**, not by session — `figures/<type>/<Day>_run-<RunN>.png` (e.g. `figures/overview/20260312_run-003.png`), so the same plot across all 22 sessions is easy to find and compare side by side. `figSavePath.m` centralizes this convention for every per-session plotting function. Cross-session summary figures (no single day/run, e.g. the pooled decoding/area figures) live in `figures/pooled/` instead. The per-channel direction timecourses (§4.11) are organized session-first instead (`figures/timecourses/<Day>_run-<RunN>/{gocue,saclocked}/ch<NNN>.png`) — there are too many files (~6400) for the by-type convention to help there, and browsing one session's channels together is more useful for that figure type specifically.

| File | Purpose |
|---|---|
| `figSavePath.m` | Centralizes the `figures/<type>/<Day>_run-<RunN>.png` naming convention. |
| `findSaccadeRuns.m` | Scans `Data_collection/m032/*/run-*` for `stim_ct_classic_8targ.m`, returns a struct array of run metadata. |
| `matchRunRecording.m` | Maps a behavioral run to its OpenEphys `experimentN/recordingM` folder. |
| `matchTTLtoLogEvents.m` | Matches TTL digital-line numbers to behavioral event names by rising-edge count. |
| `extractSaccadeSession.m` | Core extraction: MUA + eye channels, trial table, eye calibration. One `.mat` output per session. |
| `runAllExtraction.m` | Batch driver over all 22 runs, tolerant of per-run failures. |
| `loadAllSessions.m` | Loads all extracted `.mat` files into one struct array; attaches `S.channelArea` (§4.12) to each. |
| `screenResponsiveChannels.m` | Per-channel, per-session responsiveness test in 4 task epochs vs. baseline. |
| `computeDirectionTuning.m` | 8-direction tuning curve, preferred direction, permutation significance, for responsive channels. |
| `analyzeEyeData.m` | Eye calibration application, saccade detection, fixation/pupil/main-sequence analysis. |
| `checkArtifactCorrelation.m` | Trial-residual pairwise channel correlation vs. distance, to distinguish shared-artifact from local-neural signal (§6). |
| `sessionAnatomyInfo.m` | Day → burr hole / channel-selection scheme / planned target structures / estimated insertion depth lookup (§2.2, §4.13). |
| `forceLightTheme.m` | Forces white background / black text on QC figures. |
| `classifyTuning.m` | Combines the vector-sum permutation test with a Kruskal-Wallis omnibus test into 4 tuning categories (§4.7): directional / directional_complex / complex / untuned. |
| `analyzePlanningExecution.m` | Splits go-cue-relative and saccade-onset-locked windows into planning vs execution epochs, separately screened and tuned (§4.8). |
| `summarizePlanningExecution.m` | Quantifies the planning/execution split across all 22 sessions (§6.6; weak sessions excluded from the pooled totals): recruitment fractions, preferred-direction consistency. |
| `checkSharedComponentKinematics.m` | Correlates the dominant shared trial-residual component (PC1) against saccade kinematics (§6.5), and common-mode removal for kinematic-linked sessions. |
| `decodeDirection.m` | 8-way direction decoding (individual channels + population-size sweep), all 3 epochs (§4.9). |
| `decodeTemporalGeneralization.m` | Cross-epoch (train/test) decoding matrix (§4.10). |
| `summarizeDecoding.m` | Pools decoding results across all 22 sessions. |
| `plotDirectionTimecourses.m` | Per-channel mean MUA timecourse by direction, go-cue-locked and saccade-onset-locked (§4.11). |
| `decodeAfterCommonModeRemoval.m` | Decoding before/after common-mode removal, for the 4 kinematic-linked sessions (§5.3). |
| `decodeSlidingGeneralization.m` | Continuous sliding-window temporal generalization matrix (§4.10, §5.3). |
| `generateAnatomyTemplate.m` | Builds the fillable per-channel/cluster anatomy assignment template (§4.12); `opts.chamber` writes a separate per-chamber template, and it refuses to overwrite existing files. |
| `getChannelArea.m` | Looks up assigned area/hierarchical labels and estimated depth-below-surface for a session's channels (§4.12, §4.13); reads the chamber-2 file plus every `anatomy/chamber<N>/` file. |
| `extractAreaTokens.m` | Parses candidate area abbreviation tokens out of a free-text anatomy label, for liberal multi-label pooling (§6.7). |
| `poolTuningByArea.m` | Pools responsiveness/tuning across all 22 sessions by area, liberally (§6.7). |
| `decodeByArea.m` | Direction decoding restricted to each area's channels, size-matched across sessions (§6.8). |
| `poolPlanExecByArea.m` | Planning-vs-execution recruitment (§6.6) pooled by area instead of by session. |
| `runFullPipeline.m` | Runs all of the above end to end. |

## 4. Key methodological decisions (and why)

### 4.1 TTL bit codes are matched empirically, not from `par_default.m`

`par_default.m` defines `Par.StimB=1, Par.TargetB=2, Par.RewardB=3, Par.MicroB=6, Par.CorrectB=7`, with a comment that `SaccadeB=4` and `TrialB=5` are handled by a separate "DasControl" system. **None of this matches what's actually on the wire.** The NI-DAQ digital port carries exactly lines `{1, 4, 6, 8}` in every one of the 22 sessions, no exceptions, and cross-checking rising-edge counts against `Log.events` type counts shows this assignment is constant across the whole dataset:

- line 4 = `trial_start`
- line 6 = `targ_start` (the go-cue / saccade instruction)
- line 8 = `correct` / `reward_start` / `reward_stop` (these three are logged the same number of times whenever every trial that reaches the target epoch is rewarded, so they can't be told apart by count alone — harmless, since the pipeline only needs `targ_start` and `trial_start`)
- line 1 = fires more often than any single logged event type; likely a per-attempt fixation/stim marker with no 1:1 logged counterpart.

`matchTTLtoLogEvents.m` does this count-matching **per run** rather than hard-coding line numbers, as a robustness measure. On two sessions (20260317, 20260318), `targ_start`'s count ties exactly with the `correct`/`reward` cluster (all logged the same number of times that day, since every trial that reached the target was correct), making the match ambiguous by count alone. This is resolved with a secondary heuristic in `extractSaccadeSession.m`: among tied candidate lines, pick whichever has the smaller median latency after the (unambiguous) `trial_start` line — `targ_start` reliably follows `trial_start` by ~0.6s, while `correct`/`reward` follow much later, after the saccade.

### 4.2 Run-to-recording matching

`matchRunRecording.m` matches a run to a recording by **closest absolute time** between the run's folder timestamp and the recording's `sync_messages.txt` "Software Time," within a 90s tolerance — not "closest *preceding*," because operator workflow was not perfectly consistent: usually the OE recording starts a couple of seconds *before* the run script, but on at least one day the recording started *after* the run folder was created instead (20260303/run-003's matching recording, `experiment2/recording3`, starts 6.3s after its run folder timestamp).

Each behavioral run maps to its own dedicated recording — verified across every run on all 22 days, no recording is ever shared between two runs. `extractSaccadeSession.m` still time-windows the TTL stream to the run's own approximate span (its start time plus its logged duration from `Log.events`, with generous ±20/+30s margins) before count-matching, as a defensive precaution. If that run-time window clips trials (the stim-PC and OpenEphys clocks differ by ~35s on 20260819, 20260825 and 20260915), extraction falls back to the full-recording TTL stream.

### 4.3 Eye channels: empirical identification, not documented anywhere

There is no config file, log field, or channel-name metadata anywhere in the rig, log `Par` struct, or `structure.oebin` that says which of the 8 NI-DAQ analog channels (AI0–AI7) is eye X, eye Y, or pupil. This was inferred from the raw signal shape:

- **AI0, AI1**: clean, correlated step-changes of ~1–4V, time-locked to the go-cue, returning to baseline ~300–450ms later — the signature of a saccade-and-return. Assumed **Eye X, Eye Y**.
- **AI2, AI3**: sit near the ADC negative rail (−10 to −5V), with the raw value pinned at exactly −10.000V for extended stretches (classic blink artifact where the pupil tracker loses the pupil). Assumed **pupil-related**; AI2 is used as "pupil" in the pipeline, AI3 is unused.
- **AI4–AI7**: near-zero, small (~0.05V) crosstalk-looking modulation. Assumed unused/floating.

This is a data-driven inference, not a documented fact (flagged in every extracted session as `S.eyeChannelAssumption`), but well validated indirectly: the eye-position calibration below achieves R² of 0.84–0.97 against the known 8 target locations in 21 of 22 sessions (20260819 is weaker, 0.72/0.60 — see §5.1), which would not happen if AI0/AI1 weren't genuinely eye position.

### 4.4 Eye calibration and saccade detection

`Par.SCx/SCy/OFFx/OFFy` in the logs are a live GUI display scale, not a true volts→degrees calibration (identical to `par_default.m`'s hardcoded defaults in every session, never updated by a calibration routine). Instead, calibration is derived directly from task structure: since the target's position is known exactly on every trial, and on **correct** trials the animal's gaze is by definition on the target, the mean raw voltage during a post-saccade landing window (200–350ms post-go-cue — chosen because this task's dwell on target is brief, with the return saccade beginning by ~350–400ms) is regressed against the known target (x,y) using a full 2-D affine fit:

```
[dva_x, dva_y] = [V_x, V_y, 1] · M
```

A full affine (rather than two independent 1-D linear fits) is necessary because the raw channels show real cross-talk — a pure-vertical target measurably displaces the "X" channel and vice versa, consistent with a rotated/skewed camera axis in the analog eye tracker.

**Blink handling.** Blinks corrupt the raw eye **position** channels (X/Y), not just the pupil channel: during a blink, X/Y pin near a fixed hardware rail of −5.00V (identical across the original 12 sessions, regardless of each session's own raw baseline, which itself varies from about −2.3V to +0.8V; the same criterion is applied unchanged to the chamber-3 sessions) for roughly 65–180ms at a time. `extractSaccadeSession.m` detects this by proximity to the known rail voltage (`abs(raw−(−5.00))<0.5V`) and excludes these samples from the calibration fit; `analyzeEyeData.m` interpolates over them before the zero-phase velocity filter (which cannot tolerate NaNs) and then re-masks them (±10ms padding) before anything is plotted or used for saccade onset/kinematics detection.

Some **moderate** (3–7dva) landing scatter remains on plenty of correct trials — that's expected, not an error: `Stm.TargWinSz` = 7dva, i.e. the task's own online correct/error judgment accepts any landing within 7dva of the target center, and real saccades commonly under/overshoot by a few degrees before a small corrective adjustment. The trajectory plots (`*_eye_summary.png`) draw a light-grey dotted circle of this radius around each target so this is visually obvious.

**Saccade detection and kinematics** (`analyzeEyeData.m`). The calibrated, blink-interpolated position trace (X/Y, see above) is first low-pass filtered (2nd-order Butterworth, zero-phase `filtfilt`, 50Hz cutoff) — raw analog position is noise-dominated sample-to-sample (30kHz ADC noise gives spurious velocities of 10⁵–10⁶ dva/s if differentiated directly), so filtering before differentiating is standard practice for saccade detection, not an optional smoothing step. Velocity is the finite difference of the filtered X/Y trace (`hypot(diff(X),diff(Y))·Fs`, dva/s). Within a search window of [-0.05, +1.0]s relative to the go-cue, saccade **onset** is the first sample where velocity exceeds **60 dva/s**; **offset** is the first subsequent sample where it drops back below that threshold. A candidate is only accepted as a real saccade if duration ≥8ms AND amplitude ≥2dva (endpoint minus onset position) — both are noise-rejection floors, not physiological claims: velocity-threshold crossings shorter/smaller than this are far more consistent with filter ringing or micro-saccades than the large targeted saccades this task calls for. Peak velocity is the maximum filtered velocity between onset and offset (the raw "main sequence" measure in `*_eye_summary.png`). Blink-masked samples (position channels, §4.4 above) are excluded before this whole pipeline runs, propagating as `NaN` through the velocity trace and correctly preventing a blink from being detected as a spurious saccade.

The **pupil** channel gets its own, separate and simpler blink criterion — |raw voltage| > 9.5V (near the ADC's −10V negative rail) is masked to `NaN` — distinct from the position channels' rail-proximity criterion above (−5.00V ± 0.5V), because the pupil channel saturates at a different, more extreme rail value during a blink/dropout than the position channels do. The two masks are computed and applied independently; a sample can be blink-masked in position without necessarily also being masked in pupil, and vice versa.

### 4.5 MUA extraction

Reuses the filter chain from `XinyuScripts/getmua.m` (the lab's existing pipeline for this same probe/task family), all filtering **zero-phase** (`filtfilt`, 3rd-order Butterworth throughout — no time-shift artifacts introduced relative to the original AP-band data, which matters for anything downstream that reasons about precise latency, e.g. §4.8's epoch splits or §4.11's timecourses):

1. **300Hz high-pass** — removes LFP/slow drift, leaving the spiking-band signal.
2. **12-way ADC-multiplexing phase correction** — this probe's ADC samples its 12 channel-groups in an interleaved (not simultaneous) sequence; uncorrected, this shows up as a small per-channel time offset that would otherwise misalign fast transients across channels.
3. **5000Hz low-pass** — anti-alias/band-limit ahead of decimation.
4. **Spatial "destripe" high-pass across channels** — a 3rd-order Butterworth high-pass applied along the **channel axis** (not time), cutoff equivalent to a ~100-channel spatial period. This removes slowly-varying-with-depth trends common to many nearby channels at once (the spatial analogue of a common-average reference — e.g. a broadband electrical/reference offset shared by a wide swath of the probe), while preserving genuine channel-to-channel variation over shorter spans. This runs on every session as a standard part of extraction, **before** responsiveness/tuning is ever computed — it is a different, earlier step from §6.4/§6.5's shared-component check, which looks for a *residual* common-mode signal in **trial-to-trial** fluctuations (not this broadband spatial baseline) that survives destriping; the two are complementary, not redundant, and §6.4's finding (a small nonzero long-distance correlation plateau in every session) is exactly the kind of thing this filter would not be expected to catch, since it targets spatial structure, not trial-to-trial covariance.
5. **Rectify (`abs`) + 200Hz low-pass envelope** — the actual MUA "amplitude" signal.
6. **Decimate to ~1kHz.**

Channels are reordered by depth using `NP_PROBE.ELECTRODE_YPOS` from `settings.xml`. Trials are windowed 1.0s before to 1.2s after the go-cue (`targ_start`), covering the full fixation→stimulus→delay→saccade→return sequence.

### 4.6 Responsiveness screening

Per channel, per session: baseline-subtracted mean MUA in 4 windows (visual [-0.5, -0.35s], delay [-0.3, -0.05s], peri-saccadic [0, 0.25s], post-saccadic [0.3, 0.6s]) vs. a pre-stimulus baseline [−0.9, −0.6s], tested with a Wilcoxon signed-rank test across correct trials, Benjamini-Hochberg FDR-corrected across channels per window, plus a minimum effect-size requirement (|mean diff| > 0.3 baseline SD) to avoid flagging trivially small but "significant" effects given large trial counts. A channel is "responsive" if any window passes both criteria.

### 4.7 Direction tuning

For responsive channels only: mean baseline-subtracted response in the peri-saccadic window [0, 0.3s] relative to go-cue (see §4.8 for a version split into planning/execution epochs), grouped by the 8 target directions (only correct trials). Preferred direction and tuning strength are computed via **vector sum**: each direction's trial-averaged response is rectified (negative responses clamped to 0, since a negative "weight" has no sensible meaning as a vote for that direction) and used as the weight for that direction's unit vector; the resultant vector's angle is the preferred direction, its length (`resultantLen`, 0-1) is the tuning strength. Significance is tested by **permutation** (1000 shuffles of direction-group labels across trials, recomputing the resultant length each time to build a null distribution of what resultant length would look like from these same responses with no true direction structure), with the p-value computed as `(1 + #(null ≥ observed)) / (1000 + 1)` — the standard "+1" correction that guarantees a nonzero p-value regardless of how extreme the observed statistic is relative to the (finite) null sample, rather than the naive `#(null ≥ observed)/1000` which can spuriously report exactly p=0.

The vector-sum/permutation test is only **guaranteed** to fail for a perfectly symmetric bimodal response — two opposite directions, exactly equal magnitude. Any asymmetry (unequal peak heights, or peaks not exactly 180° apart) leaves a nonzero residual vector, which can easily be large enough to pass the permutation test even though the curve visibly has two lobes, not one — so "vector-sum significant" does **not**, by itself, mean "single clean peak." `classifyTuning.m` addresses this by also running a **Kruskal-Wallis omnibus test** across the 8 direction groups (does response differ across directions *at all*, regardless of shape), and — since 2026-07-14 — crossing the two tests directly (2×2) into **four** categories, used for all tuning figures, rather than letting a significant vector-sum test silently override the omnibus result:

- **`directional`** (red): vector-sum significant, Kruskal-Wallis **not** — a genuine, clean single preferred direction; `prefDir`/`resultantLen` are reliable summaries.
- **`directional_complex`** (purple): **both** significant — a real dominant peak exists (`prefDir`/`resultantLen` still computed and usable as an estimate of the *stronger* of at least two real preferences), but real secondary structure coexists with it; don't treat `prefDir` as a clean single-peak summary for these — look at `meanByDir` to see the actual shape.
- **`complex`** (orange): Kruskal-Wallis significant, vector-sum **not** — real direction-dependent structure exists with no reliable single dominant peak; `prefDir`/`resultantLen` are **not** meaningful for these, look at the actual per-direction response instead.
- **`untuned`** (gray): neither test significant.

**Why this split matters, checked directly against the data (not assumed): among channels that pass the vector-sum test, having real secondary structure too is the common case, not a rare edge case.** Before this split, every vector-sum-significant channel was labeled simply `directional` regardless of whether it also had a significant Kruskal-Wallis result — silently discarding that information. Checked against the full `tuning_results.mat`: of the 1486 channels (22 sessions) that pass the vector-sum test, **1003 (67.5%) also pass Kruskal-Wallis**; in the original 12 sessions, where this was also checked for shape, of 832 vector-sum-significant channels 576 (69%) passed Kruskal-Wallis and **411 (49%) had a secondary local peak at least half the height of the primary one, at least 90° away** — a genuinely two-lobed shape. A clean concrete example: channel 152 (20260324) has `meanByDir` = [0.63, 1.21, 1.13, **1.70**, 1.06, **−0.21**, 0.18, **1.39**] (0°→315° in 45° steps) — clear peaks at 135° and 315° with a real dip at 225° between them — and was previously labeled plain `directional` (pPerm=0.001) with no indication that pKW was *also* significant (<0.0001). It is now `directional_complex`.

**Existing counts/figures that summed over the old single `directional` category are unaffected in total**, since `directional` (old) = `directional` (new) ∪ `directional_complex` (new) by construction of the 2×2 split — any pooled "%directional" number elsewhere in this report reflects that same union unless stated otherwise. What changes is that the figures show a 4th color, and it is possible to ask (§5.1, §6.7) what fraction of "directional" channels are cleanly unimodal vs. have real secondary structure.

**Why some visually-peaked-looking curves in `tuning_examples.png` aren't classified as tuned.** With ~8-18 trials per direction and single-channel MUA noise, this is a mix of two things, not just one: (1) the figure itself, until now, min-max normalized every curve to its own [0,1] range with no error bars — this stretches any wobble, real or pure noise, to fill the same visual range, so a small noisy peak looked identical to a large reliable one; (2) genuinely borderline cases are common at this trial count. Checked directly: channel 269 (20260323) has a clean-looking peak (normalized 1.0 at 315°, next-highest 0.29) that remains visually convincing even WITH SEM error bars added (peak's lower bound clears every other direction's upper bound) — yet `pPerm` ≈ 0.055-0.061 across repeat runs, just missing the p<0.05 cutoff. This is not hidden noise or an overly strict threshold; it is a real value sitting right at the significance boundary, which is expected to happen for *some* channels given how significance testing works, not a sign of miscalibration. `tuning_examples.png` now plots SEM error bars (scaled by the same per-channel normalization factor as the curve itself) specifically so this distinction — "a visually strong peak with tight error bars" vs. "a visually strong peak whose error bars overlap the rest" — is visible directly in the figure, rather than requiring a trip to `tuning_results.mat`. The same fix applies to the `*_planexec_<mode>_curves.png` figures (§4.8) — `analyzePlanningExecution.m`'s `fitTuning` now also computes `semByDir`, plotted as error bars in raw (non-normalized) response units there.

### 4.8 Dissociating planning (premotor) from execution (perimovement) tuning

The [0, 0.3s] go-cue-relative window above mixes two different things whenever saccade latency is non-negligible (median 68–237ms across sessions, §5.1): **planning/premotor** activity (from the go-cue instruction until the eye actually starts moving) and **execution/perimovement** activity (during and immediately after the saccade itself). `analyzePlanningExecution.m` splits this two ways, producing figures separate from the main pipeline outputs:

- **`gocue` mode** — fixed split relative to the go-cue, same cutoff for every trial: planning=[0, 0.1]s, execution=[0.1, 0.3]s.
- **`saclocked` mode** (the methodologically preferred one, since latency varies so much across sessions) — anchored to each **trial's own detected saccade onset**: planning=[-0.15, -0.02]s and execution=[-0.02, +0.15]s relative to onset. Trials with no detected saccade (0–4% depending on session) are excluded from this mode.

For each mode: responsiveness (matching §4.6) and direction tuning (§4.7's three-category scheme) are computed **separately** for the planning and execution windows. Outputs per session, both modes: `*_planexec_<mode>.png` (per-channel planning-vs-execution effect map, tuning-strength scatter, preferred-direction scatter) and `*_planexec_<mode>_curves.png` (every responsive channel's actual tuning curve, planning on top/execution below, category-colored). Results quantified across all 22 sessions are in §6.6.

### 4.9 Decoding: can the recorded population actually tell the 8 directions apart?

Stimulus position and saccade direction are the same label in this task — the target position on a given trial is exactly the instructed saccade direction — so this is one 8-way classification question, asked separately of the same three epochs as §4.8: decoding from the **visual** epoch asks whether stimulus position is represented; decoding from **planning**/**execution** asks the same of the upcoming/ongoing saccade. `decodeDirection.m` implements this as a direct complement to the tuning-based analyses above (which test each channel individually against a parametric null), giving a model-free readout of how much direction information is actually present, and how it scales with the number of channels pooled together.

**Decoder**: nearest-centroid classification on per-channel z-scored responses, evaluated with leave-one-trial-out cross-validation (train on all trials but one, classify the held-out trial by its Euclidean distance to each of the 8 per-direction training centroids, repeat for every trial). Chosen over LDA/SVM because population sizes here range up to ~380 channels while trial counts are only ~85–145 per session — far more features than trials, a regime where a full-covariance classifier is singular; per-channel standardization plus nearest-centroid is the standard robust choice here (equivalent to a diagonal/naive-Bayes LDA) and needs no hyperparameters. Chance level for 8 balanced classes is 1/8 = 12.5%. **Critically, both the per-channel standardization (mean/SD) and the 8 class centroids are recomputed from the training fold only, excluding the held-out trial, for every single fold** — the held-out trial never contributes to its own classification statistics, so there is no train/test leakage inflating the reported accuracy (implemented with closed-form leave-one-out updates to the training-fold sums, not by literally refitting from scratch each fold, purely for speed — mathematically identical to the naive refit).

Two complementary views, per session per epoch:
- **Individual channels**: every channel decoded alone (not just responsive ones, so unresponsive channels serve as an internal near-chance control).
- **Population**: accuracy vs. population size, built by randomly subsampling 1, 2, 4, ..., up to the full responsive population (20 repeats per size, averaged), the standard "how many channels before performance saturates" curve.

`summarizeDecoding.m` pools this across all 22 sessions (restricted to the common population sizes shared across sessions, since each session's own full-responsive-population size is idiosyncratic and not itself poolable).

### 4.10 Temporal generalization: is it the same code in each epoch, or a different one?

Direct follow-up to §4.9: decoding well in an epoch only shows that direction information is present *at that time* — it doesn't say whether the population structure that carries it (which channels prefer which direction) is the *same* structure across epochs, or gets reconfigured. `decodeTemporalGeneralization.m` tests this directly: train the nearest-centroid decoder on one epoch's population response, test it on **another** epoch's response for the same held-out trials (leave-one-trial-out, full responsive population), for every (train, test) pair among visual/planning/execution — a 3×3 matrix per session. The diagonal (train and test the same epoch) is exactly §4.9's ordinary same-epoch accuracy; the off-diagonal cells are the new information: if a decoder trained on one epoch still classifies well on another, the code generalizes (shared structure); if it drops to chance, the code is epoch-specific (reconfigured). Standardization (z-scoring) is always computed from the **train** epoch's own statistics and applied to the test epoch's data, so this is a genuine test of whether the trained decision boundary transfers, not just whether both epochs are independently informative.

This first pass uses the same 3 fixed epoch windows as §4.8/§4.9 (a coarse 3×3 matrix). `decodeSlidingGeneralization.m` extends this to a **continuous** version: 50ms bins spanning [-0.5, +0.8]s relative to go-cue (26×26 bins), same method, full responsive population, one matrix per session. Results in §5.3.

### 4.11 Per-channel direction timecourses

A direct, model-free complement to the tuning/decoding analyses above: for every responsive channel, the mean MUA timecourse (correct trials only) split into the 8 directions and overlaid, rather than reduced to a single-window number. Two alignments, `plotDirectionTimecourses.m`:

- **go-cue-locked**: the trial's own time base, window [-0.9, +1.0]s relative to go-cue — wide enough to show fixation, stimulus flash, delay, go-cue, target-end, and the return saccade all in one plot.
- **saccade-onset-locked**: each trial's MUA linearly interpolated onto a common time axis relative to its own detected saccade onset (necessary because onset time varies trial to trial — a simple index shift would misalign different trials' actual event times), window [-0.5, +0.8]s relative to onset.

Each direction's mean trace is smoothed (zero-phase Gaussian, 50ms window) before plotting — purely a display transform on the already-trial-averaged trace, not a change to `S.MUA` or any other analysis in this pipeline. A lighter 20ms window was tried first and left the trace visually indistinguishable from unsmoothed, since the underlying MUA is itself a 200Hz-bandwidth envelope with real high-frequency content; 50ms visibly reduces trial-count noise while still preserving genuine trial-to-trial structure (checked against an unsmoothed reference before choosing this value).

This produces one figure per responsive channel per session per alignment (5467 responsive channels × 2 alignments = 10,934 figures total), organized into subfolders rather than the flat `figures/` directory used elsewhere in this pipeline:

```
data/figures/timecourses/<Day>_run-<RunN>/gocue/ch<NNN>.png
data/figures/timecourses/<Day>_run-<RunN>/saclocked/ch<NNN>.png
```

`<NNN>` is the channel index (1 = deepest / probe tip, matching every other figure in this pipeline), not a physical electrode ID. Given the volume, these are **not** embedded in the PDF appendix (§9) — browse them directly on the server if you want to inspect individual channels' raw timecourses rather than the tuning-curve summaries.

### 4.12 Anatomy assignment template

Per-channel anatomical area is not currently recoverable from the recording logs (§2.2's limitation). `generateAnatomyTemplate.m` (new) does not solve this itself — it builds a fillable template so you can assign areas by hand, informed by this pipeline's own recorded activity rather than starting from a blank per-channel list:

- Real depth (`S.CHdepthUm`), responsiveness, tuning category, preferred direction, and tuning strength for every channel in every session (`anatomy/anatomy_template_channels.csv`).
- Channels are grouped into **candidate clusters** (`anatomy/anatomy_template_clusters.csv`, the document meant to be filled in): contiguous channel ranges whose recorded activity looks similar, detected by binning channels (20 at a time) and greedily merging adjacent bins whose responsiveness rate and resultant-length-weighted circular mean preferred direction are both close. This is deliberately bin-then-merge rather than a pure jump detector: at least one session (20260312) shows preferred direction drifting smoothly and continuously with depth, with no sharp transition anywhere (checked directly) — a jump detector would have collapsed that whole session into one uninformative cluster.
- Each cluster row also carries the planned burr-hole/structure list from `sessionAnatomyInfo.m` for context, and three blank columns to fill in by hand: `AssignedArea`, `Confidence`, `Notes`.
- This is a heuristic proposal for where to look, not a certified segmentation — split, merge, or ignore clusters freely against your own reading of the atlas and elab notes.

**Filled in and wired into the pipeline (chamber 2).** Every cluster in `anatomy_template_clusters.csv` now has a hand-assigned `AssignedArea`, cross-referenced against `anatomy_template_clusters_hierarchical.csv` (a parallel hierarchical parcellation, L6=finest to L3=coarsest, per cluster). Both are propagated down to every individual channel in `anatomy_template_channels.csv` (`ClusterID`, `AssignedArea`, `L6`-`L3` columns, all 4608 channel rows across the 12 chamber-2 sessions, verified zero gaps/overlaps in cluster coverage).

`getChannelArea.m` (new) looks up these labels for a given day/run/channel(s) directly, and `loadAllSessions.m` now attaches this automatically as `S.channelArea` for every session — `S.channelArea.AssignedArea{c}` (finest tentative label) or `S.channelArea.L6/.L5/.L4/.L3{c}` (hierarchical, coarser toward L3), for channel `c` (1 = deepest, matching every other per-channel array in this pipeline, e.g. `S.CHdepthUm`). These are **tentative, hand-assigned labels compared against an atlas by eye, propagated per-cluster** — not a verified histological registration — but now available to any analysis in this pipeline without re-deriving them.

**`L6` code fix.** `L6` was originally just a verbatim copy of `AssignedArea`'s free text, with no parenthetical abbreviation code (e.g. "Habenula" rather than "Habenula (Hb)") — `L5`/`L4`/`L3` all carry the `(Code)` convention, `L6` didn't. This meant area-pooling analyses keying on `L6`'s code found essentially nothing (4.3% coverage vs. ~92% at `L5`-`L3`). Since `L6` and `L5` turned out to be near-1:1 in content (same structure, `L5` just phrased with canonical names), each of the 43 unique `L6` strings was given the matching code from its `L5` counterpart, then re-propagated into `anatomy_template_channels.csv` — `L6` now carries proper codes at ~92% coverage, matching `L5`-`L3`, while keeping its own (sometimes more colloquial, e.g. "PCC [BA23C]") phrasing.

**Chamber 3 (prefrontal, 10 sessions, 20260819–20260918).** Same approach, kept in **separate files** so the hand-edited chamber-2 files can never be overwritten, and in the same layout so the two can be merged later: `anatomy/chamber3/anatomy_template_clusters.csv` (93 clusters) and `anatomy/chamber3/anatomy_template_channels.csv` (3840 rows = 10 × 384, with a `Chamber` column). `generateAnatomyTemplate.m` gained `opts.chamber` (selects sessions by `sessionAnatomyInfo(day).chamber`, writes into `anatomy/chamber<N>/`, and adds blank `L6`–`L3` columns to the cluster file so it has the same columns as the chamber-2 hierarchical file) and now **refuses to overwrite** any existing output unless `opts.overwrite=true`; its directional counts use `directional` ∪ `directional_complex` (§4.7). Clusters per session (bin-then-merge, same parameters as chamber 2): 0819 9, 0821 9, 0825 8, 0903 12, 0904 15, 0909 9, 0910 4, 0915 10, 0916 4, 0918 13. One auto-detected 0916 cluster (channels 41–384, 0.8–7.66 mm) was split by hand at 5.9 mm into channels 41–295 and 296–384 (statistics recomputed for both halves); that session therefore has 4 clusters.

*Labels.* `AssignedArea` was hand-assigned for all 93 clusters against the atlas, using six CHARM-named cortical areas: `dorsal_periarcuate_area_8A` (36 clusters), `dorsal_area_46` (21 clusters), `ventral_area_46` (20 clusters), `area_46_in_fundus_of_the_principal_sulcus` (9 clusters), `lateral_area_13` (4 clusters), `medial_area_12` (3 clusters). `L6`–`L3` were generated mechanically from the parent chain in `CHARM_key_all.txt` (the key lists parcels in depth-first order with the levels at which each parcel exists; a parcel's ancestor at a coarser level is the most recent preceding row that spans that level), each written as `full_name (code)` — official codes from the start, avoiding the chamber-2 `L6` code problem above:

| `AssignedArea` (= CHARM L6 name) | L6 | L5 | L4 | L3 |
|---|---|---|---|---|
| dorsal_periarcuate_area_8A | area_8Ad | area_8A | area_8A | area_8A |
| dorsal_area_46 | area_46d | area_46d | rostral_dlPFC | dlPFC |
| ventral_area_46 | area_46v | area_46v/f | rostral_vlPFC | vlPFC |
| area_46_in_fundus_of_the_principal_sulcus | area_46f | area_46v/f | rostral_vlPFC | vlPFC |
| lateral_area_13 | area_13l | area_13 | areas_11/13 | lat_OFC |
| medial_area_12 | area_12m | area_12m/o | area_12m/o | lat_OFC |

(The key's full name for `area_8A` carries its own parenthetical "(Frontal_Eye_Fields)"; it is written as `periarcuate_area_8A (area_8A)` so the area-token extractor does not read it as a second code.) These labels were then propagated to every channel by channel-index range: all 3840 rows matched, no gaps or overlaps, and the per-area channel totals equal the cluster table's `NChan` sums.

*White matter.* There is no white-matter label: some channels sit in or near white matter but not in clearly separable clusters, so they were assigned to the nearest cortical area. 21 clusters (0819: 6, 0821: 5, 0903: 5, 0904: 4, 0918: 1) carry "might contain WM" in `Notes` (these notes live in the cluster table only; join on `Day`/`RunN`/`ClusterID` for the channels). The consequence is that these channels are pooled as cortical in §6.6–6.8, which can only dilute an area's % responsive/% tuned slightly; this has not been quantified.

*Code extraction.* Several official codes contain `/` (`area_46v/f`, `areas_11/13`, `area_12m/o`) or exceed 12 characters (`rostral_dlPFC`, `rostral_vlPFC`), which `extractAreaTokens.m`'s code pattern silently dropped; the pattern now allows `/` and up to 15 characters. Re-running the chamber-2 pooling after the change reproduces the earlier tuning/decoding numbers exactly, adding only one new chamber-2 token (`AREA_3A/B`, previously never extracted).

### 4.13 Estimated absolute depth below the brain surface

§2.2 flagged that absolute per-channel depth-in-brain (mm below cortical/dura surface, as opposed to `S.CHdepthUm`'s relative depth along the probe) wasn't reconstructable from the logs. Partial fix (chamber 2): for the 2 sessions where the elab channel-selection string encodes an explicit final depth (20260324 "Ch2H2@53mm", 20260325 "Ch2H1@53" → 53mm), that value tracks the burr hole's Chamber2Holes **planned** depth closely — so for the other 10 sessions (no explicit number logged), the planned depth is used as the estimated insertion depth instead. `sessionAnatomyInfo.m` now returns `estimatedInsertionDepthMm` (53mm for 20260324/20260325, else that day's burr-hole planned depth) and `depthSource` ('explicit' or 'planned').

`getChannelArea.m` uses this to compute, per channel: `DepthBelowSurfaceMm = estimatedInsertionDepthMm - CHdepthUm/1000` (channel 1, the tip, gets the full insertion depth; more distal channels get progressively less). A **negative** value means that channel is estimated to sit *above* the brain surface that day — in the guide tube, dura, or CSF, not brain tissue. Checked across the 12 chamber-2 sessions: several show a chunk of shallow channels estimated above surface (e.g. 20260224: 49/384, 20260226: 31/384, 20260317: 7/384), others show none (e.g. 20260312, 20260323, 20260324, 20260325, 20260326) — plausible given each day's online channel-selection scheme spans a different fraction of the probe relative to how far it was actually driven in.

This is still an **estimate**, not a verified measurement — it inherits whatever error exists in the Chamber2Holes plan vs. actual insertion for the 10 sessions using the fallback, and the underlying burr-hole targeting itself. Available as `S.channelArea.DepthBelowSurfaceMm(c)`, `.EstimatedInsertionDepthMm`, `.DepthSource` for every channel via `loadAllSessions.m`.

**Chamber 3.** The planned-depth fallback does not work there (§2.2: in-tissue depth is often far shallower than the Chamber3Holes plan), so `sessionAnatomyInfo.m` instead uses the **elab-logged in-tissue depth = probe final depth − signal onset depth** (`depthSource = 'elab_final_minus_onset'`, with an `onsetUncertain` flag where the log marks the onset with `?`: 20260819, 20260825, 20260903, 20260909, 20260918) and returns NaN (`depthSource = 'unknown'`) where signal onset is blank or unknown (20260904, 20260915, 20260916), so no absolute depth is assigned for those three. Estimated insertion depths (mm): 20260819 14, 20260821 8.5, 20260825 8.5, 20260903 7, 20260909 16, 20260910 8, 20260918 4. Only two chamber-3 sessions have channels estimated above the brain surface (20260819: 34/384, 20260903: 33/384); the others have none. These depths inherit the uncertainty of the logged onset (often a `?` entry, i.e. a visual/audio judgment of where signal begins) and are best treated as ±several mm; the planned hole depth is still returned separately (`plannedDepthMm`) for comparison.

### 4.14 Area-pooling methodology (§6.6-6.8)

Three analyses (`poolTuningByArea.m`, `decodeByArea.m`, `poolPlanExecByArea.m`) share a common method for turning per-channel anatomical labels (§4.12) into cross-session, by-area statistics. Documented here once rather than separately per figure.

**Token extraction (`extractAreaTokens.m`).** A channel's area membership is read from the parenthetical abbreviation code(s) in its label string (e.g. "(SC)", "(Hb)", "(fundus_IPS)"), upper-cased so inconsistent capitalization of the same area across hand-typed entries collapses together (this dataset has "(MThal)"/"(Mthal)", "(GThal)"/"(Gthal)", etc. for the same area). A descriptive full name is also captured, from the ≤3-word phrase immediately preceding each code, for human-readable reporting only (not used for pooling membership). **Code pattern (updated 2026-10-23):** the code character class allows letters, digits, `_`, `-` and `/`, up to 15 characters, so official CHARM/SARM codes such as `area_3a/b`, `area_46v/f` or `rostral_vlPFC` are extracted; previously `/` and codes over 12 characters were silently dropped.

**Liberal, single-level pooling.** Many hand-typed labels genuinely name more than one candidate area (e.g. "Posterior Commisure (or Habenula or CM?)") — a channel counts toward **every** token found in its label, not just a primary guess, so it can and does appear in more than one area's pool. This is deliberate (treat as evidence, not a forced single identity), but must operate on **exactly one hierarchy level's label at a time** (`opts.labelField`, one of `AssignedArea`/`L6`/`L5`/`L4`/`L3`) — concatenating multiple levels together before extracting tokens would conflate genuine within-level candidate uncertainty with cross-level hierarchy (the same channels appearing as two "different" areas just because they're labeled "SC" at one granularity and "colliculi (Co)" at a coarser one). All three analyses default to `L6` (finest) and are also run at `L3` (coarsest) to compare grains (§6.7-6.8 report both). Two tokens are dropped everywhere (`opts.excludeTokens`, default `{'DEEP','HIGH'}`) — these turned out to be position qualifiers within a compound cluster label (e.g. "(deep)"/"(high up)" distinguishing the two ends of a probe span assigned to different areas), not area names.

**Per-analysis specifics:**
- `poolTuningByArea.m` — pools §4.6/§4.7's responsiveness and tuning-category counts across every session, per area token: %responsive/%directional/%complex, plus a resultant-length-weighted circular mean preferred direction (for continuity with the summary table only — see §6.7 for why the preferred-direction *distribution*, not this collapsed mean, is the actual reported result). Areas with fewer than `opts.minChannels` (default 20) pooled channels are dropped from the summary/rose figures (not from the underlying `.mat`).
- `decodeByArea.m` — restricts each session's *own* decoding population (§4.9 method, same epochs) to that area's channels in that session, **size-matched at `opts.capSize` (default 20) channels** for a fair cross-area comparison — sessions with more than 20 area-X channels are randomly subsampled (`opts.nReps`=20 repeats, averaged), sessions with fewer use all available. Trial-level data can't be pooled across sessions the way per-channel tuning summaries can, so each qualifying session contributes one accuracy number per area per epoch; an area is only reported if it reaches `opts.minChannelsPerSession` (default 15) in at least `opts.minSessions` (default 2) sessions.
- `poolPlanExecByArea.m` — pools §4.8's saccade-onset-locked (`opts.mode`, default `'saclocked'`) planning/execution tuning-category counts by area, giving %execution-only and %planning-only recruitment fractions per area (§6.6's pooled numbers, broken out). The four weak sessions are excluded, as in §6.6 (`opts.excludeDays`, default `{20260224, 20260326, 20260910, 20260916}`). An area is only reported once it has at least `opts.minExecDirectional` (default 8) execution-directional channel-instances pooled, to avoid noisy small-N percentages.

## 5. Results

### 5.1 Per-session summary

| Day | Trials (correct) | Responsive | Tuned | Eye calib R² (x/y) | Saccade detected | Median latency |
|---|---|---|---|---|---|---|
| 20260224 | 85 (85) | 64/384 | 4 | 0.84 / 0.84 | 84/85 | 200ms |
| 20260226 | 145 (144) | 271/384 | 47 | 0.93 / 0.95 | 144/145 | 195ms |
| 20260303 | 114 (112) | 201/384 | 15 | 0.97 / 0.96 | 113/114 | 179ms |
| 20260310 | 117 (116) | 238/384 | 55 | 0.92 / 0.95 | 116/117 | 179ms |
| 20260312 | 85 (84) | 361/384 | 174 | 0.91 / 0.94 | 84/85 | 176ms |
| 20260317 | 88 (88) | 270/384 | 23 | 0.95 / 0.95 | 87/88 | 171ms |
| 20260318 | 85 (85) | 247/384 | 31 | 0.93 / 0.97 | 84/85 | 165ms |
| 20260320 | 100 (99) | 304/384 | 118 | 0.91 / 0.93 | 98/100 | 95ms |
| 20260323 | 95 (94) | 377/384 | 134 | 0.90 / 0.91 | 94/95 | 117ms |
| 20260324 | 99 (98) | 354/384 | 131 | 0.92 / 0.92 | 99/99 | 68ms |
| 20260325 | 99 (97) | 339/384 | 92 | 0.86 / 0.86 | 98/99 | 68ms |
| 20260326 | 101 (100) | 156/384 | 8 | 0.88 / 0.89 | 101/101 | 78ms |
| *chamber 2* | *1213 (1202)* | *3182/4608 (69%)* | *832* | | *1202/1213 (99%)* | |
| 20260819 | 102 (100) | 341/384 | 168 | 0.72 / 0.60 | 100/102 | 237ms |
| 20260821 | 102 (101) | 330/384 | 77 | 0.90 / 0.92 | 100/102 | 163ms |
| 20260825 | 115 (113) | 316/384 | 110 | 0.92 / 0.95 | 114/115 | 149ms |
| 20260903 | 129 (126) | 300/384 | 79 | 0.96 / 0.95 | 125/129 | 149ms |
| 20260904 | 128 (126) | 158/384 | 38 | 0.90 / 0.92 | 123/128 | 133ms |
| 20260909 | 110 (108) | 308/384 | 58 | 0.94 / 0.94 | 107/110 | 132ms |
| 20260910 | 107 (106) | 52/384 | 7 | 0.92 / 0.90 | 104/107 | 95ms |
| 20260915 | 106 (104) | 250/384 | 47 | 0.91 / 0.89 | 104/106 | 89ms |
| 20260916 | 117 (117) | 37/384 | 2 | 0.91 / 0.85 | 117/117 | 92ms |
| 20260918 | 101 (100) | 193/384 | 68 | 0.89 / 0.93 | 100/101 | 150ms |
| *chamber 3* | *1117 (1101)* | *2285/3840 (60%)* | *654* | | *1094/1117 (98%)* | |
| **Total** | **2330 (2303)** | **5467/8448 (65%)** | **1486** | | **2296/2330 (99%)** | |

Tuned counts can differ by a few from run to run because `computeDirectionTuning.m`'s permutation test (1000 random shuffles, no fixed seed) gives slightly different results for channels sitting right at the significance boundary — expected noise, not a bug.

**Of the 1486 "Tuned" channels, only 483 (33%) are cleanly single-peaked (`directional`); the remaining 1003 (67%) have a real dominant peak but ALSO real secondary/non-unimodal structure (`directional_complex`)** — see §4.7 for the 4-category classification this is based on. The split is similar in the two chambers (chamber 2: 256 clean / 576 complex-flagged; chamber 3: 227 / 427). It is not evenly distributed across areas; see §6.7 (some of the areas with the *highest* tuning fraction turn out to be driven mostly by `directional_complex`, not clean `directional`, channels).

**Latency.** In chamber 2, median saccade latency drops from ~170–200ms in the earlier sessions to ~70–120ms from 20260320 onward — plausibly a training/practice effect. Chamber 3 starts slow again (237ms on 20260819, the first chamber-3 session) and falls over the following weeks (163ms 20260821, 132–149ms through 20260909, 89–95ms in 20260910–20260916), with 150ms on 20260918; it does not return to the 68–78ms of late chamber 2. Whether the re-slowing reflects the new chamber/setup or the 5-month gap cannot be separated here.

**Eye calibration.** R² is 0.84–0.97 in 21 of 22 sessions. 20260819 is clearly weaker (0.72 x / 0.60 y; a fixation offset and overshoot at 180°/270° were seen in the trajectory plots) — saccade directions and onsets are fine there, but its absolute amplitudes in dva should not be trusted.

Saccade detection (velocity-threshold on a low-pass-filtered, calibrated eye trace, blink periods excluded) succeeded on 99% of trials (2296/2330); endpoint accuracy and main-sequence (peak velocity vs. amplitude) relationships were physiological in every session inspected (see `*_eye_summary.png`).

### 5.2 QC figures (per session, in `data/figures/`; all 22 sessions)

- **`*_overview.png`** — four panels, all sharing one channel-index y-axis with **channel 1 (deepest, nearest probe tip) at the bottom and the shallowest channel at the top**: (1) a real-depth ruler (electrode position in mm from the probe tip, so gaps/clusters from the day's online channel-selection scheme are visible directly), (2) the z-scored MUA heatmap with all 5 trial-epoch boundaries (fixation/stimulus/delay/go-cue/target-end) marked, computed per-session from that session's own `Stm` timing, (3) a per-channel, per-epoch significance map (visual/delay/peri/post, ms windows labeled on the axis) so you can see *which* window a channel's "responsive" flag refers to, (4) a text panel with that day's burr hole, planned target structures, and online channel-selection scheme (§2.2).
- **`*_eye_summary.png`** — 7 panels: trajectory spoke plots shown twice, once for all trials (error trials dimmed) and once for correct trials only with % correct in the title; both draw a light-grey dotted circle (radius = `Stm.TargWinSz`, the task's real hit-detection window) around each target, behind the traces; fixation scatter; main sequence; latency histogram (ms); pupil trace (ms); accuracy by direction.
- **`*_tuning_examples.png`** — all responsive channels, sorted by depth, one small tuning-curve tile per channel (min-max normalized per channel for display, with SEM error bars scaled by the same normalization factor — see §4.7), border and line color reflecting the 4-way category from §4.7: red=`directional`, purple=`directional_complex`, orange=`complex`, gray=`untuned`.
- **`*_tuning_matrix.png`** — all responsive channels, sorted by preferred direction, per-channel min-max normalized, plus a tuning-strength color strip and a 4-category color strip alongside.
- **`*_population_summary.png`** — generated per session (not pooled — each session is an independent penetration, often a different burr hole/target structure per §2.2, and the probe is not seated to the same depth across sessions). Preferred-direction rose plot and tuning-strength histogram, plus preferred direction vs. real electrode depth (mm from tip, this session only) — see §6.4.
- **`*_artifact_check.png`** — pairwise trial-residual correlation matrix across all channels (direction-tuning contribution removed), shown raw and sorted by real depth, plus mean correlation vs. inter-channel distance — see §6.4.
- **`*_shared_kinematics.png`** — per-channel PC1 loadings (the dominant shared trial-residual component) and its per-trial score plotted against saccade peak velocity and amplitude, plus (for the 4 sessions where this component is confirmed kinematics-linked) trial-residual correlation vs. distance before/after removing it — see §6.5.
- **`*_planexec_gocue.png` / `*_planexec_saclocked.png`** — per-channel effect map (visual/planning/execution, ms windows labeled), tuning-strength scatter (visual-vs-execution and planning-vs-execution), and preferred-direction scatter, for the go-cue-fixed and saccade-onset-locked window splits respectively (§4.8).
- **`*_planexec_gocue_curves.png` / `*_planexec_saclocked_curves.png`** — every responsive channel's actual tuning curve (with SEM error bars — see §4.7), stacked visual (top) / planning (middle) / execution (bottom), category-colored.
- **`*_decode_channels.png`** — individual-channel 8-way direction decoding accuracy by depth (visual/planning/execution), plus pairwise scatters (planning-vs-execution, visual-vs-execution) — see §4.9/§5.3.
- **`*_decode_population.png`** — population decoding accuracy vs. number of channels pooled, one curve per epoch, chance = 12.5% — see §4.9/§5.3.
- **`*_decode_genmatrix.png`** — 3×3 temporal generalization matrix (train epoch × test epoch), full responsive population, chance = 12.5% — see §4.10/§5.3.

### 5.3 Decoding results

Full-responsive-population decode accuracy per session (§4.9 method; chance = 12.5%):

| Day | Responsive channels | Trials | Visual | Planning | Execution |
|---|---|---|---|---|---|
| 20260224 | 64 | 85 | 15.3% | 14.1% | 17.6% |
| 20260226 | 271 | 144 | 29.2% | 28.5% | 43.8% |
| 20260303 | 201 | 112 | 14.3% | 17.9% | 18.8% |
| 20260310 | 238 | 116 | 30.2% | 25.9% | 36.2% |
| 20260312 | 361 | 84 | 29.8% | 26.2% | 33.3% |
| 20260317 | 270 | 88 | 17.0% | 8.0% | 15.9% |
| 20260318 | 247 | 85 | 27.1% | 18.8% | 15.3% |
| 20260320 | 304 | 99 | 38.4% | 32.3% | 45.5% |
| 20260323 | 377 | 94 | 31.9% | 33.0% | 44.7% |
| 20260324 | 354 | 98 | 35.7% | 35.7% | 57.1% |
| 20260325 | 339 | 97 | 53.6% | 38.1% | 44.3% |
| 20260326 | 156 | 100 | 9.0% | 14.0% | 9.0% |
| 20260819 | 341 | 100 | 56.0% | 31.0% | 56.0% |
| 20260821 | 330 | 101 | 50.5% | 22.8% | 39.6% |
| 20260825 | 316 | 113 | 45.1% | 33.6% | 42.5% |
| 20260903 | 300 | 126 | 46.0% | 27.0% | 53.2% |
| 20260904 | 158 | 126 | 38.9% | 23.0% | 36.5% |
| 20260909 | 308 | 108 | 33.3% | 21.3% | 33.3% |
| 20260910 | 52 | 106 | 19.8% | 12.3% | 7.5% |
| 20260915 | 250 | 104 | 34.6% | 29.8% | 36.5% |
| 20260916 | 37 | 117 | 12.8% | 13.7% | 18.8% |
| 20260918 | 193 | 100 | 26.0% | 24.0% | 21.0% |
| **Mean ± SEM, all 22** | | | **31.6 ± 2.9%** | **24.1 ± 1.8%** | **33.0 ± 3.2%** |
| Mean ± SEM, 18 sessions (excluding the 4 weak ones, §6.6) | | | 35.4 ± 2.7% | 26.5 ± 1.7% | 37.4 ± 3.0% |
| Mean ± SEM, chamber 2 (12) | | | 27.6 ± 3.6% | 24.4 ± 2.8% | 31.8 ± 4.5% |
| Mean ± SEM, chamber 3 (10) | | | 36.3 ± 4.4% | 23.8 ± 2.2% | 34.5 ± 4.8% |

**Which sessions decode above chance.** With 85–145 trials, an accuracy above ~20–22% is above the 99th percentile of a binomial with p=0.125, per session. By that criterion 16/22 sessions are above chance in the visual epoch and 15/22 in each of planning and execution. Six sessions are not above it in any epoch — 20260224, 20260303, 20260317, 20260326, 20260910 and 20260916 (the first four are chamber 2, the last two chamber 3; all but 20260303 and 20260317 are the weak sessions flagged in §5.1, §6.4–6.6) — and 20260318 is above it only in the visual epoch. Execution gives the highest mean accuracy (33.0% over all 22, 37.4% over the 18 non-weak sessions), consistent with a genuine, growing perimovement direction signal (matching the tuning-recruitment finding in §6.6), though visual and planning are also clearly above chance in most sessions — stimulus position and premotor direction are both decodable well before the saccade itself. Chamber 3's visual-epoch accuracy is higher than chamber 2's (36.3% vs 27.6%) while execution is similar (34.5% vs 31.8%); this comparison is across different brain regions and sample sizes and should not be over-read.

Pooling the population-size sweep across sessions (`pooled_decode_population.png`, restricted to population sizes shared by multiple sessions — see §4.9) gives a monotonic accuracy-vs-population-size curve for all three epochs: a single channel (n=1 random subsample) averages 13.6–14.1%, barely above the 12.5% chance level; the pooled mean reaches ~19% (visual/execution) and ~17% (planning) at n=8, 28%/22%/29% (visual/planning/execution) at n=64 (20 sessions), and **38.4% / 28.2% / 41.7%** at n=256 (12 sessions). Visual and execution track each other to within ~1 point up to n=128, with execution pulling ahead only at n=256; planning is the lowest curve at every population size. Individual channels are therefore weakly informative on their own and the direction information is distributed across many channels rather than carried by a few highly informative ones.

<!-- FIGURE: pooled/pooled_decode_population.png -->

**Temporal generalization matrix** (§4.10, `*_decode_genmatrix.png`), pooled mean across the 18 non-weak sessions (excluding 20260224/20260326/20260910/20260916, chance = 12.5%):

| Train \ Test | visual | planning | execution |
|---|---|---|---|
| **visual** | 35.4% | 19.2% | 21.9% |
| **planning** | 21.0% | 26.5% | 21.8% |
| **execution** | 21.2% | 20.6% | 37.4% |

Every off-diagonal cell is clearly above the 12.5% chance level — there is real shared, cross-epoch structure, a decoder trained in one epoch is not useless in another — but every off-diagonal cell is also well below its row's own diagonal (e.g. training on execution and testing on execution gets 37.4%, but training on execution and testing on visual only 21.2%). So the population code is **neither fully shared nor fully epoch-specific**: part of the direction-discriminating structure persists across visual/planning/execution, and part is reconfigured at each epoch.

**Sliding (continuous) temporal generalization matrix** (§4.10, `*_decode_slidegen.png`, 50ms bins across [-0.5, +0.8]s relative to go-cue; the diagonal of each matrix is the accuracy of decoding in that bin). In the 18-session mean, the diagonal shows **two bands**: an early one in the visual epoch (~31% at the first bins, −475 to −425ms, falling to ~25% by −325ms; the target is flashed around −500 to −400ms), a **lower plateau of ~20% through the delay period up to the go-cue** (still above chance — the target position remains decodable during the memory delay, at lower accuracy), and a **peri-saccadic band** starting right after the go-cue (28–31% from +75 to +325ms, peaking ~+125 to +225ms and tracking each session's own saccade latency, §5.1), declining to ~25% at +425ms, ~20% by +575ms and ~17% by the last bin (+775ms). Off-diagonal generalization is narrow (a decoder trained at one moment transfers well only to nearby moments), so the "reconfiguration" implied by the coarse matrix's below-diagonal off-diagonal cells is not a single sharp switch between 3 discrete codes but a population code that drifts/refines continuously through the trial, with the most stable structure concentrated around the visual response and the saccade itself. The weak sessions (20260224, 20260326, 20260910, 20260916) show diagonal values of at most ~20–25% anywhere.

**Decoding after common-mode removal** (§6.5's open item, closed): for the 4 kinematic-linked sessions (none of the 10 chamber-3 sessions qualifies, §6.5), decoding accuracy at the full responsive population, before vs. after removing that session-and-epoch's own dominant shared component (`decodeAfterCommonModeRemoval.m`, same removal method as §6.5):

| Day | Epoch | Original | Cleaned | PC1 var. explained |
|---|---|---|---|---|
| 20260310 | visual | 30.2% | 44.0% | 8.3% |
| 20260310 | planning | 25.9% | 31.0% | 22.2% |
| 20260310 | execution | 36.2% | 48.3% | 12.2% |
| 20260320 | visual | 38.4% | 59.6% | 18.6% |
| 20260320 | planning | 32.3% | 48.5% | 17.7% |
| 20260320 | execution | 45.5% | 62.6% | 14.9% |
| 20260323 | visual | 31.9% | 50.0% | 32.7% |
| 20260323 | planning | 33.0% | 44.7% | 32.1% |
| 20260323 | execution | 44.7% | 72.3% | 61.0% |
| 20260324 | visual | 35.7% | 51.0% | 11.3% |
| 20260324 | planning | 35.7% | 48.0% | 25.4% |
| 20260324 | execution | 57.1% | 65.3% | 19.0% |

**Decoding accuracy increases after removing the shared component in every epoch of every one of these 4 sessions** — often substantially (e.g. 20260323 execution: 44.7% → 72.3%). This answers the open question cleanly: these sessions' higher decoding accuracy is **not** because the shared/kinematic component was itself carrying decodable direction information (if it were, removing it would have *lowered* accuracy) — it was acting as noise that partially masked real direction-related structure. These 4 sessions are not artificially inflated; if anything they have *more* real signal than their original (pre-removal) numbers suggested. This also means common-mode removal is a genuinely useful preprocessing step, at least for sessions with a confirmed kinematic-linked component, not just a diagnostic check.

## 6. Important caveats — read before using the tuning results

Chamber-2 trajectories are mostly subcortical (§2.2) and chamber-3 trajectories are prefrontal cortex (§2.2). A per-channel trial-residual correlation check plus per-session population summaries provide evidence that at least some of the tuning is genuine, spatially-organized neural signal — but the picture is mixed, not a clean "all real" or "all artifact" result. Details below.

### 6.1 The high responsiveness fraction

65% of all channels across all 22 sessions (5467/8448) were flagged "responsive" — 69% in chamber 2 (3182/4608) and 60% in chamber 3 (2285/3840) — reaching 98% in one session (20260323, 377/384). Two things are relevant to interpreting this:

- **Online channel selection is not per-electrode, it's per bank/section**: `Kwibus_elab.pdf` logs an explicit operator-chosen scheme for every session — "quarter density," "4bankeven," "single column," "sections2/3," "4 bank," "Bank A" (§2.2 tables). The operator was already restricting acquisition toward electrode banks that looked responsive online, so a high responsive fraction *within the selected set* is an expected consequence of that selection.
- **Chamber-2 trajectories are mostly subcortical at depth** (§2.2): thalamus (MD/CL/CM), pulvinar, SC, habenula, ZI — not cortex. A depth-dependent *laminar* latency gradient is a cortical-columnar signature, not something to expect from a probe that's mostly in thalamus/SC by design. The chamber-3 trajectories are cortical, but the in-tissue span of those penetrations is short (4–16 mm estimated, §4.13) and the area labels are cluster-level (§4.12), so no laminar analysis is attempted here.

See §6.4 for a direct test of the shared-artifact possibility.

### 6.2 Preferred-direction clustering across sessions

The better-tuned areas in both chambers have dominant preferred directions in the same ~115–160° range (§6.7: SC 115°, MRT 131°, ventral area 46 129°, 8Ad 146°, area 12m 153°, MD 157°) despite being very different structures in different chambers. That is what a consistent contralateral-field bias would produce if both chambers are in the same hemisphere (the logs read for this report do not state the hemisphere of chamber 3 — to be confirmed), but a shared direction-dependent contaminant such as movement-related signal (§6.4–6.5) would produce a similar pattern, so this does not by itself establish that these are retinotopic/motor maps.

### 6.3 Cross-session depth comparison

Sessions are independent penetrations, the probe is not seated to the same depth across sessions, and (§2.2) different sessions often target different burr holes/structures. The population summary is generated **per session** (`*_population_summary.png`), using real electrode depth in mm from the probe tip rather than an abstract cross-session channel index. All screening/tuning computation is strictly within-session (`tuning_results.mat` keyed by `Day`/`RunN`/`channel`). Cross-session comparison is done only through the area labels (§6.6–6.8).

### 6.4 Per-channel trial-residual correlation check (the direct artifact test)

`checkArtifactCorrelation.m` tests the shared-artifact hypothesis directly: for each channel pair, correlate the **trial-to-trial residual** peri-saccadic response (each direction-condition's own mean subtracted out first, so genuine shared direction tuning can't inflate this) against physical inter-channel distance. A uniform artifact riding on every channel predicts high correlation that does **not** decay with distance; local neural activity predicts correlation that decays with distance.

Mean pairwise correlation per session (`*_artifact_check.png`; sessions in date order, chamber 2 then chamber 3):

| Day | Mean corr | Day | Mean corr | Day | Mean corr |
|---|---|---|---|---|---|
| 20260224 | 0.078 | 20260323 | 0.235 | 20260904 | 0.060 |
| 20260226 | 0.112 | 20260324 | 0.100 | 20260909 | 0.034 |
| 20260303 | 0.203 | 20260325 | 0.057 | 20260910 | 0.038 |
| 20260310 | 0.099 | 20260326 | 0.021 | 20260915 | 0.061 |
| 20260312 | 0.213 | 20260819 | 0.089 | 20260916 | 0.042 |
| 20260317 | 0.078 | 20260821 | 0.046 | 20260918 | 0.058 |
| 20260318 | 0.179 | 20260825 | 0.081 |  |  |
| 20260320 | 0.077 | 20260903 | 0.053 |  |  |

Two findings, looking at `*_artifact_check.png` and (for chamber 2) `*_population_summary.png`:

1. **Correlation decays with distance in all 22 sessions** (first distance bin higher than the last bin everywhere; e.g. 20260323: ~0.55 at <200µm down to ~0.18 by 1mm+). That is the signature of local neural activity rather than a probe-wide shared signal. The decay is weak in the weakest sessions (20260916: 0.06 → 0.05, essentially flat; 20260326: 0.08 → 0.02; 20260910: 0.08 → 0.03). In the better sessions it does not decay all the way to zero — it plateaus at a small positive value (~0.1–0.2 in the higher-correlation sessions) even between far-apart channels, consistent with a smaller residual shared component (common reference, EMG, or motion) riding underneath genuine local structure. Chamber 3 has uniformly low mean correlations (0.034–0.089, versus 0.021–0.235 in chamber 2); the largest chamber-3 value is 20260819 (0.089).
2. **The chamber-2 per-session population summaries show genuine, session-specific spatial organization that a shared artifact would not produce**: 20260312 (burr hole 1, SC) shows a repeating sawtooth relationship between real depth and preferred direction across ~4 distinct depth clusters; 20260323 shows a smooth, continuous, near-monotonic depth-direction relationship; 20260320 (burr hole 4, thalamus) shows depth-banded clusters of consistent preferred direction. By contrast, **20260326** — the session with the lowest mean pairwise correlation (0.021) and lowest tuned fraction (8/156) — shows no depth structure at all, just scattered noise. The chamber-3 population-summary figures were generated for all ten sessions but were not individually reviewed for this report.

**Net assessment**: several chamber-2 sessions (20260312, 20260320, 20260323 especially) show good evidence of genuine, spatially-organized tuning. A smaller shared/common-mode component likely coexists across all sessions. 20260326 and, in chamber 3, 20260910 and 20260916 (52 and 37 responsive channels, 7 and 2 tuned; 20260916's log notes no obvious signal at any depth) look like they have little usable signal. Recommend treating tuning results session-by-session rather than as one uniform "1486 tuned channels" figure.

### 6.5 Does the shared component track saccade kinematics?

Direct follow-up to §6.4. `checkSharedComponentKinematics.m` extracts the dominant shared fluctuation — the first principal component (SVD) of the same trial-residual matrix used in §6.4 — and correlates its per-trial score against saccade peak velocity and amplitude (Spearman). A real, consistent correlation with both kinematic measures points to EMG/movement contamination; no correlation argues against that specific explanation (without ruling out a different shared source, e.g. a reference issue).

| Day | PC1 var. explained | Loading sign agreement | r(peakVel), p | r(amplitude), p |
|---|---|---|---|---|
| 20260224 | 24% | 88% | 0.11, 0.33 | 0.05, 0.68 |
| 20260226 | 16% | 94% | 0.03, 0.74 | -0.02, 0.80 |
| 20260303 | 32% | 99% | -0.08, 0.39 | -0.04, 0.71 |
| 20260310 | 12% | 92% | **0.25, 0.007** | **0.25, 0.007** |
| 20260312 | 74% | 90% | -0.26, 0.017 | -0.14, 0.20 |
| 20260317 | 10% | 88% | -0.17, 0.11 | -0.13, 0.25 |
| 20260318 | 26% | 98% | 0.10, 0.35 | 0.04, 0.74 |
| 20260320 | 13% | 82% | **0.27, 0.007** | **0.31, 0.002** |
| 20260323 | 52% | 91% | **0.21, 0.044** | **0.28, 0.006** |
| 20260324 | 19% | 79% | **0.45, <0.001** | **0.30, 0.003** |
| 20260325 | 25% | 81% | -0.07, 0.49 | -0.00, 1.00 |
| 20260326 | 14% | **55%** | 0.02, 0.85 | 0.03, 0.75 |
| 20260819 | 26% | 86% | -0.07, 0.47 | 0.05, 0.61 |
| 20260821 | 15% | 81% | 0.07, 0.51 | 0.03, 0.78 |
| 20260825 | 22% | 88% | 0.08, 0.38 | 0.07, 0.44 |
| 20260903 | 7% | 93% | -0.11, 0.21 | -0.10, 0.29 |
| 20260904 | 11% | 88% | -0.19, 0.039 | -0.14, 0.13 |
| 20260909 | 9% | 83% | -0.13, 0.18 | -0.15, 0.12 |
| 20260910 | 13% | **62%** | -0.10, 0.29 | -0.15, 0.12 |
| 20260915 | 14% | 90% | -0.18, 0.073 | -0.17, 0.096 |
| 20260916 | 6% | 74% | -0.09, 0.34 | -0.10, 0.28 |
| 20260918 | 18% | 77% | -0.03, 0.76 | -0.04, 0.70 |

1. **PC1 (the dominant shared fluctuation) explains a lot of the trial-to-trial residual variance in most sessions** — from 6% (20260916) up to 74% (20260312) — and loading sign agreement across channels is high (74–99%) in nearly every session. The exceptions are **20260326 (55%, essentially chance)** and **20260910 (62%)**, meaning those sessions' residual variance isn't dominated by one common-mode axis at all, consistent with §6.4's finding that they have no organized structure.
2. **The kinematic correlation is real but not universal: exactly 4/22 sessions** (20260310, 20260320, 20260323, 20260324 — all chamber 2; none of the 10 chamber-3 sessions) show a significant, same-sign correlation with **both** peak velocity and amplitude — the signature of EMG/movement contamination. The other 18 show no reliable relationship. Two sessions show a significant correlation with peak velocity only and with a *negative* sign (20260312: r = −0.26, p = 0.017; 20260904: r = −0.19, p = 0.039; neither reaches significance for amplitude); with 22 sessions × 2 measures tested without correction, an occasional nominal p<0.05 is expected by chance.

**Interpretation**: the shared component's origin is not uniform across the dataset. In 4 chamber-2 sessions there's direct, quantitative evidence it's kinematics-linked (EMG/motion-like); in the rest (including all of chamber 3), either it has a different origin, or the correlation approach lacks power in sessions with less kinematic variance/trial count.

**Common-mode removal, for the 4 kinematic-linked sessions.** For 20260310, 20260320, 20260323, and 20260324, `checkSharedComponentKinematics.m` also subtracts the rank-1 PC1 reconstruction from the trial residual and recomputes §6.4's correlation-vs-distance relationship on the cleaned data (`*_shared_kinematics.png`, 4th panel, before/after overlay). In all 4 sessions, the long-distance correlation **plateau** (§6.4's evidence for a residual shared component) drops sharply after removal — e.g. 20260323 goes from a ~0.18–0.21 plateau at >1mm down to ~0.05–0.07 — while the short-distance correlation, though also reduced, remains clearly above zero. This confirms the plateau specifically was this kinematic component, and that at least some of the short-distance, depth-organized correlation in §6.4 is not attributable to it.

### 6.6 Planning-vs-execution tuning, across all 22 sessions

`summarizePlanningExecution.m` quantifies the saccade-onset-locked planning/execution split (§4.8) across all 22 sessions. Four sessions have at most 3 directional channels in either window and are excluded from the pooled stats below: 20260224 (0 planning / 1 execution), 20260326 (2 / 1), 20260910 (0 / 1) and 20260916 (3 / 0); the next-weakest session kept is 20260303 (3 / 11). The same four are excluded from the by-area version below (`poolPlanExecByArea.m`).

- **Pooled across the other 18 sessions: 1324 channel-instances were `directional` in execution; only 47% (619) of those were also `directional` during planning — 53% (705) only became directional during execution**, i.e. genuine recruitment of additional direction-selective activity around the movement, not just an amplitude increase on an already-tuned population.
- Conversely, of 1076 `directional`-in-planning channel-instances, 58% (619) remained `directional` into execution (42%, 457, were planning-only).
- **For channels tuned in both epochs, preferred direction is meaningfully consistent, not independent**: mean absolute angular difference (planning vs execution), pooled = **37.6°**, well below the ~90° expected for unrelated/random pairings, though far from 0°. Per-session values range from 13.0° (20260825) and 14.2° (20260325), very consistent, to 63–80° (20260903: 62.9°, 20260915: 73.1°, 20260303: 79.5°), closer to random. Note (§4.7): `directional` here pools both the `directional` and `directional_complex` sub-categories, so this average mixes clean single-peak `prefDir` estimates with dominant-peak-plus-secondary-structure ones for which `prefDir` is a less clean summary — the true consistency of the *unambiguous* single-peak population specifically has not been separately checked.

**The pooled 53% execution-only figure hides substantial area-level heterogeneity.** `poolPlanExecByArea.m` redoes this split BY AREA (same liberal, single-level pooling as §6.7, run at both `L6` and `L3`) instead of collapsing across all channels. At `L6`, %execution-only ranges from **83% (AREA_46F, n=29 execution-directional; 82% FUNDUS_IPS, 80% AREAS_1-2, 79% AREA_23)** down to **8% (CMN-PF, n=12)**; at `L3`, from **82% (VM_IPS)** and 81% (SI) down to **13% (LMID)**. Points worth calling out, consistent across both grains where the grains overlap:
- **MRT (`L6`) / LMID (`L3`) — the same adjacent-midbrain channels, and the highest tuning-fraction area in §6.7 (49.3% directional) — shows the lowest execution-only fraction of any area with at least 20 execution-directional channel-instances at both grains (`L6` MRT: 13.3%, n=60; `L3` LMID: 13.3%, n=60, identical as expected)** (CMN-PF is lower, 8.3%, but with only 12). This population's direction selectivity is mostly *already present during planning*, not newly recruited at movement time — a real, distinct finding from §6.7's tuning-fraction result: high tuning fraction here reflects a stable, persistent direction code, not a code that emerges around the saccade.
- **Posterior cortex sits at the opposite end** (FUNDUS_IPS/VM_IPS, AREAS_1-2/AREA_23 at `L6`, SI/PCGG at `L3`: 77–82% execution-only) — channels here that show direction tuning at all are almost entirely recruited during the movement itself.
- **Prefrontal cortex is intermediate and heterogeneous** (new in chamber 3): dorsal frontal areas are mostly recruited at execution (dorsal periarcuate 8A, 8Ad/AREA_8A: 68.0%, n=125; dorsal area 46, 46d/DLPFC: 60.6%, n=99), whereas ventral area 46 is closer to the thalamic/midbrain pattern (AREA_46V: 36.5%, n=203 — it has *more* planning-directional (230) than execution-directional channel-instances; VLPFC at `L3`: 42.2%, n=232). Area 46f (82.8%) and medial area 12 (40.0%) rest on only 29 and 20 execution-directional instances.
- SC (`L6`, 53.5%) / CO (`L3`, 52.5%) sits roughly in the middle, and again consistent across grains — apparently a genuine mix of persistent and newly-recruited direction-selective channels, not clearly at either extreme.
- AREA_29 (8 execution-directional instances at `L6` in the 12-session analysis) no longer meets the ≥8 threshold and is not reported.

Full per-area tables and figures saved to `data/figures/pooled/` and `data/planexec_pooled_by_area_L6.mat`/`_L3.mat` (`planexec_by_area_L6.png`, `planexec_by_area_L3.png`). Same caveats as §6.7/§6.8 apply (liberal within-level pooling, small per-area sample sizes for some areas — reported only where an area has ≥8 execution-directional channel-instances).

<!-- FIGURE: pooled/planexec_by_area_L6.png -->

<!-- FIGURE: pooled/planexec_by_area_L3.png -->

### 6.7 Pooling responsiveness/tuning by area, across all 22 sessions

Direct follow-up to §6.1/§6.3, now that channels carry an actual (tentative) area label (§4.12) rather than just an opaque depth index: `poolTuningByArea.m` pools every channel across ALL 22 sessions **by area** instead of by session — the exact thing §6.3 said couldn't be done responsibly otherwise ("sessions are independent penetrations... not seated to the same depth"). Method (liberal multi-label pooling, run separately per hierarchy level): §4.14. Results below are `L6` (finest) and `L3` (coarsest), run separately to bracket the range rather than blur it. Chamber-2 areas pool only chamber-2 channels and chamber-3 areas only chamber-3 channels (the two chambers share no area label), so the earlier chamber-2 numbers are unchanged (tuning and decoding by area reproduced exactly when re-run with the chamber-3 data present).

**Highest pooled tuning fractions are midbrain and thalamic, with ventral prefrontal cortex the highest cortical value.** At `L6`: MRT 49.3% directional (n=140, 99% responsive), **SC 40.8%** (n=620, 91.0% responsive), medial area 12 (AREA_12M) 39.3% (n=84, one session, 20260821), F (fornix label) 35.9%, CL-PC 35.0%, MD 32.3% (n=520), **ventral area 46 (AREA_46V) 30.6%** (n=941, 79.4% responsive; 288 directional channels, the most of any `L6` area; SC has 253), SG 25.7%. At `L3`: LMID (the same channels as MRT) 49.3%, CO (colliculi, same channels as SC) 42.6% (n=580), MG 42.5% (n=80), F 35.9%, MTHAL 31.2% (n=420), LAT_OFC 28.0% (n=168), **VLPFC 27.6%** (n=1165 = 46v + 46f). SC (`L6`) and CO (`L3`) agree because they are largely the same physical channels described at two granularities. Most other areas sit at 3–18% directional at either grain despite often being similarly or more responsive (e.g. `L3` PCGG: 43.6% responsive, only 5.6% directional) — high responsiveness does not, on its own, imply direction tuning.

**Within prefrontal cortex the tuning fraction is strongly area-dependent, and not what an a priori oculomotor ranking would predict.** Ventral area 46 (46v) is high (30.6%) while the dorsal areas, which include the FEF region, are low: dorsal periarcuate 8A (8Ad; n=1208, 57.4% responsive) is only 12.9% directional, dorsal area 46 (46d; n=1299, 45.3% responsive) 9.9%; area 46f (fundus of the principal sulcus) is 15.2% (n=224), lateral area 13 16.7% (n=84, one session, 20260819). Possible contributors that this analysis cannot separate: genuine area differences; which parts of each area these penetrations sampled; the shallow in-tissue spans and uncertain depths in chamber 3 (§4.13); and white-matter-contaminated channels pooled as cortical (§4.12) — 8 of the 36 `dorsal_periarcuate_area_8A` clusters carry the "might contain WM" note, and these channels are included in the 8Ad pool.

**The areas with the highest tuning fraction are driven mostly by `directional_complex`, not clean `directional`, channels — a different (per-channel, response-shape) axis from the population-level preferred-direction distribution below.** §4.7 introduced a 4th category (`directional_complex`: a real dominant peak that ALSO has real secondary/non-unimodal structure, as opposed to plain `directional`'s clean single peak), since about two thirds of all pPerm-significant channels turn out to have this property. That split is *not* evenly spread across areas: `SC`/`CO`'s 40.8%/42.6% "directional" figure breaks down to ~10% clean vs. ~31-32% `directional_complex` (i.e. **only about 1 in 4 of SC/CO's "directional" channels is cleanly single-peaked**). `MRT`/`LMID` (the highest tuning fraction overall, 49.3%) is even more skewed: ~6% clean vs. ~44% `directional_complex` (about 1 in 9 clean). Ventral area 46 resembles SC (7.0% clean vs. 23.6% complex-flagged, ~1 in 4 clean), as does medial area 12 (11.9% vs. 27.4%). The lower-fraction frontal areas are more balanced — 8Ad 6.3% clean vs. 6.6% complex-flagged, 46d 3.9% vs. 6.0%, 46f 7.6% vs. 7.6% — and smaller-n areas like `AREA_5D` (10.2% directional, but **100%** of that is clean) show the opposite pattern. So "stands out" in §6.7 should be read as "stands out for having SOME kind of real direction structure," not "stands out for a large population of cleanly single-peaked, retinotopy-like tuning" — the latter is a much smaller subset. Full breakdown (`%DirClean`, i.e. the clean sub-fraction of `%Dir`) is in `tuning_pooled_by_area_{L6,L3}.mat`, not separately plotted.

**Preferred-direction distribution, not a single "preferred direction per area."** There's no a priori reason a whole area should have one preferred direction: for a retinotopically organized structure like SC/CO, preferred direction should track WHERE in the map a given penetration happened to land, not converge to a single value — and any pooled distribution here is itself shaped by which part of an area each of the (few) penetrations sampled, not necessarily an intrinsic property of the area. `pooled_tuning_by_area_rose_L6.png`/`_L3.png` therefore show the actual **distribution** of individual channels' preferred directions per area (8 bins, matching the task's 8 directions) rather than collapsing to one resultant vector — the resultant-length-weighted circular mean is kept as a dotted reference line for continuity with the summary table, but it is not the headline result. Several areas show a clearly unimodal, skewed distribution: SC/CO strongly concentrated around ~115° (resultant concentration 0.84–0.86 over ~250 channel-instances at either grain), MRT/LMID 131° (0.91, n=69), and in prefrontal cortex ventral area 46 (46v; 129°, concentration 0.51, n=288) and 8Ad (146°, 0.66, n=156) are both concentrated in the ~120–180° range, with medial area 12 (153°, 0.88, n=33) similar. Others show a visibly **multi-lobed** shape that a single-vector summary would have misrepresented as "just noisy": dorsal area 46 (46d; concentration 0.15 over 129 channels, with lobes near 30–60°, ~150° and ~210°) and the posterior-cingulate channels (`L3` PCGG, 0.15; `L6` AREA_23, 0.23, n=31). Given only 22 sessions in total and each area typically sampled by 1–6 of them, none of these distributions should be read as a claim about the area's true population-level direction organization — only about what these particular penetrations happened to pick up.

Full per-area tables and figures saved to `data/figures/pooled/` and `data/tuning_pooled_by_area_L6.mat`/`_L3.mat`: `pooled_tuning_by_area_summary_{L6,L3}.png` (%responsive/directional/complex per area, with channel count `n=` in each label), `pooled_tuning_by_area_rose_{L6,L3}.png` (preferred-direction distribution per area, areas with ≥15 directional channels only), `pooled_tuning_by_area_legend_{L6,L3}.png` (abbreviation → full descriptive name → channel count). Given the liberal (multi-label, within-level) pooling, treat all of this as **evidence, not proof** — a channel's inclusion in a given area's pool reflects that area being one *plausible* candidate for that channel, not a confirmed identity.

<!-- FIGURE: pooled/pooled_tuning_by_area_summary_L6.png -->

<!-- FIGURE: pooled/pooled_tuning_by_area_rose_L6.png -->

<!-- FIGURE: pooled/pooled_tuning_by_area_legend_L6.png -->

<!-- FIGURE: pooled/pooled_tuning_by_area_summary_L3.png -->

<!-- FIGURE: pooled/pooled_tuning_by_area_rose_L3.png -->

<!-- FIGURE: pooled/pooled_tuning_by_area_legend_L3.png -->

### 6.8 Decoding by area: does the tuning-fraction result generalize?

Direct follow-up to §6.7: does an area's high directional-tuning fraction also show up as better population *decoding*, or is that a separate question? `decodeByArea.m` restricts each session's own decoding population (§4.9's nearest-centroid/LOO method) to just that area's responsive channels in that session, size-matched at 20 channels across sessions for a fair cross-area comparison (method: §4.14). Unlike §6.7, this can't pool raw channels across sessions (trial-level data from different sessions isn't poolable the way per-channel tuning summaries are) — each session contributes one accuracy number per area per epoch, averaged across qualifying sessions (2–6 sessions per area; chance = 12.5%). Run at both `L6` and `L3`, same reasoning as §6.7.

**The tuning-fraction result does not simply generalize, at either grain.** SC (`L6`: 40.8% directional in §6.7) decodes only moderately at a matched population size — visual 20.2%, planning 22.3%, execution 23.4% (4 sessions) — clearly above chance but below several areas; CO (the same channels at `L3`) decodes almost identically (19.9 / 22.5 / 23.7%), a good consistency check across grains. The best decoders at `L6` by execution accuracy are **MD (26.8 / 21.3 / 30.2%, 5 sessions), ventral area 46 (31.6 / 20.4 / 28.8%, 4 sessions), MRT (32.4 / 29.7 / 28.5%, 2 sessions), CL-PC (24.0 / 23.5 / 27.6%) and F (27.1 / 23.6 / 26.6%)**; MRT/LMID has the highest tuning fraction (49.3%) and the best visual- and planning-epoch decoding at both grains (`L3` LMID: 31.4 / 30.3 / 28.1%), consistent with its tuning-fraction lead, and `L3` MTHAL (26.9 / 23.2 / 31.2%) has the best execution-epoch accuracy at that grain. In prefrontal cortex, decoding follows the tuning fractions: **ventral area 46 decodes well** — and, notably, best in the visual epoch (visual 31.6% vs. planning 20.4%; `L3` VLPFC 30.7 / 20.7 / 27.9%), i.e. the stimulus position is decodable from this population more readily than the planned saccade direction — whereas the low-tuning-fraction dorsal frontal areas decode only modestly (46d 19.8 / 16.1 / 19.8%, 6 sessions; 8Ad 20.2 / 16.8 / 19.1%, 4 sessions; `L3` DLPFC 20.6 / 15.9 / 20.5%). Where tuning fraction and decoding disagree: CMN-PF decodes about as well as SC (21.3 / 22.2 / 24.6%) with roughly half SC's directional fraction (21.7%), and posterior cingulate (AREA_23, 4.8% directional) still decodes above chance (14.8 / 15.2 / 18.4%) — a channel population can apparently carry substantial pooled direction information without many of its individual channels reaching the strict single-peaked `directional` threshold, plausibly via the `complex`/non-unimodal responses `classifyTuning.m` was built to detect (§4.7) but which a simple %-directional count doesn't credit. So: %directional is a real but incomplete predictor of decodability, and this holds regardless of which hierarchy grain is used.

Full per-area tables and figures saved to `data/figures/pooled/` and `data/decode_by_area_L6.mat`/`_L3.mat` (`decode_by_area_L6.png`, `decode_by_area_L3.png`). Same liberal (within-level) multi-label pooling caveat as §6.7 applies, and each area rests on only 2–6 sessions.

<!-- FIGURE: pooled/decode_by_area_L6.png -->

<!-- FIGURE: pooled/decode_by_area_L3.png -->

## 7. Open items requiring your input

1. **Artifact investigation** (§6.1, §6.4, §6.5) — the shared-component/kinematics check gives a real but non-universal answer: 4/22 sessions (all chamber 2) show a clean EMG/movement signature, the rest don't. Not pursued further: the weakest sessions overall (chamber 2: 20260224, 20260326; chamber 3: 20260910, 20260916) are plausibly explained simply by penetration placement — these burr holes/trajectories were planned for a different, unrelated paradigm, not for saccade responsiveness specifically (§2.2), so a session landing in less saccade-sensitive territory is an expected outcome of the shared multi-paradigm penetration design, not necessarily a data-quality problem to chase. (20260916's log itself notes "no obvious signal at any depth".)
2. **Absolute per-channel depth-in-brain** (mm below cortical surface) is not directly measured in the logs (§2.2) — only relative depth along the probe. **Done, two ways, both estimates rather than verified measurements**: (a) per-channel tentative area labels, hand-assigned against an atlas (§4.12), available via `S.channelArea.AssignedArea`/`.L6`-`.L3`; (b) an estimated absolute `DepthBelowSurfaceMm` per channel (§4.13), anchored from an explicit final depth where logged (20260324/20260325), the Chamber2Holes planned depth for the other chamber-2 sessions, and the elab-logged in-tissue depth (final − signal onset) for chamber 3. **No depth is available for 20260904, 20260915 and 20260916** (signal onset blank/unknown in the elab log).
3. **Burr hole for 20260224 is genuinely ambiguous** (marked "5?" in the source elab log itself, §2.2) — confirmed not worth chasing further.
4. **Planning-vs-execution** (§4.8, §6.6) — if useful, the natural extension is checking whether the 53% execution-only recruitment fraction varies with burr hole/structure (§2.2) or with the §6.5 kinematic-link sessions.
5. **Responsiveness/tuning thresholds**: current thresholds (FDR 0.05, effect size 0.3 baseline SD for screening; permutation p<0.05 for tuning) are reasonable defaults; given §6.4's session-by-session variability, you may want to treat sessions individually rather than adjust one global threshold.
6. **Done** (§5.3): decoding after common-mode removal, for the 4 kinematic-linked sessions — accuracy *increases* in every epoch of every session, meaning the shared component was masking real signal, not itself carrying decodable direction information.
7. **Done** (§5.3): continuous sliding-window temporal generalization matrix — two bands on the diagonal (visual response, then a peri-saccadic band), with a lower, still above-chance plateau through the delay.
8. **Done** (§6.7, §6.8): pooled tuning and decoding by area across all 22 sessions, now including prefrontal cortex.
9. **Chamber-3 anatomy is hand-assigned and not yet atlas-verified.** The area labels (§4.12) come from the planning documents and by-eye atlas comparison; the trajectories are being checked against the atlas planning files separately, and any change to the hole/structure plan (`holeInfo3`/`T3` in `sessionAnatomyInfo.m`) or to cluster labels should be propagated by editing `anatomy/chamber3/` and re-running the area-pooling analyses. White-matter channels are pooled with the nearest cortical area (21 clusters flagged "might contain WM"); a sensitivity analysis excluding the flagged clusters has not been run.
10. **20260819's eye calibration is weak** (R² 0.72/0.60, §5.1): directions and onsets are fine, absolute amplitudes in dva are not trustworthy.
11. **Hemisphere/side of chamber 3** is not stated in the documents used here; §6.2's contralateral-bias reading depends on it.
12. **Chamber-3 population-summary figures** (`*_population_summary.png`, depth vs. preferred direction) were generated for all ten sessions but not individually reviewed in this report (§6.4).

## 8. Recommended next steps

1. **Given §6.4/§6.5's mixed result**, treat sessions individually rather than pooling a single "tuned channel" count — prioritize chamber-2 sessions 20260312/20260320/20260323/20260324 (clean depth-organized structure and/or confirmed kinematic-linked shared component) and the chamber-3 sessions with the most tuned channels (20260819, 20260825, 20260903, 20260821) for any further biological claims, and treat 20260224/20260326/20260910/20260916 (weakest on nearly every metric) with more skepticism.
2. **Done** (§6.5): common-mode removal, for the 4 kinematic-linked sessions — the long-distance correlation plateau drops sharply after removing the shared component, while short-distance correlation survives.
3. **Done** (§4.13): absolute per-channel depth-in-brain, estimated per session — see §7 item 2.
4. Once sessions are triaged by §6.4/§6.5, revisit direction tuning with a von Mises fit (not just vector sum) for a cleaner tuning-width estimate on the trustworthy sessions — the saccade-onset-locked alignment from §4.8 is already implemented, use it in place of go-cue alignment for this.
5. Check whether the 53% execution-only recruitment fraction (§6.6) varies systematically with burr hole/structure or with the §6.5 kinematic-link sessions; the new prefrontal result (dorsal frontal areas mostly recruited at execution, ventral area 46 closer to the midbrain/thalamic persistent-code pattern) is a natural starting point.
6. **Done** (§4.9, §5.3): direction decoding, individual-channel and population, all 3 epochs, all 22 sessions.
7. **Done** (§5.3): decoding after common-mode removal, and the continuous sliding-window generalization matrix.
8. After the atlas check of the chamber-3 trajectories (§7 item 9): update labels, re-run the area-pooling analyses (L6 and L3), and consider the white-matter sensitivity analysis.
