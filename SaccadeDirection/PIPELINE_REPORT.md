# Saccade Direction Pipeline — Report

Subject: m032 ("Kwibus"). Task: 8-target visually-guided saccade task, run under the shared curve-tracing ("TDCT"/"ct") behavioral codebase with `Stm.difficulty = 1` (single target, no curve, no distractor) and stimulus file `stim_ct_classic_8targ.m`. Data source: `/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_collection/m032`.

All code lives in `SaccadeDirection/code/`. All extracted data and figures live in `SaccadeDirection/data/` and `SaccadeDirection/data/figures/`.

---

## 1. Objective

Extract raw Neuropixels + eye-tracking data for every saccade-direction session, reduce the probe data to a multi-unit activity (MUA) envelope per channel per trial, screen every channel for task responsiveness, and compute saccade-direction tuning for responsive channels — with the eye-tracking data used both as a behavioral readout in its own right (fixation stability, saccade kinematics, pupil) and as the ground truth for calibrating and QC-ing the neural alignment.

## 2. Identifying the task and the sessions

Per your confirmation, a run is a saccade-direction run if and only if its folder contains `stim_ct_classic_8targ.m`. Scanning all of `m032`'s day folders for this marker found **12 runs**, one per day, spanning 2026-02-24 to 2026-03-26:

| Day | Run | Day | Run | Day | Run |
|---|---|---|---|---|---|
| 20260224 | run-001 | 20260312 | run-003 | 20260323 | run-003 |
| 20260226 | run-003 | 20260317 | run-003 | 20260324 | run-003 |
| 20260303 | run-003 | 20260318 | run-003 | 20260325 | run-003 |
| 20260310 | run-003 | 20260320 | run-003 | 20260326 | run-003 |

The stimulus file confirms this is a genuine 8-direction task: `Stm.Targ(1..8).xy_dva` places targets at 0°, 45°, 90°, ..., 315° at 12 dva eccentricity, with `Stm.difficulty = 1` disabling the curve/distractor machinery inherited from the curve-tracing paradigm, leaving a plain visually-guided saccade: fixate → brief target/curve flash (100ms) → 400ms delay → go-cue (target reappears / fixation point change) → saccade.

Older 2024/2025 sessions in the same folder tree use the same task family (`task-ct`) but only the curve-tracing variants (`fgdots`/`classic` without `_8targ`), so they were correctly excluded.

### 2.1 Run → OpenEphys recording → behavioral log mapping

Each behavioral run folder is matched to its OpenEphys `experimentN/recordingM` folder by `matchRunRecording.m` (closest absolute start time, §4.2) and to its own behavioral log `.mat` by `findSaccadeRuns.m` (the `sub-*.mat` file inside the run folder itself). The verified pairing for all **12** extracted sessions:

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

Each `_extracted.mat` now stores this verified path as `S.LogFile`. **Do not use `S.Par.LogFn`/`S.Par.LogPath` for this purpose** — those fields come from inside the log file's own saved `Par` struct and were found to be stale in every session checked (they point at an unrelated earlier `task-rfmap`/`stim-sparsenoise` run in the same stim-program session, not the 8-target run the data was actually extracted from). This was corrected on 2026-07-08: `extractSaccadeSession.m` now saves `S.LogFile` explicitly, and the (then-11) existing extracted files were patched in-place with the correct value (`patchExistingExtractions.m`).

**Correction (2026-07-08): `20260224/run-001` was wrongly excluded — it does have a matching recording.** The original `matchRunRecording.m` only globbed `dayPath/experiment*/recording*`, i.e. direct children of the day folder. On 20260224 the *actual* recordings live in a separate, differently-named OpenEphys output folder sitting *alongside* `experiment1` rather than inside the searched tree: `20260224/m032-2026-02-24_12-59-28/Record Node 110/experiment1/recording4`, whose `sync_messages.txt` start time is **13:45:51 — only 10.5s after** the run-001 folder timestamp (13:45:41), well inside the 90s matching tolerance. The non-recursive search never looked there, so the "closest recording" it found was 72 minutes away (from the unrelated top-level `experiment1`), which correctly (given only what it searched) reported no match. Fixed by making `matchRunRecording.m` search recursively (`dir(fullfile(dayPath,'**','recording*'))`) under the day folder; verified this doesn't change the match for any of the other 11 sessions. Re-extracted 20260224/run-001 (85 trials, 85 correct) and re-ran screening/tuning/eye-analysis across all 12 sessions — results below now include it. Only 20260224 has this alternate folder layout; the other 11 days were unaffected.

### 2.2 Burr hole, channel selection, and planned target structure per session

`Chamber2Holes.docx` (MRI-based trajectory planning) documents 5 burr holes (0, 1, 2, 4, 5), each with a planned insertion depth and a list of structures its trajectory passes through, superficial to deep. It does **not** say which day used which hole. `Kwibus_elab.pdf` (the day-by-day lab notebook, 35 pages) does: every recording day logs a burr hole and an online channel-selection scheme (e.g. "quarter density," "4bankeven," "single column," "sections2/3" — this is the operator manually restricting acquisition to part of the probe based on what looked good online, confirming your point in review: channel selection is not per-electrode, it's per bank/section). Cross-referencing the two, all 12 sessions:

| Day | Burr hole | Channel selection (online, from elab log) | Planned structures (superficial → deep, from Chamber2Holes.docx) |
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

This is directly relevant to §6: **at planned full depth, every one of these trajectories is mostly or entirely subcortical** (thalamus, pulvinar, SC, habenula, ZI), passing through only a thin S1/PCC cortical rind near the top. This matters for interpreting the "why is responsiveness so widespread and depth-independent" caveat — a clean cortical laminar latency gradient was never really the right null expectation for most of this probe's length, since most of it is likely in thalamus/SC/pulvinar, not cortex.

**Limitation (checked, not used further): absolute per-channel depth-in-brain could not be reliably reconstructed.** The elab log has fields for "guide tube final depth," "probe tip enter guide tube depth," and "signal onset depth" (the last being the actual dura/brain-surface crossing — the one reference point that would let electrode-Y-position-in-µm be converted to true mm-below-cortex). "Signal onset depth" is blank for nearly every one of these 12 days, and only 2 (20260324, 20260325) have an explicit final recording depth at all (encoded directly in the channel-selection string, e.g. "Ch2H2@53mm"). So the burr-hole/structure list above is used only as qualitative context (which region this session was *probably* in), not as a per-channel anatomical registration — `sessionAnatomyInfo.m` (used by the new figures, see §5.2) documents this limitation in its header. If you have the actual per-session final depths from elsewhere (e.g. your own notes or the Nandrive logs), that would let this be tightened considerably.

## 3. Pipeline architecture

| File | Purpose |
|---|---|
| `findSaccadeRuns.m` | Scans `Data_collection/m032/*/run-*` for `stim_ct_classic_8targ.m`, returns a struct array of run metadata. |
| `matchRunRecording.m` | Maps a behavioral run to its OpenEphys `experimentN/recordingM` folder. |
| `matchTTLtoLogEvents.m` | Matches TTL digital-line numbers to behavioral event names by rising-edge count. |
| `extractSaccadeSession.m` | Core extraction: MUA + eye channels, trial table, eye calibration. One `.mat` output per session. |
| `runAllExtraction.m` | Batch driver over all 12 runs, tolerant of per-run failures. |
| `loadAllSessions.m` | Loads all extracted `.mat` files into one struct array. |
| `screenResponsiveChannels.m` | Per-channel, per-session responsiveness test in 4 task epochs vs. baseline. |
| `computeDirectionTuning.m` | 8-direction tuning curve, preferred direction, permutation significance, for responsive channels. |
| `analyzeEyeData.m` | Eye calibration application, saccade detection, fixation/pupil/main-sequence analysis. |
| `checkArtifactCorrelation.m` | Trial-residual pairwise channel correlation vs. distance, to distinguish shared-artifact from local-neural signal (§6). |
| `sessionAnatomyInfo.m` | Day → burr hole / channel-selection scheme / planned target structures lookup (§2.2), from `Kwibus_elab.pdf` + `Chamber2Holes.docx`. |
| `forceLightTheme.m` | Forces white background / black text on QC figures regardless of the MATLAB app's light/dark theme setting (fixes illegible dark-mode figures). |
| `classifyTuning.m` | Combines the vector-sum permutation test with a Kruskal-Wallis omnibus test into 3 tuning categories (§4.7): directional / complex / untuned. |
| `analyzePlanningExecution.m` | Splits go-cue-relative and saccade-onset-locked windows into planning vs execution epochs, separately screened and tuned (§4.8). |
| `summarizePlanningExecution.m` | Quantifies the planning/execution split across all 12 sessions (§6.6): recruitment fractions, preferred-direction consistency. |
| `checkSharedComponentKinematics.m` | Correlates the dominant shared trial-residual component (PC1) against saccade kinematics, to test the §6.4 artifact hypothesis directly (§6.5). |
| `runFullPipeline.m` | Runs all of the above end to end. |

## 4. Key methodological decisions (and why)

### 4.1 TTL bit codes are matched empirically, not from `par_default.m`

`par_default.m` defines `Par.StimB=1, Par.TargetB=2, Par.RewardB=3, Par.MicroB=6, Par.CorrectB=7`, with a comment that `SaccadeB=4` and `TrialB=5` are handled by a separate "DasControl" system. **None of this matches what's actually on the wire.** The NI-DAQ digital port in every session examined only ever carries lines `{1, 4, 6, 8}` (occasionally `10`), and cross-checking rising-edge counts against `Log.events` type counts shows:

- line 4 = `trial_start`
- line 6 = `targ_start` (the go-cue / saccade instruction)
- line 8 = `correct` / `reward_start` / `reward_stop` (these three are logged the same number of times whenever every trial that reaches the target epoch is rewarded, so they can't be told apart by count alone — harmless, since the pipeline only needs `targ_start` and `trial_start`)
- line 1 = fires more often than any single logged event type; likely a per-attempt fixation/stim marker that doesn't have a 1:1 logged counterpart.

`matchTTLtoLogEvents.m` does this count-matching **per run**, rather than hard-coding line numbers, specifically because bit assignment is apparently not guaranteed stable across the hardware/software configuration. This paid off: on two sessions (20260317, 20260318), `targ_start`'s count tied exactly with the `correct`/`reward` cluster (all logged 85–88 times that day, since every trial that reached the target was correct), making the match ambiguous by count alone. This was resolved with a secondary heuristic in `extractSaccadeSession.m`: among tied candidate lines, pick whichever has the smaller median latency after the (unambiguous) `trial_start` line — `targ_start` reliably follows `trial_start` by ~0.6s, while `correct`/`reward` follow much later, after the saccade. This tie-break was verified to pick the same line (6) that direct trial-by-trial data matching confirms.

**Checked (2026-07-08): TTL assignment is in fact constant across all 12 sessions, not just per-run-robust.** Reading the raw digital-line inventory directly (before any count-matching) for every session shows the NI-DAQ port carries exactly lines `{1, 4, 6, 8}` in all 12 sessions, no exceptions, with line 4 = `trial_start` and line 6 = `targ_start` in every single one — the two lines the pipeline actually uses never vary. The occasional extra map entries (e.g. `delay_start`/`iti_start`/`iti_stop` shown as "line 4" or "line 8" in some sessions but not others) are incidental count-coincidences for event types the pipeline doesn't use, not evidence of real wiring drift — e.g. `delay_start` and `trial_start` just happen to be logged the same number of times on some days. Per-run empirical matching was kept anyway (cheap, and correctly resolves the 20260317/20260318 tie above), but there's no evidence line assignment ever actually changes across sessions.

### 4.2 Run-to-recording matching, and a within-recording windowing problem

`matchRunRecording.m` matches a run to a recording by **closest absolute time** between the run's folder timestamp and the recording's `sync_messages.txt` "Software Time," within a 90s tolerance — not "closest *preceding*" as originally implemented, because operator workflow was not consistent: usually the OE recording starts a couple of seconds *before* the run script (as expected), but on at least one day the recording started **after** the run folder was created instead (confirmed: 20260303/run-003's matching recording, `experiment2/recording3`, starts 6.3s after its run folder timestamp). A pure "must precede" rule silently matched the wrong (much earlier) recording in that case.

**Correction (2026-07-08): the claim that a single recording can span multiple behavioral runs was checked directly and is not supported by the data.** An earlier version of this report asserted `20260326/experiment1/recording1` covers runs 1–3 as a confirmed example. Re-checking this directly (every run folder on that day vs. every recording's actual start time) shows no such thing: on 20260326, runs 1–4 map to `experiment2/recording1` through `recording4` **one-to-one**, each within ~1–19s of its own run's folder timestamp — there is no shared recording. Extending the same check to **every run on all 12 saccade-direction days** (not just run-003) confirms this holds everywhere: every behavioral run that has a match gets its own dedicated recording, never one shared with a neighboring run. The original claim appears to have been a mistaken or misremembered example and has been retracted. This also means the `20260224/run-001` "no matching recording" conclusion mentioned in an earlier version of this section was itself wrong for an unrelated reason — see the §2.1 correction; it does have a match, just not one the original (non-recursive) search could find.

Given the above, `extractSaccadeSession.m` still **time-windows the TTL stream to the run's own approximate span** before count-matching (the run's start time plus its logged duration from `Log.events`, with generous ±20/+30s margins) purely as a cheap defensive precaution against a multi-run recording occurring on some day not yet checked — not because it was needed to fix any of the 12 sessions actually processed. (The 20260317/20260318 TTL ambiguity was a separate issue — a genuine count-tie between `targ_start` and the `correct`/`reward` cluster, resolved by the latency tie-break in §4.1 — not a multi-run windowing problem.)

### 4.3 Eye channels: empirical identification, not documented anywhere

There is no config file, log field, or channel-name metadata anywhere in the rig, log `Par` struct, or `structure.oebin` that says which of the 8 NI-DAQ analog channels (AI0–AI7) is eye X, eye Y, or pupil. This was inferred from the raw signal shape:

- **AI0, AI1**: clean, correlated step-changes of ~1–4V, time-locked to the go-cue, returning to baseline ~300–450ms later — the signature of a saccade-and-return. Assumed **Eye X, Eye Y**.
- **AI2, AI3**: sit near the ADC negative rail (−10 to −5V), with the raw value pinned at exactly −10.000V for extended stretches (classic blink artifact where the pupil tracker loses the pupil). Assumed **pupil-related**; AI2 is used as "pupil" in the pipeline, AI3 is unused.
- **AI4–AI7**: near-zero, small (~0.05V) crosstalk-looking modulation. Assumed unused/floating.

**This is a data-driven guess, not a confirmed fact, and it is flagged in every extracted session as `S.eyeChannelAssumption`.** It is, however, well validated indirectly: the eye-position calibration (below) achieves R² of 0.75–0.97 against the known 8 target locations, which would not happen if AI0/AI1 weren't genuinely eye position.

### 4.4 Eye calibration: derived per-session from task structure, because no calibration file exists

`Par.SCx/SCy/OFFx/OFFy` in the logs are almost certainly a live GUI display scale, not a true volts→degrees calibration (the values are identical to `par_default.m`'s hardcoded defaults in every session checked, i.e. never actually updated by a calibration routine). Instead, calibration is derived directly from task structure: since the target's position is known exactly on every trial, and on **correct** trials the animal's gaze is by definition on the target, the mean raw voltage during a post-saccade landing window is regressed against the known target (x,y) using a full 2-D affine fit (not independent per-axis gain):

```
[dva_x, dva_y] = [V_x, V_y, 1] · M
```

A full affine (rather than two independent 1-D linear fits) was necessary because the raw channels show real cross-talk — a pure-vertical target measurably displaces the "X" channel and vice versa, consistent with a rotated/skewed camera axis in the analog eye tracker. This raised R² substantially (e.g. one session went from R²=0.18/0.56 with independent per-axis fits to R²=0.76/0.82 with the affine fit).

**The landing window itself required empirical discovery and was a real bug along the way**: this task's dwell on target is brief. Checking raw traces directly (mean X position for rightward vs. leftward target trials, at 20ms resolution) showed gaze reaches the target by ~100–150ms post-go-cue, stays until ~350ms, then a return saccade begins — round-trip complete by ~450–500ms. An initial landing window of 250–550ms (and a later, worse, 350–750ms) partially or mostly captured the *return* saccade instead of the landing fixation, degrading calibration substantially. The final window is **200–350ms** post-go-cue, chosen from directly inspecting the rise/plateau/fall of raw traces (documented in code comments).

**Second real bug, found 2026-07-08 from your question about trajectory plots showing "correct" saccades landing at the wrong location**: blinks corrupt the raw eye **position** channels (X/Y), not just the pupil channel. During a blink, X/Y pin near a **fixed hardware rail of -5.00V** — confirmed empirically identical across all 12 sessions — for anywhere from ~65ms to ~180ms at a time. This was never masked (only the pupil channel had blink-masking, via `opts.blinkThresh`), so:
- Any trial whose 200–350ms landing window overlapped a blink got a garbage landing position averaged into the calibration fit.
- The same contamination fed straight into saccade-onset/velocity detection and the trajectory plots, producing the huge spurious "saccades" to nonsense coordinates you noticed (verified directly: trial 22 of 20260323 showed a real, on-target landing at (7,6.5) at 150–250ms, then a single ~2-frame excursion to (30,24) — a value that recurs near-identically across unrelated trials/targets, confirming it's a fixed rail artifact, not real gaze).

**First fix attempt had its own bug, caught by your follow-up question ("could be temporal clipping? saccade latency is also larger in these runs?") about why upward/upper-diagonal saccades looked collapsed to a few degrees specifically in the earlier sessions.** The first fix flagged blinks with an *absolute* threshold (`abs(raw)>3V`), reasoning that real saccades stay within ~1.5V of baseline. That reasoning about excursion *size* was right, but applying it as an absolute cutoff wasn't: this rig has no fixed zero point — raw baseline varies from about -2.3V to +0.8V across the 12 sessions. In 20260224 specifically (baseline Y ≈ -2.3V), a genuine upward saccade only had to cross ~0.7V further to trip an absolute |raw|>3V test — confirmed directly by tracing raw Y during several "short" upward trials: a clean, physiological rise from baseline to the target with no blink-like plateau and no coincident pupil change, exactly what a real saccade looks like, masked out anyway. This collapsed measured amplitude specifically for whichever directions happened to push that session's baseline past the threshold (here: upward/upper-diagonal, since 20260224's Y baseline was already negative) — it wasn't a general "detection window" or latency-related clipping problem, it was this one session-relative-vs-absolute threshold bug. **Fixed properly**: detect proximity to the actual, fixed rail voltage (-5.00V ± 0.5V tolerance) instead of a generic magnitude — this generalizes correctly across every session's own baseline, since the rail doesn't move with it.

`extractSaccadeSession.m` computes `S.eyeBlinkMask` this way and excludes these samples from the calibration fit; `analyzeEyeData.m` interpolates over them before the zero-phase velocity filter (which cannot tolerate NaNs) and then re-masks them (±10ms padding) before anything is plotted or used for onset/kinematics detection. All 12 sessions were patched twice (`patchEyeBlinkAndCalib.m` — cheap, only touches eye channels, no need to redo the raw AP extraction) and eye analysis re-run after each fix. Net effect of both fixes together: calibration R² improved in every session that had meaningful baseline offset or blink contamination (e.g. 20260224 R²x/y 0.47/0.77 → 0.84/0.84; 20260303 0.97/0.96 → 0.97/0.96 recovered after dipping mid-fix; table below), saccade detection rate recovered to ~99% overall (were briefly as low as 63/85 in 20260224 mid-fix, now 84/85), and per-direction saccade amplitude is now physiological (11-17dva) in every direction in every session — 20260224's 90° amplitude specifically went from an apparent 2.7dva (artifact) to 13.6dva (real) once the rail-based mask replaced the absolute-threshold one.

Some **moderate** (3–7dva) landing scatter remains on plenty of correct trials even after this fix — that's expected, not a bug: `Stm.TargWinSz` = 7dva, i.e. the task's own online correct/error judgment accepts any landing within 7dva of the target center, and real saccades commonly under/overshoot by a few degrees before a small corrective adjustment. A "correct" trial landing several degrees off dead-center is the task's tolerance window working as designed, not a labeling error — the trajectory plots (`*_eye_summary.png`) now draw a light-grey dotted circle of this radius around each target so this is visually obvious rather than needing to be taken on faith.

### 4.5 MUA extraction

Reuses the exact filter chain from `XinyuScripts/getmua.m` (the lab's existing, working pipeline for this same probe/task family): 300Hz high-pass → 12-way ADC-multiplexing phase correction → 5000Hz low-pass → spatial "destripe" high-pass across channels → rectify + 200Hz low-pass envelope → decimate to ~1kHz. Channels are reordered by depth using `NP_PROBE.ELECTRODE_YPOS` from `settings.xml`, same as the existing scripts. Trials are windowed 1.0s before to 1.2s after the go-cue (`targ_start`), long enough to cover the full fixation→stimulus→delay→saccade→return sequence.

### 4.6 Responsiveness screening

Per channel, per session: baseline-subtracted mean MUA in 4 windows (stim, delay, peri-saccadic [0, 0.25s], post-saccadic [0.3, 0.6s]) vs. a pre-stimulus baseline [−0.9, −0.6s], tested with a Wilcoxon signed-rank test across correct trials, Benjamini-Hochberg FDR-corrected across channels per window, plus a minimum effect-size requirement (|mean diff| > 0.3 baseline SD) to avoid flagging trivially small but "significant" effects given large trial counts. A channel is "responsive" if any window passes both criteria.

### 4.7 Direction tuning

For responsive channels only: mean baseline-subtracted response in the peri-saccadic window [0, 0.3s] (relative to go-cue — see §4.8 for why this conflates planning and execution, and the separate analysis that splits them), grouped by the 8 target directions (only correct trials). Preferred direction and tuning strength via vector sum (resultant length, using rectified responses as weights); "tuned" significance via a **permutation test** (1000 shuffles of direction labels, comparing observed resultant length to the shuffled null).

**What "tuned" means, and its blind spot (raised in review, fixed 2026-07-08).** The vector-sum/permutation test above is only sensitive to a **unimodal** direction effect — one clear preferred direction. A channel that responds to two roughly opposite directions (bimodal), or broadly across several adjacent directions without a single dominant peak, can have a real, statistically robust direction effect and still fail this test, because opposing vector contributions cancel each other out in the sum. `classifyTuning.m` (new) addresses this by also running a **Kruskal-Wallis omnibus test** across the 8 direction groups (does response differ across directions *at all*, regardless of shape) and combining the two into three categories, now used for all tuning figures (`tuning_examples.png`, `tuning_matrix.png`, and the new planning/execution curve figures, §4.8):

- **`directional`** (red): vector-sum permutation test significant — has a genuine single preferred direction; `prefDir`/`resultantLen` are meaningful summaries.
- **`complex`** (orange): Kruskal-Wallis significant but the vector-sum test is not — real direction-dependent structure exists, but not a single-peaked one (visibly bimodal or multi-peaked in the small-multiple curve figures); `prefDir`/`resultantLen` are **not** a meaningful summary for these, look at the actual per-direction response instead.
- **`untuned`** (gray): neither test significant.

In `tuning_results.mat`, `tuning(i).tuned` (boolean, `directional`-only) is kept for backward compatibility; `tuning(i).category` and `tuning(i).pKW` are new fields carrying the fuller picture.

### 4.8 New: dissociating planning (premotor) from execution (perimovement) tuning

Raised in review: the [0, 0.3s] go-cue-relative window above mixes two different things whenever saccade latency is non-negligible (median 68-218ms across sessions, §5.1) — **planning/premotor** activity (from the go-cue instruction until the eye actually starts moving) and **execution/perimovement** activity (during and immediately after the saccade itself). `analyzePlanningExecution.m` (new) splits this two ways, producing separate figures without touching the main pipeline's existing outputs:

- **`gocue` mode** — simple fixed split relative to the go-cue, same cutoff for every trial: planning=[0, 0.1]s, execution=[0.1, 0.3]s. Cheap, but the cutoff means something different in every session since latency varies so much (in a 68ms-median-latency session the animal is often already moving well before 0.1s; in a 200ms-median-latency session it's still planning well past 0.1s).
- **`saclocked` mode** (the methodologically correct one) — anchored to each **trial's own detected saccade onset** (from `analyzeEyeData.m`): planning=[-0.15, -0.02]s and execution=[-0.02, +0.15]s relative to onset, not go-cue. Trials with no detected saccade (~1-11% depending on session) are excluded from this mode only.

For each mode: responsiveness (signed-rank + BH-FDR + minimum effect size, matching §4.6) and direction tuning (§4.7's three-category scheme) are computed **separately** for the planning and execution windows. Outputs per session, both modes: `*_planexec_<mode>.png` (per-channel planning-vs-execution effect map, tuning-strength scatter, preferred-direction scatter) and `*_planexec_<mode>_curves.png` (every responsive channel's actual tuning curve, planning on top/execution below, category-colored — the most direct way to see whether/how a channel's tuning shape changes between the two epochs). See `planning_execution_results.mat`.

**First look (saccade-locked mode, run on all 12 sessions): a substantial fraction of channels show little or no direction tuning during planning but become clearly `directional` during execution** — visible directly in the curve figures as a gray/flat top row paired with a red/peaked bottom row for the same channel, and in the `*_planexec_saclocked.png` tuning-strength scatter as most points sitting above the diagonal (execution stronger than planning). This is consistent with genuine perimovement (execution-locked) direction selectivity rather than a sustained premotor direction signal. **For channels that ARE tuned in both epochs, the preferred direction itself is highly consistent between planning and execution** (e.g. 20260312: 66 channels tuned in both, tightly clustered near the diagonal in the preferred-direction scatter, not scattered) — i.e. where a premotor direction signal does exist, it appears to already encode the same direction the execution burst will show, just weaker, rather than encoding something different that later changes. Worth a closer quantitative look across all 12 sessions (e.g. what fraction of `directional`-in-execution channels were `untuned` in planning, and whether the planning-vs-execution preferred-direction consistency holds everywhere) if this distinction matters for your interpretation.

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
| **Total** | **1213 (1192)** | **3182/4608 (69%)** | **832** | | **1202/1213 (99%)** | |

Tuned counts differ by a few from any earlier version of this table because `computeDirectionTuning.m`'s permutation test (1000 random shuffles, no fixed seed) gives slightly different results each run for channels sitting right at the significance boundary — expected run-to-run noise, not a bug.

**Correction (2026-07-08): eye calibration R² and saccade-detection numbers reflect the FINAL blink-masking fix (§4.4, both rounds).** All 12 sessions now have R² ≥ 0.84 on both axes and ≥98% saccade detection except 20260224 (84/85, still the weakest session overall on multiple measures — see §6.4). An intermediate version of this table (now superseded) briefly showed much worse numbers for several sessions after the *first* blink-masking attempt turned out to have its own bug (an absolute rather than session-relative voltage threshold, §4.4) — that intermediate state was never a real regression in the underlying data, just a bug in a bug-fix that's now corrected properly.

**Correction (2026-07-08):** `20260224/run-001` was previously excluded as having no matching recording — wrong, see §2.1. It's now included above (4/64 responsive channels tuned — the lowest tuned fraction and, notably, the worst eye-calibration R²x=0.47 of any session, both worth keeping in mind if this session gets used for anything beyond the population count).

**Correction (2026-07-08):** `extraction_summary.mat` previously marked 20260317 and 20260318 as `FAILED` ("could not identify targ_start/trial_start TTL lines") — stale, from before the tie-break fix in §4.1 was added. Both extracted successfully (88 and 85 trials respectively, included in the table above) and this is reflected in the per-session table here. `extraction_summary.mat` on the server was rebuilt from the actual extracted files (`patchExistingExtractions.m`) and now shows `ok` for all 12.

Note the clear **latency drop** from ~170–195ms in the earlier sessions to ~70–120ms from 20260320 onward — plausibly a training/practice effect (faster, more automatic responses later in the recording series), worth checking against the session order/history in `Kwibus_elab.pdf`.

Saccade detection (velocity-threshold on a low-pass-filtered, calibrated eye trace, blink periods excluded — §4.4) succeeded on 99% of trials (1202/1213); endpoint accuracy and main-sequence (peak velocity vs. amplitude) relationships were physiological in every session inspected, with no more of the wild outlier "saccades" seen before the blink fix (see `*_eye_summary.png`).

### 5.2 QC figures (per session, in `data/figures/`)

**Reworked 2026-07-08** based on review feedback. All figures now force a white background/black text regardless of the MATLAB app's light/dark theme (`forceLightTheme.m`) — the previous versions were illegible on a dark-themed installation. A first pass at this missed several object types that carry their own color independent of the axes: **colorbars** (tick labels/axis line have their own `Color` property, untouched by an axes/text sweep), **tiledlayout titles** (`title(tl,...)` sets a property on the layout object itself, not an axes `Title`), and a **line explicitly colored white** (the go-cue marker on the overview figure, originally chosen to be visible on a dark background — literally invisible once the background is white). `forceLightTheme.m` now explicitly handles colorbars and tiledlayout titles/subtitles, and auto-repaints any leftover pure-white lines to black; the overview figure's anatomy text panel also now word-wraps instead of clipping past its tile edge.

- **`*_overview.png`** — four panels, all sharing one channel-index y-axis with **channel 1 (deepest, nearest probe tip) at the bottom and the shallowest channel at the top**, matching anatomical up/down: (1) a real-depth ruler (electrode position in mm from the probe tip, so gaps/clusters from the day's online channel-selection scheme are visible directly, not just implied by channel index), (2) the z-scored MUA heatmap with all 5 trial-epoch boundaries (fixation/stimulus/delay/go-cue/target-end) drawn as thick colored lines with a legend, computed per-session from that session's own `Stm` timing rather than hardcoded, (3) a per-channel, per-epoch significance map (which of stim/delay/peri/post each channel is significant in, signed by effect direction, gray where not significant) so you can see *which* window a channel's "responsive" flag refers to instead of just a binary map, (4) a text panel with that day's burr hole, planned target structures, and online channel-selection scheme (§2.2).
- **`*_eye_summary.png`** — 7 panels: trajectory spoke plots shown **twice**, once for all trials (error trials dimmed) and once for correct trials only with **% correct in the title**; both now draw a **light-grey dotted circle** (radius = `Stm.TargWinSz`, the task's real hit-detection window) around each target, behind the traces, so a saccade landing off dead-center but inside the circle reads as the legitimate hit it is; fixation scatter; main sequence; latency histogram **in ms**; pupil trace; accuracy by direction. Direction colors use a fixed hand-picked palette, not `hsv(8)` — the latter put a pale yellow-green at the 90° slot that nearly vanished at low alpha on white, which looked like "fewer trials" for that direction but was a color-visibility artifact (verified against the actual trial counts).
- **`*_tuning_examples.png`** — **all responsive channels** (not a capped/arbitrary subset — the old version showed `min(nTuned,12)`, which is why the count looked inconsistent across sessions), sorted by depth, one small tuning-curve tile per channel (min-max normalized per channel for display). Border **and line color** now reflect the 3-way category from §4.7: red=`directional`, orange=`complex`, gray=`untuned` (previously just a binary red/gray border).
- **`*_tuning_matrix.png`** — all responsive channels, sorted by preferred direction, with a **bounded per-channel min-max normalization** instead of dividing by the row max (the old version could blow up the shared color axis from a single noisy, weakly-modulated channel and wash out every well-behaved row — see code comment in `computeDirectionTuning.m`), plus a tuning-strength color strip and a 3-category color strip (§4.7) alongside so amplitude and tuning shape aren't lost to the normalization.
- **`*_population_summary.png`** — **now generated once per session, not pooled across all 12** (each session is an independent penetration, often a different burr hole/target structure per §2.2, and the probe isn't seated to the same depth across sessions, so cross-session pooling implied false comparability). Preferred-direction rose plot and tuning-strength histogram for that session, plus preferred direction vs. **real electrode depth (mm from tip, this session only)** — see §6.4, this reveals a striking within-session depth relationship in at least one session that the old pooled/index-based version could never have shown.
- **`*_artifact_check.png`** — pairwise trial-residual correlation matrix across all channels (direction-tuning contribution removed), shown raw and sorted by real depth, plus mean correlation vs. inter-channel distance — see §6.4.
- **`*_planexec_gocue.png` / `*_planexec_saclocked.png`** (new, §4.8) — per-channel effect map (planning vs execution), tuning-strength scatter, and preferred-direction scatter, for the go-cue-fixed and saccade-onset-locked window splits respectively.
- **`*_planexec_gocue_curves.png` / `*_planexec_saclocked_curves.png`** (new, §4.8) — every responsive channel's actual tuning curve, planning epoch on top and execution epoch below, category-colored (§4.7) — the most direct visual answer to "is this channel's direction tuning different between planning and execution."
- **`*_shared_kinematics.png`** (new, §6.5) — per-channel PC1 loadings (the dominant shared trial-residual component) and its per-trial score plotted against saccade peak velocity and amplitude.

## 6. Important caveats — read before using the tuning results

**Updated substantially 2026-07-08** following review: most recording sites here are subcortical (§2.2), and a new per-channel trial-residual correlation check plus per-session (not pooled) population summaries provide much stronger evidence than before that at least some of this is genuine, spatially-organized neural signal — not purely a shared artifact. The picture is now genuinely mixed rather than "probably all artifact." Details below.

### 6.1 The high responsiveness fraction is now better explained by online channel selection + subcortical targets than by assuming a uniform artifact

74% of all channels across all sessions were flagged "responsive," reaching 98% in one session (20260323). Two things temper how alarming this is, both raised in review and confirmed by checking the source records:

- **Online channel selection is not per-electrode, it's per bank/section** (your point, confirmed directly): `Kwibus_elab.pdf` logs an explicit operator-chosen scheme for every session — "quarter density," "4bankeven," "single column," "sections2/3" (§2.2 table). The operator was already restricting acquisition toward electrode banks that looked responsive online, so a high responsive fraction *within the selected set* is an expected consequence of that selection, not evidence of an indiscriminate artifact.
- **Most of these trajectories are subcortical at depth** (§2.2): thalamus (MD/CL/CM), pulvinar, SC, habenula, ZI — not cortex. A depth-dependent *laminar* latency gradient (the original null expectation in this section) is a cortical-columnar signature; it was never really the right thing to expect from a probe that's mostly in thalamus/SC by design. Absence of a laminar gradient is therefore much weaker evidence of artifact than originally stated.

This softens, but does not eliminate, the concern — see §6.4 for the direct check.

### 6.2 Preferred-direction clustering across sessions is less suspicious once you account for same-hemisphere sampling

All 12 sessions are the same hemisphere (implicit in a single chamber), so a consistent contralateral-field bias in preferred direction across sessions is exactly what a real retinotopic/motor map (SC, pulvinar) would produce — it is not on its own evidence of a shared artifact the way it would be if sessions sampled both hemispheres or unrelated brain regions. This doesn't fully explain the tight ~90–150° clustering, but it removes "different unrelated sites shouldn't agree" as a strong argument, since these sites are not anatomically unrelated in the way independent penetrations in different regions/hemispheres would be.

### 6.3 Cross-session depth comparison — fixed

Per your point: sessions are independent penetrations, the probe is not seated to the same depth across sessions, and (per §2.2) different sessions often target different burr holes/structures entirely. The population summary is now generated **per session** (`*_population_summary.png`, one per session, not one pooled figure), using **real electrode depth in mm from the probe tip** rather than an abstract cross-session channel index. Every screening/tuning computation was already strictly within-session (`tuning_results.mat` keyed by `Day`/`RunN`/`channel`) — only the pooled-plot presentation was fixed.

### 6.4 New: per-channel trial-residual correlation check (the direct artifact test)

`checkArtifactCorrelation.m` (new) tests the shared-artifact hypothesis directly: for each channel pair, correlate the **trial-to-trial residual** peri-saccadic response (each direction-condition's own mean subtracted out first, so genuine shared direction tuning can't inflate this) and look at how that correlation depends on physical inter-channel distance. A uniform artifact riding on every channel predicts high correlation that does **not** decay with distance; local neural activity predicts correlation that decays with distance.

Mean pairwise correlation per session (`*_artifact_check.png`):

| Day | Mean corr | Day | Mean corr | Day | Mean corr |
|---|---|---|---|---|---|
| 20260224 | 0.078 | 20260312 | 0.213 | 20260323 | 0.235 |
| 20260226 | 0.112 | 20260317 | 0.078 | 20260324 | 0.100 |
| 20260303 | 0.203 | 20260318 | 0.179 | 20260325 | 0.057 |
| 20260310 | 0.099 | 20260320 | 0.077 | 20260326 | 0.021 |

Two findings, looking at `*_artifact_check.png` and `*_population_summary.png` together:

1. **Correlation decays sharply with distance in every session** (e.g. 20260323: ~0.55 at <200µm down to ~0.18 by 1mm+), which is the signature of local neural activity, not a probe-wide shared signal. But it does not decay all the way to zero — it plateaus at a small positive value (~0.1–0.2 in the higher-correlation sessions) even between far-apart channels, consistent with a **smaller residual shared component riding underneath genuine local structure** (common reference, EMG, or motion), not a clean "all real" result either.
2. **The per-session population summaries show genuine, session-specific spatial organization that a shared artifact would not produce**: 20260312 (burr hole 1, SC) shows a strikingly clean, repeating sawtooth relationship between real depth and preferred direction across ~4 distinct depth clusters; 20260323 shows a smooth, continuous, near-monotonic depth-direction relationship; 20260320 (burr hole 4, thalamus) shows depth-banded clusters of consistent preferred direction. By contrast, **20260326 — the session with the lowest mean pairwise correlation (0.021) and lowest tuned fraction (8/156)** — shows no depth structure at all, just scattered noise. This session-to-session variation (highly organized vs. apparently unstructured, tracking the correlation number) is itself an argument against a single fixed artifact explaining everything: a universal shared-artifact source should contaminate every session similarly, not appear strongly organized in some and absent in others.

**Net assessment**: this is a mixed result, not a clean resolution. Several sessions (20260312, 20260320, 20260323 especially) show good evidence of genuine, spatially-organized tuning. A smaller shared/common-mode component likely coexists across all sessions (the nonzero correlation plateau at large distances). 20260326 looks like it may have little real signal at all. Recommend treating tuning results session-by-session rather than as one uniform "832 tuned channels" — some sessions look considerably more trustworthy than others by this metric.

### 6.5 New: does the shared component track saccade kinematics? (mixed answer, not resolved)

Direct follow-up to §6.4's open question. `checkSharedComponentKinematics.m` (new) extracts the dominant shared fluctuation directly — the first principal component (SVD) of the same trial-residual matrix used in §6.4 — and correlates its per-trial score against saccade peak velocity and amplitude (Spearman). A real, consistent correlation with both kinematic measures would point to EMG/movement contamination; no correlation would argue against that specific explanation (though not rule out a different shared source, e.g. a reference issue).

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

Findings:
1. **PC1 (the dominant shared fluctuation) explains a lot of the trial-to-trial residual variance in most sessions** — from 10% up to 74% (20260312) — and loading sign agreement across channels is high (79-99%) in every session except **20260326 (55%, essentially chance)**, meaning that session's residual variance isn't dominated by one common-mode axis at all — consistent with §6.4's finding that 20260326 shows no organized structure of any kind.
2. **The kinematic correlation is real but not universal: exactly 4/12 sessions** (20260310, 20260320, 20260323, 20260324) show a **significant, same-sign correlation with BOTH peak velocity and amplitude** — the clean signature of EMG/movement contamination this check was designed to find. The other 8 sessions show no reliable relationship (one, 20260312, shows a significant correlation with peak velocity only, not amplitude — weaker/inconsistent evidence).

**Interpretation**: the shared component's origin is not uniform across the dataset. In at least 4 sessions there's direct, quantitative evidence it's genuinely kinematics-linked (EMG/motion-like). In the rest, either the shared component has a different origin (reference/hardware), or the correlation approach lacks power in sessions with less kinematic variance/trial count — this isn't distinguishable from the current analysis. Doesn't change the session-by-session recommendation above.

### 6.6 New: planning-vs-execution tuning, quantified across all 12 sessions

§4.8's "first look" was 1-3 sessions; `summarizePlanningExecution.m` (new) quantifies the saccade-onset-locked split across all 12. Two sessions (20260224, 20260326 — already flagged as the weakest/noisiest by every other metric in this report) had **zero or one** directional channels in either window and are excluded from the pooled stats below.

- **Pooled across the other 10 sessions: 835 channel-instances were `directional` in execution; only 47% of those were also `directional` during planning — 53% (439) only became directional during execution**, i.e. genuine recruitment of additional direction-selective activity specifically around the movement, not just an amplitude increase on an already-tuned population.
- Conversely, of 606 `directional`-in-planning channel-instances, 65% remained `directional` into execution (35%, 210, were planning-only — tuned before the movement but not during/after it).
- **For channels tuned in both epochs, preferred direction is meaningfully consistent, not independent**: mean absolute angular difference (planning vs execution) pooled = **40.5°**, well below the ~90° expected for unrelated/random pairings, though far from 0° — i.e. where a premotor signal exists, it tends to already point toward the direction the execution burst will show, with some genuine drift/refinement, not a totally different or unrelated direction. Per-session values range from 14.5° (20260325, very consistent) to ~50-80° (20260226, 20260303, 20260317 — closer to random, though these sessions also have the fewest both-tuned channels, so noisier estimates).

## 7. Open items requiring your input

1. **AI channel mapping** (§4.3) — per your review, looks correct; no further action needed.
2. **Artifact investigation** (§6.1, §6.4, §6.5) — the shared-component/kinematics check (§6.5) gives a real but non-universal answer: 4/12 sessions show a clean EMG/movement signature, the rest don't. Not fully resolved; if you want to push further, the natural next step is checking whether the 8 non-significant sessions differ systematically (e.g. less kinematic variance, fewer trials, different burr hole/structure) from the 4 significant ones.
3. **Cross-session depth comparison** (§6.3) — fixed: population overviews are now per-session using real depth. Absolute per-channel depth-*in-brain* (mm below cortical surface) is still not reconstructable from the current logs (§2.2) — only relative depth along the probe.
4. **Burr hole for 20260224 is genuinely ambiguous** (marked "5?" in the source elab log itself, §2.2) — noted, not resolved; this session is independently the weakest on nearly every other metric in this report (§5.1, §6.4, §6.5, §6.6), so it's likely not worth chasing further unless it specifically matters for your interpretation.
5. **Planning-vs-execution** (§4.8, §6.6) — now quantified across all 12 sessions; if useful, the natural extension is checking whether the `execution-only` recruitment fraction (53% pooled) varies with burr hole/structure (§2.2) or with the artifact metrics above.
6. **Responsiveness/tuning thresholds**: current thresholds (FDR 0.05, effect size 0.3 baseline SD for screening; permutation p<0.05 for tuning) are reasonable defaults; given §6.4's session-by-session variability, you may want to treat sessions individually rather than adjust one global threshold.

## 8. Recommended next steps

1. **Given §6.4/§6.5's mixed result**, treat sessions individually rather than pooling a single "tuned channel" count — prioritize 20260312/20260320/20260323/20260324 (clean depth-organized structure and/or confirmed kinematic-linked shared component) for any further biological claims, and treat 20260224/20260326 (weakest on nearly every metric) with more skepticism.
2. **Done** (§6.5): shared component correlated against saccade kinematics — real in 4/12 sessions, absent in the rest. If you want to act on this: a common-mode removal step (median/PCA across the array before MUA extraction) could be added to `extractSaccadeSession.m` and re-run, at least for the 4 sessions with a confirmed kinematic link, to see whether the depth-organized tuning in §6.4 survives removing it.
3. Absolute per-channel depth-in-brain (mm below cortical surface, needed to sharpen the §2.2 burr-hole/structure overlay into an actual per-channel anatomical label) needs a numeric anchor not currently in the logs — see §2.2's limitation note. If you have the true final recording depths from elsewhere, `sessionAnatomyInfo.m` is the place to add them.
4. Once sessions are triaged by §6.4/§6.5, revisit direction tuning with a von Mises fit (not just vector sum) for a cleaner tuning-width estimate on the trustworthy sessions — the saccade-onset-locked alignment from §4.8 is already implemented, use it in place of go-cue alignment for this.
5. **Done** (§6.6): planning-vs-execution quantified across all 12 sessions. If useful, check whether the 53% execution-only recruitment fraction varies systematically with burr hole/structure or with the §6.5 kinematic-link sessions.
