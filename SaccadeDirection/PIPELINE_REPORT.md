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
| `runFullPipeline.m` | Runs all of the above end to end. |

## 4. Key methodological decisions (and why)

### 4.1 TTL bit codes are matched empirically, not from `par_default.m`

`par_default.m` defines `Par.StimB=1, Par.TargetB=2, Par.RewardB=3, Par.MicroB=6, Par.CorrectB=7`, with a comment that `SaccadeB=4` and `TrialB=5` are handled by a separate "DasControl" system. **None of this matches what's actually on the wire.** The NI-DAQ digital port in every session examined only ever carries lines `{1, 4, 6, 8}` (occasionally `10`), and cross-checking rising-edge counts against `Log.events` type counts shows:

- line 4 = `trial_start`
- line 6 = `targ_start` (the go-cue / saccade instruction)
- line 8 = `correct` / `reward_start` / `reward_stop` (these three are logged the same number of times whenever every trial that reaches the target epoch is rewarded, so they can't be told apart by count alone — harmless, since the pipeline only needs `targ_start` and `trial_start`)
- line 1 = fires more often than any single logged event type; likely a per-attempt fixation/stim marker that doesn't have a 1:1 logged counterpart.

`matchTTLtoLogEvents.m` does this count-matching **per run**, rather than hard-coding line numbers, specifically because bit assignment is apparently not guaranteed stable across the hardware/software configuration. This paid off: on two sessions (20260317, 20260318), `targ_start`'s count tied exactly with the `correct`/`reward` cluster (all logged 85–88 times that day, since every trial that reached the target was correct), making the match ambiguous by count alone. This was resolved with a secondary heuristic in `extractSaccadeSession.m`: among tied candidate lines, pick whichever has the smaller median latency after the (unambiguous) `trial_start` line — `targ_start` reliably follows `trial_start` by ~0.6s, while `correct`/`reward` follow much later, after the saccade. This tie-break was verified to pick the same line (6) that direct trial-by-trial data matching confirms.

### 4.2 Run-to-recording matching, and a within-recording windowing problem

A single OpenEphys "recording" can span **multiple behavioral runs** if acquisition was never paused between them (confirmed directly for several days, e.g. `20260326/experiment1/recording1` covers runs 1–3). `matchRunRecording.m` matches a run to a recording by **closest absolute time** between the run's folder timestamp and the recording's `sync_messages.txt` "Software Time," within a 90s tolerance — not "closest *preceding*" as originally implemented, because operator workflow was not consistent: usually the OE recording starts a couple of seconds *before* the run script (as expected), but on at least one day (20260303) the recording started **6–11 seconds after** the run folder was created. A pure "must precede" rule silently matched the wrong (much earlier) recording in that case. The 90s tolerance is tight enough that `20260224/run-001` — where the actual first recording of the day starts 72 minutes after the run — is correctly reported as having **no matching neural recording** rather than being forced onto the wrong one.

Because a matched recording can contain other runs' trials too, the TTL stream is **time-windowed to just this run's span** before any count-matching is attempted: the run's approximate start time (converted into the recording's own timestamp clock via the recording's absolute start time) plus the run's own logged duration (from `Log.events`), with generous ±20/+30s margins. Skipping this step was the direct cause of the two failures above resolving incorrectly on the first pass — the unwindowed TTL stream on multi-run recordings contains events from neighboring runs too, and count-matching against a single run's `Log` then fails outright or (worse, though not observed here) silently matches the wrong slice.

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

### 4.5 MUA extraction

Reuses the exact filter chain from `XinyuScripts/getmua.m` (the lab's existing, working pipeline for this same probe/task family): 300Hz high-pass → 12-way ADC-multiplexing phase correction → 5000Hz low-pass → spatial "destripe" high-pass across channels → rectify + 200Hz low-pass envelope → decimate to ~1kHz. Channels are reordered by depth using `NP_PROBE.ELECTRODE_YPOS` from `settings.xml`, same as the existing scripts. Trials are windowed 1.0s before to 1.2s after the go-cue (`targ_start`), long enough to cover the full fixation→stimulus→delay→saccade→return sequence.

### 4.6 Responsiveness screening

Per channel, per session: baseline-subtracted mean MUA in 4 windows (stim, delay, peri-saccadic [0, 0.25s], post-saccadic [0.3, 0.6s]) vs. a pre-stimulus baseline [−0.9, −0.6s], tested with a Wilcoxon signed-rank test across correct trials, Benjamini-Hochberg FDR-corrected across channels per window, plus a minimum effect-size requirement (|mean diff| > 0.3 baseline SD) to avoid flagging trivially small but "significant" effects given large trial counts. A channel is "responsive" if any window passes both criteria.

### 4.7 Direction tuning

For responsive channels only: mean baseline-subtracted response in the peri-saccadic window [0, 0.3s], grouped by the 8 target directions (only correct trials). Preferred direction and tuning strength via vector sum (resultant length, using rectified responses as weights); significance via a **permutation test** (1000 shuffles of direction labels, comparing observed resultant length to the shuffled null).

## 5. Results

### 5.1 Per-session summary

| Day | Trials (correct) | Responsive | Tuned | Eye calib R² (x/y) | Saccade detected | Median latency |
|---|---|---|---|---|---|---|
| 20260226 | 145 (144) | 271/384 | 44 | 0.89 / 0.93 | 144/145 | 195ms |
| 20260303 | 114 (112) | 201/384 | 16 | 0.97 / 0.96 | 113/114 | 179ms |
| 20260310 | 117 (116) | 238/384 | 54 | 0.85 / 0.91 | 116/117 | 179ms |
| 20260312 | 85 (84) | 361/384 | 172 | 0.87 / 0.91 | 84/85 | 176ms |
| 20260317 | 88 (88) | 270/384 | 23 | 0.83 / 0.82 | 87/88 | 171ms |
| 20260318 | 85 (85) | 247/384 | 29 | 0.85 / 0.94 | 85/85 | 165ms |
| 20260320 | 100 (99) | 304/384 | 117 | 0.80 / 0.86 | 98/100 | 100ms |
| 20260323 | 95 (94) | 377/384 | 134 | 0.75 / 0.80 | 93/95 | 117ms |
| 20260324 | 99 (98) | 354/384 | 132 | 0.88 / 0.89 | 99/99 | 68ms |
| 20260325 | 99 (97) | 339/384 | 91 | 0.76 / 0.83 | 99/99 | 69ms |
| 20260326 | 101 (100) | 156/384 | 9 | 0.76 / 0.82 | 101/101 | 79ms |
| **Total** | **1128 (1107)** | **3118/4224 (74%)** | **821** | | **1319/1329 (99%)** | |

`20260224/run-001` was excluded: the closest OpenEphys recording that day starts 72 minutes after the run, i.e. acquisition had not started yet — there is no neural data for that behavioral run.

Note the clear **latency drop** from ~170–195ms in the earlier sessions to ~70–120ms from 20260320 onward — plausibly a training/practice effect (faster, more automatic responses later in the recording series), worth checking against the session order/history in `Kwibus_elab.pdf`.

Saccade detection (velocity-threshold on a low-pass-filtered, calibrated eye trace) succeeded on 99% of trials; endpoint accuracy and main-sequence (peak velocity vs. amplitude) relationships were physiological in every session inspected (see `*_eye_summary.png`).

### 5.2 QC figures (per session, in `data/figures/`)

- `*_overview.png` — full-probe MUA heatmap (z-scored to baseline) and responsive-channel map
- `*_eye_summary.png` — trajectory "spoke" plot, fixation scatter, main sequence, latency histogram, pupil trace, accuracy-by-direction
- `*_tuning_examples.png`, `*_tuning_matrix.png` — polar tuning curves and direction-sorted tuning matrix for tuned channels
- `population_tuning_summary.png` — preferred-direction rose plot, tuning-strength histogram, and preferred-direction-vs-depth scatter, pooled across all 11 sessions

## 6. Important caveats — read before using the tuning results

### 6.1 The high responsiveness fraction likely includes a shared, non-specific signal

74% of all channels across all sessions were flagged "responsive," reaching 98% in one session (20260323). Looking at `20260323_run-003_overview.png` directly: the MUA increase around the go-cue is **near-uniform across the entire depth of the probe**, essentially simultaneous, with no depth-dependent latency gradient. Genuine layer/columnar neural responses to a saccade would be expected to vary in latency and magnitude across depth; a signal that turns on at the same instant on (nearly) every one of 384 channels looks like a shared, common-mode contamination — most plausibly saccade-related EMG (eye muscle activity) or a cable-movement/reference artifact synchronized to the eye movement, riding on top of whatever genuine neural signal is also present. This does not mean the screening or the extraction is broken; it means "responsive" as currently defined cannot yet distinguish genuine local neural modulation from a shared movement artifact.

### 6.2 The population preferred-direction distribution is suspicious in the same way

`population_tuning_summary.png` (left panel) shows preferred directions pooled across all 821 tuned channel-sessions cluster tightly around ~90–150°, with almost none at 0° or 180–330°. If tuning were driven by genuinely diverse neural populations recorded across 11 independent penetrations (different days, and per your note below, likely different burr holes/depths), preferred direction should **not** consistently cluster the same way across unrelated recording sites — different sites should sample different parts of any motor/visual map. A single dominant direction reproduced across independent sessions is much more consistent with a systematic, direction-dependent artifact (e.g. eye-muscle EMG or cable strain that is inherently larger for saccades toward upper-left than other directions) than with 821 genuinely tuned neurons.

**Recommendation before trusting either number**: check whether the same locked, non-specific response is present on any clearly-non-neural or far-out-of-brain channel (e.g. a channel known to be in white matter or above the brain surface for that session, if identifiable from `Chamber2Holes.docx` / histology), and/or check whether trial-by-trial response amplitude in the "peri" window correlates with saccade peak velocity or amplitude in a way that is essentially identical across all channels (a hallmark of shared kinematic/EMG contamination rather than distinct neural tuning).

### 6.3 Your point: sessions are independent penetrations — channels are not comparable across sessions

This is correct and **already respected in the underlying data structure** — every row in `tuning_results.mat` is keyed by `(Day, RunN, channel)`, and every screening/tuning computation is done strictly within one session. No trial data, baselines, or statistics are pooled across sessions.

However, **one figure currently misrepresents this**: `population_tuning_summary.png`'s right-hand panel plots "preferred direction vs. channel depth index," pooling the raw depth-ordered channel index (1–384) across all 11 sessions as if position 200 in one session's probe geometry means the same tissue location as position 200 in another. Since each session is a distinct penetration (different burr hole and/or different insertion depth per `Chamber2Holes.docx`), **this comparison is not meaningful as currently plotted** and should not be used to infer any spatial/laminar organization of preferred direction across the dataset. It should either be dropped, or re-expressed as depth *relative to a session-specific anatomical reference* (e.g. relative to estimated brain entry point, if that can be recovered per session from the chamber/depth log), and even then should be interpreted per-session/per-burr-hole rather than pooled, unless there's reason to believe multiple sessions targeted the same location.

The **left and center panels** (pooled preferred-direction rose plot and tuning-strength histogram) are not subject to this specific problem — they only pool "one number per tuned channel-session" without implying cross-session positional correspondence — but the *interpretation* of the clustered rose plot (6.2 above) still needs to account for each session being an independent sample.

## 7. Open items requiring your input

1. **Confirm the AI channel mapping** (§4.3): AI0=EyeX, AI1=EyeY, AI2=Pupil. This is inferred, not documented, though indirectly well-supported by calibration R².
2. **Artifact investigation** (§6.1–6.2): decide whether to pursue re-referencing, single-channel differential analysis, or correlating response amplitude with saccade kinematics to separate genuine neural signal from shared movement artifact before treating "821 tuned channels" as a biological result.
3. **Cross-session depth comparison** (§6.3): decide whether/how you want per-session depth related to a common anatomical reference (would need the burr-hole/depth estimates in `Chamber2Holes.docx`, not yet used by this pipeline), or whether population summaries should simply stay per-session.
4. **Responsiveness/tuning thresholds**: current thresholds (FDR 0.05, effect size 0.3 baseline SD for screening; permutation p<0.05 for tuning) were chosen to be reasonable defaults, not tuned against the artifact concern above — may want revisiting once §6.1–6.2 are resolved.

## 8. Recommended next steps

1. Investigate the shared-artifact hypothesis (§6.1–6.2) before any further biological interpretation.
2. If confirmed as artifact: consider re-referencing (e.g. median/PCA-based common-noise removal across the array before MUA extraction) and re-running screening/tuning.
3. Incorporate `Chamber2Holes.docx` burr-hole/depth estimates so each session's channels can be related to a consistent anatomical depth axis, enabling meaningful cross-session comparison per your point in §6.3.
4. Once genuine responsive/tuned channels are isolated, revisit direction tuning with a von Mises fit (not just vector sum) for a cleaner tuning-width estimate, and consider aligning to detected saccade onset (from `analyzeEyeData.m`, already computed per trial) rather than the go-cue, since saccade latency varies session-to-session (§5.1) and trial-to-trial.
