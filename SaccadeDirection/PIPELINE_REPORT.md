# Saccade Direction Pipeline — Report

Subject: m032 ("Kwibus"). Task: 8-target visually-guided saccade task, run under the shared curve-tracing ("TDCT"/"ct") behavioral codebase with `Stm.difficulty = 1` (single target, no curve, no distractor) and stimulus file `stim_ct_classic_8targ.m`. Data source: `/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_collection/m032`.

All code lives in `SaccadeDirection/code/`. All extracted data and figures live in `SaccadeDirection/data/` and `SaccadeDirection/data/figures/`.

---

## 1. Objective

Extract raw Neuropixels + eye-tracking data for every saccade-direction session, reduce the probe data to a multi-unit activity (MUA) envelope per channel per trial, screen every channel for task responsiveness, and compute saccade-direction tuning for responsive channels — with the eye-tracking data used both as a behavioral readout in its own right (fixation stability, saccade kinematics, pupil) and as the ground truth for calibrating and QC-ing the neural alignment.

## 2. Identifying the task and the sessions

A run is a saccade-direction run if and only if its folder contains `stim_ct_classic_8targ.m`. Scanning all of `m032`'s day folders for this marker found **12 runs**, one per day, spanning 2026-02-24 to 2026-03-26:

| Day | Run | Day | Run | Day | Run |
|---|---|---|---|---|---|
| 20260224 | run-001 | 20260312 | run-003 | 20260323 | run-003 |
| 20260226 | run-003 | 20260317 | run-003 | 20260324 | run-003 |
| 20260303 | run-003 | 20260318 | run-003 | 20260325 | run-003 |
| 20260310 | run-003 | 20260320 | run-003 | 20260326 | run-003 |

The stimulus file confirms this is a genuine 8-direction task: `Stm.Targ(1..8).xy_dva` places targets at 0°, 45°, 90°, ..., 315° at 12 dva eccentricity, with `Stm.difficulty = 1` disabling the curve/distractor machinery inherited from the curve-tracing paradigm, leaving a plain visually-guided saccade: fixate → brief target flash (100ms) → 400ms memory delay → go-cue → saccade.

Older 2024/2025 sessions in the same folder tree use the same task family (`task-ct`) but only the curve-tracing variants (`fgdots`/`classic` without `_8targ`), and are excluded.

### 2.1 Run → OpenEphys recording → behavioral log mapping

Each behavioral run folder is matched to its OpenEphys `experimentN/recordingM` folder by `matchRunRecording.m` (closest absolute start time between the run's folder timestamp and the recording's `sync_messages.txt` "Software Time," within a 90s tolerance) and to its own behavioral log `.mat` by `findSaccadeRuns.m` (the `sub-*.mat` file inside the run folder). Recording folders are searched recursively under the day folder, since they are not always in the default `experimentN/recordingM` location directly under it — on 20260224, for instance, the recordings are nested inside an additional `m032-2026-02-24_12-59-28/Record Node 110/` folder. The verified pairing for all 12 sessions:

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

Each `_extracted.mat` stores this verified path as `S.LogFile`. **Do not use `S.Par.LogFn`/`S.Par.LogPath` for this purpose** — those fields come from inside the log file's own saved `Par` struct and are stale in every session checked (they point at an unrelated earlier `task-rfmap`/`stim-sparsenoise` run in the same stim-program session, not the 8-target run the data was actually extracted from).

### 2.2 Burr hole, channel selection, and planned target structure per session

`Chamber2Holes.docx` (MRI-based trajectory planning) documents 5 burr holes (0, 1, 2, 4, 5), each with a planned insertion depth and a list of structures its trajectory passes through, superficial to deep. `Kwibus_elab.pdf` (the day-by-day lab notebook) gives the burr hole and online channel-selection scheme actually used on each recording day (e.g. "quarter density," "4bankeven," "single column," "sections2/3" — the operator manually restricting acquisition to part of the probe based on what looked good online; this is a per-bank/section choice, not per-electrode). Cross-referencing the two, all 12 sessions:

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

**At planned full depth, every one of these trajectories is mostly or entirely subcortical** (thalamus, pulvinar, SC, habenula, ZI), passing through only a thin S1/PCC cortical rind near the top — relevant context for §6.1.

**Limitation: absolute per-channel depth-in-brain could not be reliably reconstructed.** The elab log has fields for "guide tube final depth," "probe tip enter guide tube depth," and "signal onset depth" (the last being the actual dura/brain-surface crossing — the reference point needed to convert electrode-Y-position-in-µm to true mm-below-cortex). "Signal onset depth" is blank for nearly every one of these 12 days, and only 2 (20260324, 20260325) have an explicit final recording depth at all (encoded directly in the channel-selection string, e.g. "Ch2H2@53mm"). The burr-hole/structure list above is therefore qualitative context (which region a session was *probably* in), not a per-channel anatomical registration — `sessionAnatomyInfo.m` documents this limitation in its header.

## 3. Pipeline architecture

**Figure organization**: `data/figures/` is organized **by figure type**, not by session — `figures/<type>/<Day>_run-<RunN>.png` (e.g. `figures/overview/20260312_run-003.png`), so the same plot across all 12 sessions is easy to find and compare side by side. `figSavePath.m` centralizes this convention for every per-session plotting function. Cross-session summary figures (no single day/run, e.g. the pooled decoding/area figures) live in `figures/pooled/` instead. The per-channel direction timecourses (§4.11) are organized session-first instead (`figures/timecourses/<Day>_run-<RunN>/{gocue,saclocked}/ch<NNN>.png`) — there are too many files (~6400) for the by-type convention to help there, and browsing one session's channels together is more useful for that figure type specifically.

| File | Purpose |
|---|---|
| `figSavePath.m` | Centralizes the `figures/<type>/<Day>_run-<RunN>.png` naming convention. |
| `findSaccadeRuns.m` | Scans `Data_collection/m032/*/run-*` for `stim_ct_classic_8targ.m`, returns a struct array of run metadata. |
| `matchRunRecording.m` | Maps a behavioral run to its OpenEphys `experimentN/recordingM` folder. |
| `matchTTLtoLogEvents.m` | Matches TTL digital-line numbers to behavioral event names by rising-edge count. |
| `extractSaccadeSession.m` | Core extraction: MUA + eye channels, trial table, eye calibration. One `.mat` output per session. |
| `runAllExtraction.m` | Batch driver over all 12 runs, tolerant of per-run failures. |
| `loadAllSessions.m` | Loads all extracted `.mat` files into one struct array; attaches `S.channelArea` (§4.12) to each. |
| `screenResponsiveChannels.m` | Per-channel, per-session responsiveness test in 4 task epochs vs. baseline. |
| `computeDirectionTuning.m` | 8-direction tuning curve, preferred direction, permutation significance, for responsive channels. |
| `analyzeEyeData.m` | Eye calibration application, saccade detection, fixation/pupil/main-sequence analysis. |
| `checkArtifactCorrelation.m` | Trial-residual pairwise channel correlation vs. distance, to distinguish shared-artifact from local-neural signal (§6). |
| `sessionAnatomyInfo.m` | Day → burr hole / channel-selection scheme / planned target structures / estimated insertion depth lookup (§2.2, §4.13). |
| `forceLightTheme.m` | Forces white background / black text on QC figures. |
| `classifyTuning.m` | Combines the vector-sum permutation test with a Kruskal-Wallis omnibus test into 3 tuning categories (§4.7): directional / complex / untuned. |
| `analyzePlanningExecution.m` | Splits go-cue-relative and saccade-onset-locked windows into planning vs execution epochs, separately screened and tuned (§4.8). |
| `summarizePlanningExecution.m` | Quantifies the planning/execution split across all 12 sessions (§6.6): recruitment fractions, preferred-direction consistency. |
| `checkSharedComponentKinematics.m` | Correlates the dominant shared trial-residual component (PC1) against saccade kinematics (§6.5), and common-mode removal for kinematic-linked sessions. |
| `decodeDirection.m` | 8-way direction decoding (individual channels + population-size sweep), all 3 epochs (§4.9). |
| `decodeTemporalGeneralization.m` | Cross-epoch (train/test) decoding matrix (§4.10). |
| `summarizeDecoding.m` | Pools decoding results across all 12 sessions. |
| `plotDirectionTimecourses.m` | Per-channel mean MUA timecourse by direction, go-cue-locked and saccade-onset-locked (§4.11). |
| `decodeAfterCommonModeRemoval.m` | Decoding before/after common-mode removal, for the 4 kinematic-linked sessions (§5.3). |
| `decodeSlidingGeneralization.m` | Continuous sliding-window temporal generalization matrix (§4.10, §5.3). |
| `generateAnatomyTemplate.m` | Builds the fillable per-channel/cluster anatomy assignment template (§4.12). |
| `getChannelArea.m` | Looks up assigned area/hierarchical labels and estimated depth-below-surface for a session's channels (§4.12, §4.13). |
| `extractAreaTokens.m` | Parses candidate area abbreviation tokens out of a free-text anatomy label, for liberal multi-label pooling (§6.7). |
| `poolTuningByArea.m` | Pools responsiveness/tuning across all 12 sessions by area, liberally (§6.7). |
| `runFullPipeline.m` | Runs all of the above end to end. |

## 4. Key methodological decisions (and why)

### 4.1 TTL bit codes are matched empirically, not from `par_default.m`

`par_default.m` defines `Par.StimB=1, Par.TargetB=2, Par.RewardB=3, Par.MicroB=6, Par.CorrectB=7`, with a comment that `SaccadeB=4` and `TrialB=5` are handled by a separate "DasControl" system. **None of this matches what's actually on the wire.** The NI-DAQ digital port carries exactly lines `{1, 4, 6, 8}` in every one of the 12 sessions, no exceptions, and cross-checking rising-edge counts against `Log.events` type counts shows this assignment is constant across the whole dataset:

- line 4 = `trial_start`
- line 6 = `targ_start` (the go-cue / saccade instruction)
- line 8 = `correct` / `reward_start` / `reward_stop` (these three are logged the same number of times whenever every trial that reaches the target epoch is rewarded, so they can't be told apart by count alone — harmless, since the pipeline only needs `targ_start` and `trial_start`)
- line 1 = fires more often than any single logged event type; likely a per-attempt fixation/stim marker with no 1:1 logged counterpart.

`matchTTLtoLogEvents.m` does this count-matching **per run** rather than hard-coding line numbers, as a robustness measure. On two sessions (20260317, 20260318), `targ_start`'s count ties exactly with the `correct`/`reward` cluster (all logged the same number of times that day, since every trial that reached the target was correct), making the match ambiguous by count alone. This is resolved with a secondary heuristic in `extractSaccadeSession.m`: among tied candidate lines, pick whichever has the smaller median latency after the (unambiguous) `trial_start` line — `targ_start` reliably follows `trial_start` by ~0.6s, while `correct`/`reward` follow much later, after the saccade.

### 4.2 Run-to-recording matching

`matchRunRecording.m` matches a run to a recording by **closest absolute time** between the run's folder timestamp and the recording's `sync_messages.txt` "Software Time," within a 90s tolerance — not "closest *preceding*," because operator workflow was not perfectly consistent: usually the OE recording starts a couple of seconds *before* the run script, but on at least one day the recording started *after* the run folder was created instead (20260303/run-003's matching recording, `experiment2/recording3`, starts 6.3s after its run folder timestamp).

Each behavioral run maps to its own dedicated recording — verified across every run on all 12 days, no recording is ever shared between two runs. `extractSaccadeSession.m` still time-windows the TTL stream to the run's own approximate span (its start time plus its logged duration from `Log.events`, with generous ±20/+30s margins) before count-matching, as a defensive precaution.

### 4.3 Eye channels: empirical identification, not documented anywhere

There is no config file, log field, or channel-name metadata anywhere in the rig, log `Par` struct, or `structure.oebin` that says which of the 8 NI-DAQ analog channels (AI0–AI7) is eye X, eye Y, or pupil. This was inferred from the raw signal shape:

- **AI0, AI1**: clean, correlated step-changes of ~1–4V, time-locked to the go-cue, returning to baseline ~300–450ms later — the signature of a saccade-and-return. Assumed **Eye X, Eye Y**.
- **AI2, AI3**: sit near the ADC negative rail (−10 to −5V), with the raw value pinned at exactly −10.000V for extended stretches (classic blink artifact where the pupil tracker loses the pupil). Assumed **pupil-related**; AI2 is used as "pupil" in the pipeline, AI3 is unused.
- **AI4–AI7**: near-zero, small (~0.05V) crosstalk-looking modulation. Assumed unused/floating.

This is a data-driven inference, not a documented fact (flagged in every extracted session as `S.eyeChannelAssumption`), but well validated indirectly: the eye-position calibration below achieves R² of 0.84–0.97 against the known 8 target locations, which would not happen if AI0/AI1 weren't genuinely eye position.

### 4.4 Eye calibration: derived per-session from task structure

`Par.SCx/SCy/OFFx/OFFy` in the logs are a live GUI display scale, not a true volts→degrees calibration (identical to `par_default.m`'s hardcoded defaults in every session, never updated by a calibration routine). Instead, calibration is derived directly from task structure: since the target's position is known exactly on every trial, and on **correct** trials the animal's gaze is by definition on the target, the mean raw voltage during a post-saccade landing window (200–350ms post-go-cue — chosen because this task's dwell on target is brief, with the return saccade beginning by ~350–400ms) is regressed against the known target (x,y) using a full 2-D affine fit:

```
[dva_x, dva_y] = [V_x, V_y, 1] · M
```

A full affine (rather than two independent 1-D linear fits) is necessary because the raw channels show real cross-talk — a pure-vertical target measurably displaces the "X" channel and vice versa, consistent with a rotated/skewed camera axis in the analog eye tracker.

**Blink handling.** Blinks corrupt the raw eye **position** channels (X/Y), not just the pupil channel: during a blink, X/Y pin near a fixed hardware rail of −5.00V (identical across all 12 sessions, regardless of each session's own raw baseline, which itself varies from about −2.3V to +0.8V) for roughly 65–180ms at a time. `extractSaccadeSession.m` detects this by proximity to the known rail voltage (`abs(raw−(−5.00))<0.5V`) and excludes these samples from the calibration fit; `analyzeEyeData.m` interpolates over them before the zero-phase velocity filter (which cannot tolerate NaNs) and then re-masks them (±10ms padding) before anything is plotted or used for saccade onset/kinematics detection.

Some **moderate** (3–7dva) landing scatter remains on plenty of correct trials — that's expected, not an error: `Stm.TargWinSz` = 7dva, i.e. the task's own online correct/error judgment accepts any landing within 7dva of the target center, and real saccades commonly under/overshoot by a few degrees before a small corrective adjustment. The trajectory plots (`*_eye_summary.png`) draw a light-grey dotted circle of this radius around each target so this is visually obvious.

### 4.5 MUA extraction

Reuses the filter chain from `XinyuScripts/getmua.m` (the lab's existing pipeline for this same probe/task family): 300Hz high-pass → 12-way ADC-multiplexing phase correction → 5000Hz low-pass → spatial "destripe" high-pass across channels → rectify + 200Hz low-pass envelope → decimate to ~1kHz. Channels are reordered by depth using `NP_PROBE.ELECTRODE_YPOS` from `settings.xml`. Trials are windowed 1.0s before to 1.2s after the go-cue (`targ_start`), covering the full fixation→stimulus→delay→saccade→return sequence.

### 4.6 Responsiveness screening

Per channel, per session: baseline-subtracted mean MUA in 4 windows (visual [-0.5, -0.35s], delay [-0.3, -0.05s], peri-saccadic [0, 0.25s], post-saccadic [0.3, 0.6s]) vs. a pre-stimulus baseline [−0.9, −0.6s], tested with a Wilcoxon signed-rank test across correct trials, Benjamini-Hochberg FDR-corrected across channels per window, plus a minimum effect-size requirement (|mean diff| > 0.3 baseline SD) to avoid flagging trivially small but "significant" effects given large trial counts. A channel is "responsive" if any window passes both criteria.

### 4.7 Direction tuning

For responsive channels only: mean baseline-subtracted response in the peri-saccadic window [0, 0.3s] relative to go-cue (see §4.8 for a version split into planning/execution epochs), grouped by the 8 target directions (only correct trials). Preferred direction and tuning strength via vector sum (resultant length, using rectified responses as weights); significance via a **permutation test** (1000 shuffles of direction labels, comparing observed resultant length to the shuffled null).

The vector-sum/permutation test is only sensitive to a **unimodal** direction effect — one clear preferred direction. A channel that responds to two roughly opposite directions (bimodal), or broadly across several adjacent directions without a single dominant peak, can have a real, statistically robust direction effect and still fail this test, because opposing vector contributions cancel in the sum. `classifyTuning.m` addresses this by also running a **Kruskal-Wallis omnibus test** across the 8 direction groups (does response differ across directions *at all*, regardless of shape) and combining the two into three categories, used for all tuning figures:

- **`directional`** (red): vector-sum permutation test significant — has a genuine single preferred direction; `prefDir`/`resultantLen` are meaningful summaries.
- **`complex`** (orange): Kruskal-Wallis significant but the vector-sum test is not — real direction-dependent structure exists, but not a single-peaked one (visibly bimodal or multi-peaked in the curve figures); `prefDir`/`resultantLen` are **not** meaningful for these, look at the actual per-direction response instead.
- **`untuned`** (gray): neither test significant.

### 4.8 Dissociating planning (premotor) from execution (perimovement) tuning

The [0, 0.3s] go-cue-relative window above mixes two different things whenever saccade latency is non-negligible (median 68–218ms across sessions, §5.1): **planning/premotor** activity (from the go-cue instruction until the eye actually starts moving) and **execution/perimovement** activity (during and immediately after the saccade itself). `analyzePlanningExecution.m` splits this two ways, producing figures separate from the main pipeline outputs:

- **`gocue` mode** — fixed split relative to the go-cue, same cutoff for every trial: planning=[0, 0.1]s, execution=[0.1, 0.3]s.
- **`saclocked` mode** (the methodologically preferred one, since latency varies so much across sessions) — anchored to each **trial's own detected saccade onset**: planning=[-0.15, -0.02]s and execution=[-0.02, +0.15]s relative to onset. Trials with no detected saccade (~1–11% depending on session) are excluded from this mode.

For each mode: responsiveness (matching §4.6) and direction tuning (§4.7's three-category scheme) are computed **separately** for the planning and execution windows. Outputs per session, both modes: `*_planexec_<mode>.png` (per-channel planning-vs-execution effect map, tuning-strength scatter, preferred-direction scatter) and `*_planexec_<mode>_curves.png` (every responsive channel's actual tuning curve, planning on top/execution below, category-colored). Results quantified across all 12 sessions are in §6.6.

### 4.9 Decoding: can the recorded population actually tell the 8 directions apart?

Stimulus position and saccade direction are the same label in this task — the target position on a given trial is exactly the instructed saccade direction — so this is one 8-way classification question, asked separately of the same three epochs as §4.8: decoding from the **visual** epoch asks whether stimulus position is represented; decoding from **planning**/**execution** asks the same of the upcoming/ongoing saccade. `decodeDirection.m` implements this as a direct complement to the tuning-based analyses above (which test each channel individually against a parametric null), giving a model-free readout of how much direction information is actually present, and how it scales with the number of channels pooled together.

**Decoder**: nearest-centroid classification on per-channel z-scored responses, evaluated with leave-one-trial-out cross-validation (train on all trials but one, classify the held-out trial by its Euclidean distance to each of the 8 per-direction training centroids, repeat for every trial). Chosen over LDA/SVM because population sizes here range up to ~380 channels while trial counts are only ~85–145 per session — far more features than trials, a regime where a full-covariance classifier is singular; per-channel standardization plus nearest-centroid is the standard robust choice here (equivalent to a diagonal/naive-Bayes LDA) and needs no hyperparameters. Chance level for 8 balanced classes is 1/8 = 12.5%.

Two complementary views, per session per epoch:
- **Individual channels**: every channel decoded alone (not just responsive ones, so unresponsive channels serve as an internal near-chance control).
- **Population**: accuracy vs. population size, built by randomly subsampling 1, 2, 4, ..., up to the full responsive population (20 repeats per size, averaged), the standard "how many channels before performance saturates" curve.

`summarizeDecoding.m` pools this across all 12 sessions (restricted to the common population sizes shared across sessions, since each session's own full-responsive-population size is idiosyncratic and not itself poolable).

### 4.10 Temporal generalization: is it the same code in each epoch, or a different one?

Direct follow-up to §4.9: decoding well in an epoch only shows that direction information is present *at that time* — it doesn't say whether the population structure that carries it (which channels prefer which direction) is the *same* structure across epochs, or gets reconfigured. `decodeTemporalGeneralization.m` tests this directly: train the nearest-centroid decoder on one epoch's population response, test it on **another** epoch's response for the same held-out trials (leave-one-trial-out, full responsive population), for every (train, test) pair among visual/planning/execution — a 3×3 matrix per session. The diagonal (train and test the same epoch) is exactly §4.9's ordinary same-epoch accuracy; the off-diagonal cells are the new information: if a decoder trained on one epoch still classifies well on another, the code generalizes (shared structure); if it drops to chance, the code is epoch-specific (reconfigured). Standardization (z-scoring) is always computed from the **train** epoch's own statistics and applied to the test epoch's data, so this is a genuine test of whether the trained decision boundary transfers, not just whether both epochs are independently informative.

This first pass uses the same 3 fixed epoch windows as §4.8/§4.9 (a coarse 3×3 matrix). `decodeSlidingGeneralization.m` extends this to a **continuous** version: 50ms bins spanning [-0.5, +0.8]s relative to go-cue (26×26 bins), same method, full responsive population, one matrix per session. Results in §5.3.

### 4.11 Per-channel direction timecourses

A direct, model-free complement to the tuning/decoding analyses above: for every responsive channel, the mean MUA timecourse (correct trials only) split into the 8 directions and overlaid, rather than reduced to a single-window number. Two alignments, `plotDirectionTimecourses.m`:

- **go-cue-locked**: the trial's own time base, window [-0.9, +1.0]s relative to go-cue — wide enough to show fixation, stimulus flash, delay, go-cue, target-end, and the return saccade all in one plot.
- **saccade-onset-locked**: each trial's MUA linearly interpolated onto a common time axis relative to its own detected saccade onset (necessary because onset time varies trial to trial — a simple index shift would misalign different trials' actual event times), window [-0.5, +0.8]s relative to onset.

Each direction's mean trace is smoothed (zero-phase Gaussian, 50ms window) before plotting — purely a display transform on the already-trial-averaged trace, not a change to `S.MUA` or any other analysis in this pipeline. A lighter 20ms window was tried first and left the trace visually indistinguishable from unsmoothed, since the underlying MUA is itself a 200Hz-bandwidth envelope with real high-frequency content; 50ms visibly reduces trial-count noise while still preserving genuine trial-to-trial structure (checked against an unsmoothed reference before choosing this value).

This produces one figure per responsive channel per session per alignment (~3200 channels × 2 ≈ 6400 figures total), organized into subfolders rather than the flat `figures/` directory used elsewhere in this pipeline:

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

**Filled in and wired into the pipeline.** Every cluster in `anatomy_template_clusters.csv` now has a hand-assigned `AssignedArea`, cross-referenced against `anatomy_template_clusters_hierarchical.csv` (a parallel hierarchical parcellation, L6=finest to L3=coarsest, per cluster). Both are propagated down to every individual channel in `anatomy_template_channels.csv` (`ClusterID`, `AssignedArea`, `L6`-`L3` columns, all 4608 channel rows across the 12 sessions, verified zero gaps/overlaps in cluster coverage).

`getChannelArea.m` (new) looks up these labels for a given day/run/channel(s) directly, and `loadAllSessions.m` now attaches this automatically as `S.channelArea` for every session — `S.channelArea.AssignedArea{c}` (finest tentative label) or `S.channelArea.L6/.L5/.L4/.L3{c}` (hierarchical, coarser toward L3), for channel `c` (1 = deepest, matching every other per-channel array in this pipeline, e.g. `S.CHdepthUm`). These are **tentative, hand-assigned labels compared against an atlas by eye, propagated per-cluster** — not a verified histological registration — but now available to any analysis in this pipeline without re-deriving them.

### 4.13 Estimated absolute depth below the brain surface

§2.2 flagged that absolute per-channel depth-in-brain (mm below cortical/dura surface, as opposed to `S.CHdepthUm`'s relative depth along the probe) wasn't reconstructable from the logs. Partial fix: for the 2 sessions where the elab channel-selection string encodes an explicit final depth (20260324 "Ch2H2@53mm", 20260325 "Ch2H1@53" → 53mm), that value tracks the burr hole's Chamber2Holes **planned** depth closely — so for the other 10 sessions (no explicit number logged), the planned depth is used as the estimated insertion depth instead. `sessionAnatomyInfo.m` now returns `estimatedInsertionDepthMm` (53mm for 20260324/20260325, else that day's burr-hole planned depth) and `depthSource` ('explicit' or 'planned').

`getChannelArea.m` uses this to compute, per channel: `DepthBelowSurfaceMm = estimatedInsertionDepthMm - CHdepthUm/1000` (channel 1, the tip, gets the full insertion depth; more distal channels get progressively less). A **negative** value means that channel is estimated to sit *above* the brain surface that day — in the guide tube, dura, or CSF, not brain tissue. Checked across all 12 sessions: several show a chunk of shallow channels estimated above surface (e.g. 20260224: 49/384, 20260226: 31/384, 20260317: 7/384), others show none (e.g. 20260312, 20260323, 20260324, 20260325, 20260326) — plausible given each day's online channel-selection scheme spans a different fraction of the probe relative to how far it was actually driven in.

This is still an **estimate**, not a verified measurement — it inherits whatever error exists in the Chamber2Holes plan vs. actual insertion for the 10 sessions using the fallback, and the underlying burr-hole targeting itself. Available as `S.channelArea.DepthBelowSurfaceMm(c)`, `.EstimatedInsertionDepthMm`, `.DepthSource` for every channel via `loadAllSessions.m`.

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

Tuned counts can differ by a few from run to run because `computeDirectionTuning.m`'s permutation test (1000 random shuffles, no fixed seed) gives slightly different results for channels sitting right at the significance boundary — expected noise, not a bug.

Note the clear **latency drop** from ~170–195ms in the earlier sessions to ~70–120ms from 20260320 onward — plausibly a training/practice effect (faster, more automatic responses later in the recording series).

Saccade detection (velocity-threshold on a low-pass-filtered, calibrated eye trace, blink periods excluded) succeeded on 99% of trials (1202/1213); endpoint accuracy and main-sequence (peak velocity vs. amplitude) relationships were physiological in every session inspected (see `*_eye_summary.png`).

### 5.2 QC figures (per session, in `data/figures/`)

- **`*_overview.png`** — four panels, all sharing one channel-index y-axis with **channel 1 (deepest, nearest probe tip) at the bottom and the shallowest channel at the top**: (1) a real-depth ruler (electrode position in mm from the probe tip, so gaps/clusters from the day's online channel-selection scheme are visible directly), (2) the z-scored MUA heatmap with all 5 trial-epoch boundaries (fixation/stimulus/delay/go-cue/target-end) marked, computed per-session from that session's own `Stm` timing, (3) a per-channel, per-epoch significance map (visual/delay/peri/post, ms windows labeled on the axis) so you can see *which* window a channel's "responsive" flag refers to, (4) a text panel with that day's burr hole, planned target structures, and online channel-selection scheme (§2.2).
- **`*_eye_summary.png`** — 7 panels: trajectory spoke plots shown twice, once for all trials (error trials dimmed) and once for correct trials only with % correct in the title; both draw a light-grey dotted circle (radius = `Stm.TargWinSz`, the task's real hit-detection window) around each target, behind the traces; fixation scatter; main sequence; latency histogram (ms); pupil trace (ms); accuracy by direction.
- **`*_tuning_examples.png`** — all responsive channels, sorted by depth, one small tuning-curve tile per channel (min-max normalized per channel for display), border and line color reflecting the 3-way category from §4.7: red=`directional`, orange=`complex`, gray=`untuned`.
- **`*_tuning_matrix.png`** — all responsive channels, sorted by preferred direction, per-channel min-max normalized, plus a tuning-strength color strip and a 3-category color strip alongside.
- **`*_population_summary.png`** — generated per session (not pooled — each session is an independent penetration, often a different burr hole/target structure per §2.2, and the probe is not seated to the same depth across sessions). Preferred-direction rose plot and tuning-strength histogram, plus preferred direction vs. real electrode depth (mm from tip, this session only) — see §6.4.
- **`*_artifact_check.png`** — pairwise trial-residual correlation matrix across all channels (direction-tuning contribution removed), shown raw and sorted by real depth, plus mean correlation vs. inter-channel distance — see §6.4.
- **`*_shared_kinematics.png`** — per-channel PC1 loadings (the dominant shared trial-residual component) and its per-trial score plotted against saccade peak velocity and amplitude, plus (for the 4 sessions where this component is confirmed kinematics-linked) trial-residual correlation vs. distance before/after removing it — see §6.5.
- **`*_planexec_gocue.png` / `*_planexec_saclocked.png`** — per-channel effect map (visual/planning/execution, ms windows labeled), tuning-strength scatter (visual-vs-execution and planning-vs-execution), and preferred-direction scatter, for the go-cue-fixed and saccade-onset-locked window splits respectively (§4.8).
- **`*_planexec_gocue_curves.png` / `*_planexec_saclocked_curves.png`** — every responsive channel's actual tuning curve, stacked visual (top) / planning (middle) / execution (bottom), category-colored.
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
| **Mean ± SEM** | | | **27.6 ± 3.6%** | **24.4 ± 2.8%** | **31.8 ± 4.5%** |

**11/12 sessions decode direction clearly above chance in every epoch** — the two exceptions (20260224, 20260326) are the same two sessions already flagged as the weakest on nearly every other metric in this report (§5.1, §6.4, §6.5, §6.6). Execution consistently gives the highest accuracy, consistent with a genuine, growing perimovement direction signal (matching the tuning-recruitment finding in §6.6), though visual and planning are also clearly above chance in most sessions — stimulus position and premotor direction are both decodable well before the saccade itself.

Pooling the population-size sweep across sessions (`pooled_decode_population.png`, restricted to population sizes shared by multiple sessions — see §4.9) gives a clean, monotonic accuracy-vs-population-size curve for all three epochs, rising from chance at n=1 to ~29–40% by n=256 pooled channels, with execution above visual/planning at every population size from ~n=8 upward. Individual channels are informative but weak on their own (most single channels sit in the 10–25% range, well below the ceiling reached by pooling), so the direction information in this population is real but genuinely distributed across many channels rather than carried by a small number of highly informative ones.

<!-- POOLED_DECODE_FIGURE -->

**Temporal generalization matrix** (§4.10, `*_decode_genmatrix.png`), pooled mean across the same 10 sessions (excluding 20260224/20260326, chance = 12.5%):

| Train \ Test | visual | planning | execution |
|---|---|---|---|
| **visual** | 30.7% | 18.0% | 20.0% |
| **planning** | 19.3% | 26.4% | 20.9% |
| **execution** | 18.3% | 20.2% | 35.5% |

Every off-diagonal cell is clearly above the 12.5% chance level — there is real shared, cross-epoch structure, a decoder trained in one epoch is not useless in another — but every off-diagonal cell is also well below its row's own diagonal (e.g. training on execution and testing on execution gets 35.5%, but training on execution and testing on visual only gets 18.3%). So the population code is **neither fully shared nor fully epoch-specific**: part of the direction-discriminating structure persists across visual/planning/execution, and part is reconfigured at each epoch.

**Sliding (continuous) temporal generalization matrix** (§4.10, `*_decode_slidegen.png`, 50ms bins across [-0.5, +0.8]s relative to go-cue): every session shows the same qualitative pattern — a compact, elevated **diagonal band** concentrated roughly 100-450ms post-go-cue (peaking near 150-400ms, tracking each session's own saccade latency, §5.1), with narrow off-diagonal generalization (a decoder trained at one moment in this band transfers well only to nearby moments, not distant ones) and accuracy at or near chance everywhere before the go-cue and in the return/ITI period after ~500-600ms. This is a much more precise picture than the coarse 3-epoch matrix: the "reconfiguration" implied by that matrix's below-diagonal off-diagonal cells is not a single sharp switch between 3 discrete codes, but a population code that continuously drifts/refines through the peri-saccadic window, with the most stable (best-generalizing) structure concentrated in a few-hundred-ms band around the saccade itself.

**Decoding after common-mode removal** (§6.5's open item, closed): for the 4 kinematic-linked sessions, decoding accuracy at the full responsive population, before vs. after removing that session-and-epoch's own dominant shared component (`decodeAfterCommonModeRemoval.m`, same removal method as §6.5):

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

Most recording sites here are subcortical (§2.2). A per-channel trial-residual correlation check plus per-session population summaries provide evidence that at least some of the tuning is genuine, spatially-organized neural signal — but the picture is mixed, not a clean "all real" or "all artifact" result. Details below.

### 6.1 The high responsiveness fraction

74% of all channels across all sessions were flagged "responsive," reaching 98% in one session (20260323). Two things are relevant to interpreting this:

- **Online channel selection is not per-electrode, it's per bank/section**: `Kwibus_elab.pdf` logs an explicit operator-chosen scheme for every session — "quarter density," "4bankeven," "single column," "sections2/3" (§2.2 table). The operator was already restricting acquisition toward electrode banks that looked responsive online, so a high responsive fraction *within the selected set* is an expected consequence of that selection.
- **Most of these trajectories are subcortical at depth** (§2.2): thalamus (MD/CL/CM), pulvinar, SC, habenula, ZI — not cortex. A depth-dependent *laminar* latency gradient is a cortical-columnar signature, not something to expect from a probe that's mostly in thalamus/SC by design.

See §6.4 for a direct test of the shared-artifact possibility.

### 6.2 Preferred-direction clustering across sessions

All 12 sessions are the same hemisphere (a single chamber), so a consistent contralateral-field bias in preferred direction across sessions is what a real retinotopic/motor map (SC, pulvinar) would produce. This doesn't fully explain the tight ~90–150° clustering seen across sessions, but it means the sites are not anatomically unrelated in the way independent penetrations across hemispheres/regions would be.

### 6.3 Cross-session depth comparison

Sessions are independent penetrations, the probe is not seated to the same depth across sessions, and (§2.2) different sessions often target different burr holes/structures. The population summary is generated **per session** (`*_population_summary.png`), using real electrode depth in mm from the probe tip rather than an abstract cross-session channel index. All screening/tuning computation is strictly within-session (`tuning_results.mat` keyed by `Day`/`RunN`/`channel`).

### 6.4 Per-channel trial-residual correlation check (the direct artifact test)

`checkArtifactCorrelation.m` tests the shared-artifact hypothesis directly: for each channel pair, correlate the **trial-to-trial residual** peri-saccadic response (each direction-condition's own mean subtracted out first, so genuine shared direction tuning can't inflate this) against physical inter-channel distance. A uniform artifact riding on every channel predicts high correlation that does **not** decay with distance; local neural activity predicts correlation that decays with distance.

Mean pairwise correlation per session (`*_artifact_check.png`):

| Day | Mean corr | Day | Mean corr | Day | Mean corr |
|---|---|---|---|---|---|
| 20260224 | 0.078 | 20260312 | 0.213 | 20260323 | 0.235 |
| 20260226 | 0.112 | 20260317 | 0.078 | 20260324 | 0.100 |
| 20260303 | 0.203 | 20260318 | 0.179 | 20260325 | 0.057 |
| 20260310 | 0.099 | 20260320 | 0.077 | 20260326 | 0.021 |

Two findings, looking at `*_artifact_check.png` and `*_population_summary.png` together:

1. **Correlation decays sharply with distance in every session** (e.g. 20260323: ~0.55 at <200µm down to ~0.18 by 1mm+), the signature of local neural activity rather than a probe-wide shared signal. It does not decay all the way to zero — it plateaus at a small positive value (~0.1–0.2 in the higher-correlation sessions) even between far-apart channels, consistent with a smaller residual shared component (common reference, EMG, or motion) riding underneath genuine local structure.
2. **The per-session population summaries show genuine, session-specific spatial organization that a shared artifact would not produce**: 20260312 (burr hole 1, SC) shows a repeating sawtooth relationship between real depth and preferred direction across ~4 distinct depth clusters; 20260323 shows a smooth, continuous, near-monotonic depth-direction relationship; 20260320 (burr hole 4, thalamus) shows depth-banded clusters of consistent preferred direction. By contrast, **20260326** — the session with the lowest mean pairwise correlation (0.021) and lowest tuned fraction (8/156) — shows no depth structure at all, just scattered noise.

**Net assessment**: several sessions (20260312, 20260320, 20260323 especially) show good evidence of genuine, spatially-organized tuning. A smaller shared/common-mode component likely coexists across all sessions. 20260326 looks like it may have little real signal. Recommend treating tuning results session-by-session rather than as one uniform "832 tuned channels" figure.

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

1. **PC1 (the dominant shared fluctuation) explains a lot of the trial-to-trial residual variance in most sessions** — from 10% up to 74% (20260312) — and loading sign agreement across channels is high (79–99%) in every session except **20260326 (55%, essentially chance)**, meaning that session's residual variance isn't dominated by one common-mode axis at all, consistent with §6.4's finding that 20260326 has no organized structure.
2. **The kinematic correlation is real but not universal: exactly 4/12 sessions** (20260310, 20260320, 20260323, 20260324) show a significant, same-sign correlation with **both** peak velocity and amplitude — the signature of EMG/movement contamination. The other 8 sessions show no reliable relationship (20260312 shows a significant correlation with peak velocity only, not amplitude).

**Interpretation**: the shared component's origin is not uniform across the dataset. In at least 4 sessions there's direct, quantitative evidence it's kinematics-linked (EMG/motion-like); in the rest, either it has a different origin, or the correlation approach lacks power in sessions with less kinematic variance/trial count.

**Common-mode removal, for the 4 kinematic-linked sessions.** For 20260310, 20260320, 20260323, and 20260324, `checkSharedComponentKinematics.m` also subtracts the rank-1 PC1 reconstruction from the trial residual and recomputes §6.4's correlation-vs-distance relationship on the cleaned data (`*_shared_kinematics.png`, 4th panel, before/after overlay). In all 4 sessions, the long-distance correlation **plateau** (§6.4's evidence for a residual shared component) drops sharply after removal — e.g. 20260323 goes from a ~0.18–0.21 plateau at >1mm down to ~0.05–0.07 — while the short-distance correlation, though also reduced, remains clearly above zero. This confirms the plateau specifically was this kinematic component, and that at least some of the short-distance, depth-organized correlation in §6.4 is not attributable to it.

### 6.6 Planning-vs-execution tuning, across all 12 sessions

`summarizePlanningExecution.m` quantifies the saccade-onset-locked planning/execution split (§4.8) across all 12 sessions. Two sessions (20260224, 20260326 — already the weakest/noisiest sessions by every other metric here) have zero or one directional channels in either window and are excluded from the pooled stats below.

- **Pooled across the other 10 sessions: 835 channel-instances were `directional` in execution; only 47% of those were also `directional` during planning — 53% (439) only became directional during execution**, i.e. genuine recruitment of additional direction-selective activity around the movement, not just an amplitude increase on an already-tuned population.
- Conversely, of 606 `directional`-in-planning channel-instances, 65% remained `directional` into execution (35%, 210, were planning-only).
- **For channels tuned in both epochs, preferred direction is meaningfully consistent, not independent**: mean absolute angular difference (planning vs execution), pooled = **40.5°**, well below the ~90° expected for unrelated/random pairings, though far from 0°. Per-session values range from 14.5° (20260325, very consistent) to ~50–80° (20260226, 20260303, 20260317 — closer to random, though these sessions also have the fewest both-tuned channels).

### 6.7 Pooling responsiveness/tuning by area, across all 12 sessions

Direct follow-up to §6.1/§6.3, now that channels carry an actual (tentative) area label (§4.12) rather than just an opaque depth index: `poolTuningByArea.m` pools every channel across ALL 12 sessions **by area** instead of by session — the exact thing §6.3 said couldn't be done responsibly before ("sessions are independent penetrations... not seated to the same depth"). Pooling is deliberately **liberal**: many hand-assigned labels genuinely span multiple candidate areas (e.g. "Posterior Commisure (or Habenula or CM?)"), so a channel counts toward *every* area token found anywhere in its label (`AssignedArea`/`L6`/`L5`/`L4`/`L3`), not just its single primary guess — a channel can and does appear in more than one area's pool (`extractAreaTokens.m`, matching on the parenthetical abbreviation codes, e.g. "(SC)", "(Pul)", to collapse inconsistent capitalization of the same area across entries).

**SC (superior colliculus) stands out sharply from every other pooled area**: 90.5% responsive, **42.6% directional** — by a wide margin the highest directional-tuning fraction of any area with reasonable sample size (next highest: MRT/LMID/RTMID — likely-synonymous labels for the same adjacent midbrain reticular region — at 49.3%, actually higher, but on a smaller n=140; SG at 25.7%). Most other pooled areas sit at 5-20% directional despite often being similarly or more responsive (e.g. PCC 37.1% responsive but only 4.7% directional; PCGG 43.6% responsive, 5.6% directional) — high responsiveness does not, on its own, imply direction tuning, and the areas that DO show strong tuning are exactly the ones with an a priori motor/saccade role (SC, adjacent midbrain). This is independent evidence — from pooled tuning fraction, not just the a priori "these trajectories are mostly subcortical" argument in §6.1 — that responsiveness and direction-tuning are not uniformly distributed by chance/artifact, but concentrate in anatomically sensible places.

**Preferred-direction distribution, not a single "preferred direction per area."** There's no a priori reason a whole area should have one preferred direction: for a retinotopically organized structure like SC, preferred direction should track WHERE in the map a given penetration happened to land, not converge to a single value — and any pooled distribution here is itself shaped by which part of an area each of the (few) penetrations sampled, not necessarily an intrinsic property of the area. `pooled_tuning_by_area_rose.png` therefore shows the actual **distribution** of individual channels' preferred directions per area (8 bins, matching the task's 8 directions) rather than collapsing to one resultant vector — the resultant-length-weighted circular mean is kept as a dotted reference line for continuity with the summary table, but it is not the headline result. Some areas show a clearly unimodal, skewed distribution (SC/CO: strongly concentrated around ~90-135°, resultant concentration 0.86 over 247 channel-instances; MRT/LMID/RTMID similarly concentrated, 0.91, n=69 each); others show a visibly **bimodal** or multi-peaked shape that a single-vector summary would have hidden or misrepresented as "just noisy" (e.g. PCGG, AREA_23 — low resultant concentration 0.15-0.23 turns out to reflect two opposing lobes, not a uniform/random spread). Given only 12 sessions total and each area typically sampled by just 1-4 of them, none of these distributions should be read as a claim about the area's true population-level direction organization — only about what these particular penetrations happened to pick up.

Full per-area table and three figures saved to `data/figures/pooled/` and `data/tuning_pooled_by_area.mat`: `pooled_tuning_by_area_summary.png` (%responsive/directional/complex per area, with channel count `n=` in each label), `pooled_tuning_by_area_rose.png` (preferred-direction distribution per area, areas with ≥15 directional channels only), `pooled_tuning_by_area_legend.png` (abbreviation → full descriptive name → channel count, for the ~45 area tokens — too many to fit as readable labels on the other two figures). Given the liberal (multi-label) pooling, treat all of this as **evidence, not proof** — a channel's inclusion in a given area's pool reflects that area being one *plausible* candidate for that channel, not a confirmed identity.

## 7. Open items requiring your input

1. **Artifact investigation** (§6.1, §6.4, §6.5) — the shared-component/kinematics check gives a real but non-universal answer: 4/12 sessions show a clean EMG/movement signature, the rest don't. Not pursued further: the two weakest sessions overall (20260224, 20260326) are plausibly explained simply by penetration placement — these burr holes/trajectories were planned for a different, unrelated paradigm, not for saccade responsiveness specifically (§2.2), so a session landing in less saccade-sensitive territory is an expected outcome of the shared multi-paradigm penetration design, not necessarily a data-quality problem to chase.
2. **Absolute per-channel depth-in-brain** (mm below cortical surface) is not directly measured in the logs (§2.2) — only relative depth along the probe. **Done, two ways, both estimates rather than verified measurements**: (a) per-channel tentative area labels, hand-assigned against an atlas (§4.12), available via `S.channelArea.AssignedArea`/`.L6`-`.L3`; (b) an estimated absolute `DepthBelowSurfaceMm` per channel (§4.13), anchored from an explicit final depth where logged (20260324/20260325) or the Chamber2Holes planned depth otherwise (the other 10 sessions).
3. **Burr hole for 20260224 is genuinely ambiguous** (marked "5?" in the source elab log itself, §2.2) — confirmed not worth chasing further.
4. **Planning-vs-execution** (§4.8, §6.6) — if useful, the natural extension is checking whether the 53% execution-only recruitment fraction varies with burr hole/structure (§2.2) or with the §6.5 kinematic-link sessions.
5. **Responsiveness/tuning thresholds**: current thresholds (FDR 0.05, effect size 0.3 baseline SD for screening; permutation p<0.05 for tuning) are reasonable defaults; given §6.4's session-by-session variability, you may want to treat sessions individually rather than adjust one global threshold.
6. **Done** (§5.3): decoding after common-mode removal, for the 4 kinematic-linked sessions — accuracy *increases* in every epoch of every session, meaning the shared component was masking real signal, not itself carrying decodable direction information.
7. **Done** (§5.3): continuous sliding-window temporal generalization matrix — shows a compact diagonal band of strong, narrowly-generalizing decoding concentrated ~100-450ms post-go-cue, not a sharp switch between discrete codes.
8. **Done** (§6.7): pooled tuning by area across all 12 sessions — SC stands out sharply (42.6% directional vs. 5-20% for most other areas). Natural follow-up, not yet built: the same area-pooling idea applied to *decoding* rather than tuning is less straightforward (trial-level data can't be naively concatenated across sessions the way per-channel tuning summaries can), but could be approximated by restricting each session's own decoding population to just its area-X channels and comparing accuracy across sessions/areas — worth doing if the tuning-by-area result is interesting enough to chase further.

## 8. Recommended next steps

1. **Given §6.4/§6.5's mixed result**, treat sessions individually rather than pooling a single "tuned channel" count — prioritize 20260312/20260320/20260323/20260324 (clean depth-organized structure and/or confirmed kinematic-linked shared component) for any further biological claims, and treat 20260224/20260326 (weakest on nearly every metric, plausibly explained by penetration placement — see §7 item 1) with more skepticism.
2. **Done** (§6.5): common-mode removal, for the 4 kinematic-linked sessions — the long-distance correlation plateau drops sharply after removing the shared component, while short-distance correlation survives.
3. **Done** (§4.13): absolute per-channel depth-in-brain, estimated per session from an explicit or planned insertion depth anchor — see §7 item 2.
4. Once sessions are triaged by §6.4/§6.5, revisit direction tuning with a von Mises fit (not just vector sum) for a cleaner tuning-width estimate on the trustworthy sessions — the saccade-onset-locked alignment from §4.8 is already implemented, use it in place of go-cue alignment for this.
5. Check whether the 53% execution-only recruitment fraction (§6.6) varies systematically with burr hole/structure or with the §6.5 kinematic-link sessions.
6. **Done** (§4.9, §5.3): direction decoding, individual-channel and population, all 3 epochs, all 12 sessions.
7. **Done** (§5.3): decoding after common-mode removal, and the continuous sliding-window generalization matrix.
