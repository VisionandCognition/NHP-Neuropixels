function S = extractSaccadeSession(runInfo, opts)
% EXTRACTSACCADESESSION Extract trial-aligned MUA + eye data for one
% 8-target saccade-direction run.
%
% S = extractSaccadeSession(runInfo, opts)
%   runInfo : one element of findSaccadeRuns() output
%   opts    : struct with optional fields
%       .preWin  (default 1.0 s before targ_start / go-cue)
%       .postWin (default 1.2 s after targ_start / go-cue)
%       .eyeChans   (default [0 1 2], zero-based AI index for [X Y Pupil])
%       .save    (default true) -- save result to saccadeDataDir()
%       .outDir  (default saccadeDataDir(), on the server -- see saccadeDataDir.m)
%
% Alignment strategy:
%   - TTL digital lines are matched to Log.events types by rising-edge
%     COUNT per run (see matchTTLtoLogEvents.m), not by the bit numbers
%     in par_default.m, which do not match what's actually wired here.
%   - Trials are anchored on 'targ_start' (the go-cue / target onset,
%     i.e. the saccade instruction), matched in temporal order to
%     Log.events entries of the same type to recover trial identity
%     (targ_idx -> direction) and outcome (correct / no 'correct' event).
%   - Eye channels (AI0/AI1 assumed EyeX/EyeY, AI2 assumed Pupil -- this
%     is an EMPIRICAL guess from raw trace shape, not documented anywhere
%     in the rig config; flagged in S.eyeChannelAssumption) are extracted
%     raw, then calibrated in degrees-of-visual-angle per session using
%     the known 8 target locations and each correct trial's post-saccade
%     landing voltage (see calibrateEyeSession subfunction).

if nargin < 2, opts = struct(); end
opts = setDefault(opts, 'preWin', 1.0);
opts = setDefault(opts, 'postWin', 1.2);
opts = setDefault(opts, 'eyeChans', [0 1 2]);
opts = setDefault(opts, 'save', true);
opts = setDefault(opts, 'outDir', saccadeDataDir());

fprintf('=== Extracting %s run-%03d ===\n', runInfo.Day, runInfo.RunN);

%% Load behavioral log
L = load(runInfo.LogFile, 'Log', 'Stm', 'Par');
Log = L.Log; Stm = L.Stm; Par = L.Par;

%% Locate matching OpenEphys recording
dayPath = fileparts(runInfo.RunPath);
rec = matchRunRecording(dayPath, runInfo.Day, runInfo.RunTimeStr);
fprintf('  matched to experiment%d/recording%d\n', rec.ExpN, rec.RecN);

apFs = rec.APInfo.sample_rate;
niFs = rec.NIDAQInfo.sample_rate;
apNChan = rec.APInfo.num_channels;
niNChan = rec.NIDAQInfo.num_channels;

%% Load TTL events (behavioral bit codes are on the NI-DAQ digital port)
ttl.lines = readNPY(fullfile(rec.NIDAQEventsTTLDir, 'states.npy'));
ttl.sampleNumbers = readNPY(fullfile(rec.NIDAQEventsTTLDir, 'sample_numbers.npy'));
ttl.timestamps = readNPY(fullfile(rec.NIDAQEventsTTLDir, 'timestamps.npy'));

%% Restrict TTL events to this run's own time span within the recording.
% A single OpenEphys recording can span MULTIPLE behavioral runs (observed
% directly: e.g. 20260326 experiment1/recording1 covers runs 1-3), so the
% raw TTL stream for the matched recording can contain events from other
% runs too. Localize using the run's start time (from its folder name,
% relative to the matched recording's absolute start) and its own Log
% duration, then apply the count-matching only within that window.
niFileTs = readNPY(fullfile(rec.NIDAQContinuousDir, 'timestamps.npy'));
niT0forWindow = niFileTs(1);
runDT = datetime([runInfo.Day runInfo.RunTimeStr], 'InputFormat', 'yyyyMMddHHmmss', ...
    'TimeZone', 'Europe/Amsterdam');
runStartT = niT0forWindow + seconds(runDT - rec.RecordingStartDT);
logEventT1 = arrayfun(@(e) e.t(1), Log.events);
runWindow = [runStartT - 20, runStartT + max(logEventT1) + 30];
inWindow = ttl.timestamps >= runWindow(1) & ttl.timestamps <= runWindow(2);
ttl.lines = ttl.lines(inWindow);
ttl.sampleNumbers = ttl.sampleNumbers(inWindow);
ttl.timestamps = ttl.timestamps(inWindow);
fprintf('  run window [%.1f %.1f]s within recording, %d TTL events in window\n', ...
    runWindow(1), runWindow(2), numel(ttl.lines));

%% Count Log.events types and match to TTL lines
eventTypes = {Log.events.type};
uTypes = unique(eventTypes);
eventTypeCounts = containers.Map('KeyType','char','ValueType','double');
for i = 1:numel(uTypes)
    eventTypeCounts(uTypes{i}) = sum(strcmp(eventTypes, uTypes{i}));
end
fprintf('  Log event types: %s\n', strjoin(cellfun(@(t) sprintf('%s=%d',t,eventTypeCounts(t)), uTypes, 'uni', 0), ', '));

lineMap = matchTTLtoLogEvents(ttl.lines, eventTypeCounts);

% targ_start's count sometimes ties with correct/reward_start/reward_stop
% (all logged the same number of times when every trial that reaches the
% target epoch is correct), which makes the pure count-match ambiguous.
% Break the tie using trial timing order instead: targ_start always
% follows trial_start soon after (within ~1s), whereas correct/reward
% follow much later (after the saccade completes) -- pick whichever tied
% candidate line has the smallest median latency after trial_start.
if ~isKey(lineMap, 'targ_start') && isKey(lineMap, 'trial_start')
    trialLineTmp = lineMap('trial_start');
    trialTimesTmp = ttl.timestamps(ttl.lines == trialLineTmp);
    nTarg = eventTypeCounts('targ_start');
    posLinesTmp = ttl.lines(ttl.lines > 0);
    uLinesTmp = unique(posLinesTmp);
    countsTmp = arrayfun(@(L) sum(posLinesTmp == L), uLinesTmp);
    candidateLines = uLinesTmp(countsTmp == nTarg);
    if numel(candidateLines) > 1
        medLatency = nan(size(candidateLines));
        for ci = 1:numel(candidateLines)
            cTimes = ttl.timestamps(ttl.lines == candidateLines(ci));
            lat = nan(size(cTimes));
            for k = 1:numel(cTimes)
                prior = trialTimesTmp(trialTimesTmp <= cTimes(k));
                if ~isempty(prior)
                    lat(k) = cTimes(k) - prior(end);
                end
            end
            medLatency(ci) = median(lat, 'omitnan');
        end
        [~, pick] = min(medLatency);
        lineMap('targ_start') = candidateLines(pick);
        fprintf('  TTL match: line %d <-> ''targ_start'' (tie-broken by timing, count=%d, median latency %.3fs vs %s)\n', ...
            candidateLines(pick), nTarg, medLatency(pick), mat2str(medLatency));
    end
end

if ~isKey(lineMap, 'targ_start') || ~isKey(lineMap, 'trial_start')
    error('extractSaccadeSession:noAlignmentBit', ...
        'Could not identify targ_start/trial_start TTL lines for %s run-%03d', ...
        runInfo.Day, runInfo.RunN);
end

targLine = lineMap('targ_start');
trialLine = lineMap('trial_start');
targTTLtimes = ttl.timestamps(ttl.lines == targLine);
trialTTLtimes = ttl.timestamps(ttl.lines == trialLine);

%% Match targ_start TTL events (in order) to Log.events 'targ_start' entries
logTargIdx = find(strcmp(eventTypes, 'targ_start'));
logCorrIdx = find(strcmp(eventTypes, 'correct'));
if numel(logTargIdx) ~= numel(targTTLtimes)
    error('extractSaccadeSession:countMismatch', ...
        'targ_start TTL count (%d) does not match Log.events count (%d)', ...
        numel(targTTLtimes), numel(logTargIdx));
end
targT1 = arrayfun(@(e) e.t(1), Log.events(logTargIdx));
[~, order] = sort(targT1);
logTargIdx = logTargIdx(order);
logTrialn = arrayfun(@(e) e.trialn(1), Log.events(logTargIdx));
correctTrialn = arrayfun(@(e) e.trialn(1), Log.events(logCorrIdx));

% map Trlnum -> trial struct (last write wins for duplicate/placeholder rows)
trialByNum = containers.Map('KeyType','double','ValueType','any');
for i = 1:numel(Log.trial)
    if isempty(Log.trial(i).Trlnum)
        continue % placeholder row, no trial number assigned yet
    end
    trialByNum(Log.trial(i).Trlnum) = Log.trial(i);
end

nTrials = numel(targTTLtimes);
targIdx = nan(nTrials,1);
distIdx = nan(nTrials,1);
outcome = repmat({'incomplete'}, nTrials, 1);
for i = 1:nTrials
    tn = logTrialn(i);
    if isKey(trialByNum, tn)
        trl = trialByNum(tn);
        targIdx(i) = trl.targ_idx;
        distIdx(i) = trl.dist_idx;
    end
    if any(correctTrialn == tn)
        outcome{i} = 'correct';
    else
        outcome{i} = 'error';
    end
end

% direction (deg) of each of the 8 target locations
targXY = cat(1, Stm.Targ.xy_dva);
targDirDeg = mod(atan2d(targXY(:,2), targXY(:,1)), 360);
directionDeg = nan(nTrials,1);
valid = ~isnan(targIdx);
directionDeg(valid) = targDirDeg(targIdx(valid));

fprintf('  %d trials (%d correct, %d error)\n', nTrials, ...
    sum(strcmp(outcome,'correct')), sum(strcmp(outcome,'error')));

%% Channel y-position order (probe geometry)
s = xml2struct(rec.SettingsFile);
sigChain = s.SETTINGS.SIGNALCHAIN;
if ~iscell(sigChain), sigChain = {sigChain}; end
CHorder = [];
for c = 1:numel(sigChain)
    procs = sigChain{c}.PROCESSOR;
    if ~iscell(procs), procs = {procs}; end
    for p = 1:numel(procs)
        if isfield(procs{p}, 'EDITOR') && isfield(procs{p}.EDITOR, 'CUSTOM_PARAMETERS') ...
                && isfield(procs{p}.EDITOR.CUSTOM_PARAMETERS, 'NP_PROBE')
            ypos = procs{p}.EDITOR.CUSTOM_PARAMETERS.NP_PROBE.ELECTRODE_YPOS.Attributes;
            ch = fieldnames(ypos);
            yp = zeros(numel(ch), 2);
            for k = 1:numel(ch)
                yp(k,2) = str2double(getfield(ypos, ch{k}));
                yp(k,1) = str2double(extractAfter(ch{k}, 'CH'));
            end
            yp = sortrows(yp, 2);
            CHorder = yp(:,1); % 0-based channel numbers, ordered by depth
            CHdepthUm = yp(:,2); % electrode Y position (um from probe tip), same order as CHorder
            break
        end
    end
    if ~isempty(CHorder), break; end
end
if isempty(CHorder)
    error('extractSaccadeSession:noProbeGeometry', 'Could not find NP_PROBE.ELECTRODE_YPOS in %s', rec.SettingsFile);
end

%% Memory-map continuous data
apFile = fullfile(rec.APContinuousDir, 'continuous.dat');
apBuf = memmapfile(apFile, 'Format', 'int16');
apSampleNumbers = readNPY(fullfile(rec.APContinuousDir, 'sample_numbers.npy'));
apTimestamps = readNPY(fullfile(rec.APContinuousDir, 'timestamps.npy'));
apT0 = apTimestamps(1);
apNSamp = numel(apTimestamps);

niFile = fullfile(rec.NIDAQContinuousDir, 'continuous.dat');
niBuf = memmapfile(niFile, 'Format', 'int16');
niTimestamps = readNPY(fullfile(rec.NIDAQContinuousDir, 'timestamps.npy'));
niT0 = niTimestamps(1);
niNSamp = numel(niTimestamps);
niBitVolts = 0.0003051851; % NI-DAQ AI, fixed per structure.oebin across sessions checked

%% Filters for MUA (as in XinyuScripts/getmua.m)
[b_hp, a_hp] = butter(3, 300/(apFs/2), 'high');
[b_lp, a_lp] = butter(3, 5000/(apFs/2), 'low');
[b_mualp, a_mualp] = butter(3, 200/(apFs/2), 'low');
[b_ds, a_ds] = butter(3, 0.01/1, 'high'); % spatial destripe filter across channels

% ADC multiplexing phase-shift correction (12-way interleave, as in getmua.m)
a = 1:24:384; b = 2:24:384;
cycle = nan(1,384);
for cshift = 1:12
    cycle(a+(cshift-1)*2) = cshift;
    cycle(b+(cshift-1)*2) = cshift;
end
ph = (rem(cycle-1,12)./13).*(1/apFs);

decimateFactor = round(apFs/1000); % ~1kHz MUA
preSampAP = round(opts.preWin*apFs);
postSampAP = round(opts.postWin*apFs);
sampsPerTrialAP = preSampAP + postSampAP + 1;
tb = (-preSampAP:postSampAP)/apFs;
tbds = tb(1:decimateFactor:end);
nDS = numel(tbds);

fmat = repmat((0:sampsPerTrialAP-1)'/(sampsPerTrialAP/apFs), 1, 384);
phmat = repmat(ph, sampsPerTrialAP, 1);

preSampNI = round(opts.preWin*niFs);
postSampNI = round(opts.postWin*niFs);
sampsPerTrialNI = preSampNI + postSampNI + 1;
tbEye = (-preSampNI:postSampNI)/niFs;

MUA = nan(384, nTrials, nDS, 'single');
eyeRaw = nan(numel(opts.eyeChans), nTrials, sampsPerTrialNI, 'single');
trialStartRel = nan(nTrials,1); % trial_start time relative to targ_start (s)

fprintf('  processing trials:');
for i = 1:nTrials
    if mod(i,20)==0, fprintf(' %d', i); end
    tEvent = targTTLtimes(i);

    % nearest preceding trial_start, for reference
    priorTrial = trialTTLtimes(trialTTLtimes <= tEvent);
    if ~isempty(priorTrial)
        trialStartRel(i) = priorTrial(end) - tEvent;
    end

    %% AP snippet
    apIdx0 = round((tEvent - apT0)*apFs) + 1 - preSampAP;
    apIdx1 = apIdx0 + sampsPerTrialAP - 1;
    if apIdx0 < 1 || apIdx1 > apNSamp
        continue % trial too close to recording edge
    end
    byteOffset = (apIdx0-1)*apNChan;
    nSampRead = sampsPerTrialAP*apNChan;
    raw = double(apBuf.Data(byteOffset+1 : byteOffset+nSampRead));
    cutdata = reshape(raw, apNChan, sampsPerTrialAP)';
    cutdata = cutdata(:,1:384)*0.195; % uV, drop sync channel

    cutdata = filtfilt(b_hp, a_hp, cutdata);
    fbuf = fft(cutdata);
    fbuf = exp(-1j*2*pi.*fmat.*phmat).*fbuf;
    cutdata = ifft(fbuf, 'symmetric');
    cutdata = filtfilt(b_lp, a_lp, cutdata);
    cutdata = filtfilt(b_ds, a_ds, cutdata')';
    env = filtfilt(b_mualp, a_mualp, abs(cutdata));
    MUA(:,i,:) = single(env(1:decimateFactor:end,:)');

    %% NI-DAQ eye snippet
    niIdx0 = round((tEvent - niT0)*niFs) + 1 - preSampNI;
    niIdx1 = niIdx0 + sampsPerTrialNI - 1;
    if niIdx0 < 1 || niIdx1 > niNSamp
        continue
    end
    byteOffsetNI = (niIdx0-1)*niNChan;
    nSampReadNI = sampsPerTrialNI*niNChan;
    rawNI = double(niBuf.Data(byteOffsetNI+1 : byteOffsetNI+nSampReadNI));
    niChunk = reshape(rawNI, niNChan, sampsPerTrialNI)'*niBitVolts;
    eyeRaw(:,i,:) = single(niChunk(:, opts.eyeChans+1)');
end
fprintf(' done\n');

MUA = MUA(CHorder+1, :, :); % reorder by depth

%% Reassemble trial table
trialTable = table((1:nTrials)', logTrialn(:), targIdx, distIdx, directionDeg, ...
    outcome(:), targTTLtimes(:), trialStartRel, ...
    'VariableNames', {'trialIdx','trialn','targ_idx','dist_idx','direction_deg', ...
    'outcome','targStartTime','trialStartRelTime'});

%% Blink/dropout mask for the eye channels (see calibrateEyeSession header
% for why this is needed): during a blink, the pupil channel is known to
% pin near the ADC rail (S.eyeChannelAssumption note above), but the
% raw X/Y position channels ALSO corrupt during the same event, pinning
% near a FIXED ~-5.00V hardware rail -- confirmed empirically identical
% (within noise) across all 12 sessions checked, regardless of each
% session's own raw baseline offset (which itself ranges from about -2.3V
% to +0.8V across sessions -- this rig has no fixed zero point).
%
% IMPORTANT: an earlier version of this mask used an ABSOLUTE threshold
% (|raw|>3V), reasoning that real saccades stay within ~1.5V of baseline.
% That reasoning about the excursion SIZE was fine, but applying it as an
% absolute cutoff was not: in sessions where the raw baseline itself sits
% close to +-3V (e.g. 20260224, baseline Y around -2.3V), a genuine
% target-directed saccade only needs to cross ~0.7V further to trip an
% absolute |raw|>3V test, causing real saccades (confirmed directly: a
% clean, physiological rise from baseline to the target, no blink-like
% plateau, pupil channel unchanged at the same moment) to be misidentified
% as blinks and masked out -- collapsing measured amplitude for whichever
% directions happened to push that session's baseline past the threshold
% (here, upward/upper-diagonal targets, since this session's Y baseline
% was already negative). Detecting proximity to the known, fixed rail
% value instead of a generic magnitude avoids this entirely, since the
% rail sits at the same voltage no matter what the session's baseline is.
railV = -5.00;
railTol = 0.5; % generous tolerance around the rail; real data in every session checked stays well outside [-5.5,-4.5]
eyeBlinkMask = squeeze(abs(eyeRaw(1,:,:) - railV) < railTol | abs(eyeRaw(2,:,:) - railV) < railTol); % nTrial x nSamp

%% Eye calibration (session-specific, data-driven -- see subfunction)
eyeCalib = calibrateEyeSession(eyeRaw, tbEye, trialTable, targXY, eyeBlinkMask);

%% Package output
S = struct();
S.Day = runInfo.Day;
S.RunN = runInfo.RunN;
S.RunPath = runInfo.RunPath;
% NB: runInfo.LogFile (below) is the verified behavioral log for this run
% (matched via findSaccadeRuns by RunPath). Do NOT use S.Par.LogFn /
% S.Par.LogPath for this -- those come from inside the log file's own
% saved Par struct and were found to be stale, reflecting a different
% (earlier, different-task) run in the same stim-program session, not
% the log this data was actually extracted from.
S.LogFile = runInfo.LogFile;
S.ExpN = rec.ExpN;
S.RecN = rec.RecN;
S.apFs = apFs;
S.niFs = niFs;
S.CHorder = CHorder;
S.CHdepthUm = CHdepthUm; % electrode Y position (um from probe tip), row-aligned with S.MUA / S.CHorder (row 1 = deepest/nearest tip)
S.MUA = MUA;
S.tb = tbds; % time (s) relative to targ_start (go-cue)
S.trialTable = trialTable;
S.eyeRaw = eyeRaw; % [chan x trial x time], chan order = opts.eyeChans
S.eyeChans = opts.eyeChans;
S.eyeChannelAssumption = 'AI0=EyeX, AI1=EyeY, AI2=Pupil -- inferred from raw trace shape, not documented on rig; verify before trusting absolute values';
S.tbEye = tbEye;
S.eyeCalib = eyeCalib;
S.eyeBlinkMask = eyeBlinkMask; % nTrial x nSamp logical, true = blink/dropout on the eye X/Y channels (see comment above); analyzeEyeData.m NaNs these samples before any downstream use
S.Stm = Stm;
S.Par = Par;
S.lineMap = lineMap;

if opts.save
    if ~isfolder(opts.outDir), mkdir(opts.outDir); end
    outFile = fullfile(opts.outDir, sprintf('%s_run-%03d_extracted.mat', S.Day, S.RunN));
    save(outFile, '-struct', 'S', '-v7.3');
    fprintf('  saved %s\n', outFile);
end

end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function calib = calibrateEyeSession(eyeRaw, tbEye, trialTable, targXY, eyeBlinkMask)
% Fit raw-voltage -> degrees-of-visual-angle calibration per session using
% correct trials' post-saccade landing voltage vs. the known target
% location for that trial (no independent calibration data exists in the
% logs, so this is inferred directly from task structure: on correct
% trials the animal's gaze is, by definition, on/near the instructed
% target during the landing window).
%
% Fits a full 2-D affine transform (not independent per-axis gains):
%   [tgtX tgtY] = [landX landY 1] * M
% because the raw AI0/AI1 channels show substantial cross-talk between
% axes (a pure-Y target displaces the "X" channel and vice versa) --
% consistent with a rotated/skewed analog eye-tracker camera axis.
%
% Blink/dropout samples (eyeBlinkMask, see caller) are excluded from the
% landing-window mean. Without this, a trial whose landing window happens
% to overlap a blink (common -- ~1/3 of trials in some sessions) averages
% in several ~-5V rail samples alongside genuine gaze samples, which can
% shift that trial's fitted "landing position" by many dva and, because
% these bad trials are part of the fit itself, distort the affine
% calibration for the WHOLE session, not just the contaminated trials.

landWin = tbEye > 0.20 & tbEye < 0.35; % post go-cue landing window
% NB: this task's saccade dwell is brief -- gaze reaches the target by
% ~150-200ms and often starts returning to center by ~350-400ms (checked
% empirically from raw traces), so the landing window must be narrow.
isCorrect = strcmp(trialTable.outcome, 'correct') & ~isnan(trialTable.targ_idx);
trls = find(isCorrect);

landX = nan(numel(trls),1);
landY = nan(numel(trls),1);
tgtX = nan(numel(trls),1);
tgtY = nan(numel(trls),1);
for k = 1:numel(trls)
    t = trls(k);
    blink = eyeBlinkMask(t, landWin);
    x = squeeze(eyeRaw(1,t,landWin)); x(blink) = nan;
    y = squeeze(eyeRaw(2,t,landWin)); y(blink) = nan;
    landX(k) = mean(x, 'omitnan');
    landY(k) = mean(y, 'omitnan');
    tgtX(k) = targXY(trialTable.targ_idx(t),1);
    tgtY(k) = targXY(trialTable.targ_idx(t),2);
end
ok = ~isnan(landX) & ~isnan(landY);

calib = struct();
if sum(ok) < 8
    calib.valid = false;
    calib.note = 'insufficient correct trials for calibration';
    return
end

A = [landX(ok), landY(ok), ones(sum(ok),1)];
Mx = A \ tgtX(ok); % tgtX = Mx(1)*landX + Mx(2)*landY + Mx(3)
My = A \ tgtY(ok);
predX = A*Mx;
predY = A*My;
r2X = 1 - sum((tgtX(ok)-predX).^2)/sum((tgtX(ok)-mean(tgtX(ok))).^2);
r2Y = 1 - sum((tgtY(ok)-predY).^2)/sum((tgtY(ok)-mean(tgtY(ok))).^2);

calib.valid = true;
calib.M = [Mx, My]; % 3x2: [landX landY 1] * M = [dva_x dva_y]
calib.r2X = r2X; calib.r2Y = r2Y;
calib.nCalibTrials = sum(ok);
calib.landWin = [0.20 0.35];
fprintf('  eye calib (2D affine): R2x=%.3f R2y=%.3f (n=%d correct trials)\n', r2X, r2Y, sum(ok));

end
