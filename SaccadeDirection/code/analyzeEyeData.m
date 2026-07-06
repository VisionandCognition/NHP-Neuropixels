function eyeSummary = analyzeEyeData(sessions, opts)
% ANALYZEEYEDATA Calibrate raw eye channels, detect saccades, and
% characterize saccade dynamics, fixation stability, and pupil size for
% every extracted session.
%
% eyeSummary = analyzeEyeData(sessions, opts)
%   sessions : struct array from loadAllSessions()
%   opts.velThresh  : saccade-detection velocity threshold, dva/s (default 60)
%   opts.minDur     : minimum saccade duration, s (default 0.008)
%   opts.minAmp     : minimum saccade amplitude, dva (default 2)
%   opts.searchWin  : window (s, rel. to targ_start) to search for the
%                     response saccade (default [-0.05 1.0])
%   opts.blinkThresh: |pupil raw V| above which samples are treated as a
%                     blink/dropout and set to NaN (default 9.5, close to
%                     the ADC negative rail seen in the raw traces)
%   opts.figDir     : figure output dir
%
% Returns a struct array (one per session) with per-trial saccade
% kinematics (onsetTime, duration, peakVel, amplitude, endpoint, error)
% plus fixation-stability and pupil summaries, and saves QC figures.

if nargin < 2, opts = struct(); end
opts = setDefault(opts, 'velThresh', 60);
opts = setDefault(opts, 'minDur', 0.008);
opts = setDefault(opts, 'minAmp', 2);
opts = setDefault(opts, 'searchWin', [-0.05 1.0]);
opts = setDefault(opts, 'blinkThresh', 9.5);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));

if ~isfolder(opts.figDir), mkdir(opts.figDir); end

dirsDeg = 0:45:315;
eyeSummary = struct([]);

for s = 1:numel(sessions)
    S = sessions(s);
    if ~S.eyeCalib.valid
        fprintf('%s run-%03d: no valid eye calibration, skipping\n', S.Day, S.RunN);
        continue
    end

    nTrials = size(S.eyeRaw, 2);
    nSamp = size(S.eyeRaw, 3);
    Fs = S.niFs;

    rawX = squeeze(S.eyeRaw(1,:,:));
    rawY = squeeze(S.eyeRaw(2,:,:));
    rawP = squeeze(S.eyeRaw(3,:,:));

    M = S.eyeCalib.M;
    dvaX = rawX*M(1,1) + rawY*M(2,1) + M(3,1);
    dvaY = rawX*M(1,2) + rawY*M(2,2) + M(3,2);

    pupil = rawP;
    pupil(abs(pupil) > opts.blinkThresh) = nan;

    % raw analog eye channels are noise-dominated sample-to-sample (30kHz
    % ADC noise gives spurious velocities of 1e5-1e6 dva/s if
    % differentiated directly) -- low-pass filter before computing
    % velocity, as is standard for saccade onset/offset detection.
    [b_eye, a_eye] = butter(2, 50/(Fs/2), 'low');
    dvaXs = filtfilt(b_eye, a_eye, dvaX')';
    dvaYs = filtfilt(b_eye, a_eye, dvaY')';

    vel = hypot(diff(dvaXs,1,2), diff(dvaYs,1,2)) * Fs; % dva/s, nTrial x (nSamp-1)
    tbVel = S.tbEye(1:end-1) + diff(S.tbEye)/2;

    searchIdx = tbVel >= opts.searchWin(1) & tbVel < opts.searchWin(2);
    baseIdx = S.tbEye >= -0.9 & S.tbEye < -0.6;

    onsetTime = nan(nTrials,1); offsetTime = nan(nTrials,1);
    peakVel = nan(nTrials,1); amplitude = nan(nTrials,1);
    endX = nan(nTrials,1); endY = nan(nTrials,1);
    fixSD = nan(nTrials,1);

    targXY = cat(1, S.Stm.Targ.xy_dva);
    endpointErr = nan(nTrials,1);

    for t = 1:nTrials
        v = vel(t,:);
        fixSD(t) = mean(std([dvaXs(t,baseIdx); dvaYs(t,baseIdx)], 0, 2), 'omitnan');

        vs = v; vs(~searchIdx) = -inf;
        aboveThresh = vs > opts.velThresh;
        onsetSamp = find(aboveThresh, 1, 'first');
        if isempty(onsetSamp), continue; end

        % extend onset backward to where velocity actually left baseline,
        % and find offset as first return below threshold after onset
        offsetSamp = onsetSamp + find(v(onsetSamp:end) < opts.velThresh, 1, 'first') - 1;
        if isempty(offsetSamp), offsetSamp = numel(v); end

        dur = tbVel(offsetSamp) - tbVel(onsetSamp);
        x0 = dvaXs(t, find(S.tbEye<=tbVel(onsetSamp),1,'last'));
        y0 = dvaYs(t, find(S.tbEye<=tbVel(onsetSamp),1,'last'));
        x1 = dvaXs(t, find(S.tbEye<=tbVel(offsetSamp),1,'last'));
        y1 = dvaYs(t, find(S.tbEye<=tbVel(offsetSamp),1,'last'));
        amp = hypot(x1-x0, y1-y0);

        if dur < opts.minDur || amp < opts.minAmp
            continue
        end

        onsetTime(t) = tbVel(onsetSamp);
        offsetTime(t) = tbVel(offsetSamp);
        peakVel(t) = max(v(onsetSamp:offsetSamp));
        amplitude(t) = amp;
        endX(t) = x1; endY(t) = y1;

        if ~isnan(S.trialTable.targ_idx(t))
            tgt = targXY(S.trialTable.targ_idx(t),:);
            endpointErr(t) = hypot(x1-tgt(1), y1-tgt(2));
        end
    end

    es = struct();
    es.Day = S.Day; es.RunN = S.RunN;
    es.onsetTime = onsetTime; es.offsetTime = offsetTime;
    es.duration = offsetTime - onsetTime;
    es.peakVel = peakVel; es.amplitude = amplitude;
    es.endX = endX; es.endY = endY; es.endpointErr = endpointErr;
    es.fixSD = fixSD;
    es.dvaX = dvaXs; es.dvaY = dvaYs; es.pupil = pupil;
    eyeSummary = [eyeSummary, es]; %#ok<AGROW>

    nDetected = sum(~isnan(onsetTime));
    fprintf('%s run-%03d: saccade detected in %d/%d trials, median latency %.3fs, median endpoint err %.2f dva\n', ...
        S.Day, S.RunN, nDetected, nTrials, median(onsetTime,'omitnan'), median(endpointErr,'omitnan'));

    plotEyeSummary(S, es, dirsDeg, opts);
end

outFile = fullfile(fileparts(opts.figDir), 'eye_summary.mat');
save(outFile, 'eyeSummary', 'opts');
fprintf('saved %s\n', outFile);

end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function plotEyeSummary(S, es, dirsDeg, opts)
isCorrect = strcmp(S.trialTable.outcome, 'correct');
dirDeg = S.trialTable.direction_deg;
cmap = hsv(8);

f = figure('Visible','off','Position',[0 0 1500 900]);

% 1) trajectory spoke plot
subplot(2,3,1); hold on;
for d = 1:8
    trl = find(isCorrect & dirDeg==dirsDeg(d));
    for k = trl'
        plot(es.dvaX(k,:), es.dvaY(k,:), 'Color', [cmap(d,:) 0.2]);
    end
end
targXY = cat(1, S.Stm.Targ.xy_dva);
scatter(targXY(:,1), targXY(:,2), 80, 'k', 'filled');
axis equal; xlim([-16 16]); ylim([-16 16]);
xlabel('X (dva)'); ylabel('Y (dva)');
title(sprintf('%s run-%03d: eye trajectories', S.Day, S.RunN));

% 2) fixation scatter (baseline window)
subplot(2,3,2);
baseIdx = S.tbEye >= -0.9 & S.tbEye < -0.6;
bx = mean(es.dvaX(:,baseIdx),2,'omitnan'); by = mean(es.dvaY(:,baseIdx),2,'omitnan');
scatter(bx, by, 10, 'filled'); axis equal;
xlim([-4 4]); ylim([-4 4]);
xlabel('X (dva)'); ylabel('Y (dva)');
title(sprintf('Fixation scatter (median SD=%.2f dva)', median(es.fixSD,'omitnan')));

% 3) main sequence
subplot(2,3,3);
scatter(es.amplitude, es.peakVel, 15, 'filled');
xlabel('Amplitude (dva)'); ylabel('Peak velocity (dva/s)');
title('Saccade main sequence');

% 4) latency histogram
subplot(2,3,4);
histogram(es.onsetTime, 30);
xlabel('Saccade onset (s from go-cue)'); ylabel('# trials');
title(sprintf('Latency (median %.3fs)', median(es.onsetTime,'omitnan')));

% 5) pupil trace aligned to go-cue
subplot(2,3,5);
pm = mean(es.pupil(isCorrect,:), 1, 'omitnan');
ps = std(es.pupil(isCorrect,:), 0, 1, 'omitnan') ./ sqrt(sum(isCorrect));
tbEye = S.tbEye;
plot(tbEye, pm, 'k'); hold on;
fill([tbEye fliplr(tbEye)], [pm+ps fliplr(pm-ps)], 'k', 'FaceAlpha', 0.2, 'EdgeColor','none');
xline(0,'r--');
xlabel('Time from go-cue (s)'); ylabel('Pupil channel (V, raw)');
title('Pupil (raw, blink-masked)');

% 6) endpoint error by direction
subplot(2,3,6);
errByDir = nan(1,8);
for d = 1:8
    errByDir(d) = mean(es.endpointErr(isCorrect & dirDeg==dirsDeg(d)), 'omitnan');
end
bar(dirsDeg, errByDir);
xlabel('Direction (deg)'); ylabel('Endpoint error (dva)');
title('Saccade accuracy by direction');

saveas(f, fullfile(opts.figDir, sprintf('%s_run-%03d_eye_summary.png', S.Day, S.RunN)));
close(f);
end
