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

    % Blink/dropout mask on the eye POSITION channels themselves (see
    % extractSaccadeSession.m/S.eyeBlinkMask for why this exists and is
    % separate from the pupil-only masking above): during a blink, X/Y
    % pin near a ~-5V rail well outside any real saccade's raw excursion,
    % producing huge spurious "saccades" to nonsense dva coordinates if
    % left in. filtfilt cannot tolerate NaNs (a single NaN corrupts the
    % ENTIRE zero-phase-filtered trial, not just nearby samples), so blink
    % samples are linearly interpolated before filtering, then re-masked
    % to NaN afterward everywhere they're used for onset search, plotting,
    % and kinematics -- i.e. filtering sees a continuous signal, but nothing
    % downstream ever treats an interpolated/blink sample as real gaze.
    if isfield(S, 'eyeBlinkMask')
        blinkMask = S.eyeBlinkMask;
    else
        blinkMask = false(nTrials, nSamp); % older extractions without this field: no masking possible
    end
    dvaXi = dvaX; dvaYi = dvaY;
    for t = 1:nTrials
        if any(blinkMask(t,:))
            dvaXi(t,:) = fillmissing(dvaX(t,:), 'linear', 'EndValues', 'nearest', 'SamplePoints', find(true(1,nSamp)), ...
                'MissingLocations', blinkMask(t,:));
            dvaYi(t,:) = fillmissing(dvaY(t,:), 'linear', 'EndValues', 'nearest', 'SamplePoints', find(true(1,nSamp)), ...
                'MissingLocations', blinkMask(t,:));
        end
    end

    % raw analog eye channels are noise-dominated sample-to-sample (30kHz
    % ADC noise gives spurious velocities of 1e5-1e6 dva/s if
    % differentiated directly) -- low-pass filter before computing
    % velocity, as is standard for saccade onset/offset detection.
    [b_eye, a_eye] = butter(2, 50/(Fs/2), 'low');
    dvaXs = filtfilt(b_eye, a_eye, dvaXi')';
    dvaYs = filtfilt(b_eye, a_eye, dvaYi')';

    % Re-mask blink samples (plus 10ms padding either side, to cover the
    % filter's transient response at the edge of a blink) to NaN in the
    % position traces actually used from here on -- nothing downstream
    % should see interpolated values as if they were real gaze.
    padSamp = round(0.010*Fs);
    blinkMaskPad = blinkMask;
    if padSamp > 0
        blinkMaskPad = conv2(double(blinkMask), ones(1,2*padSamp+1), 'same') > 0;
    end
    dvaXs(blinkMaskPad) = nan;
    dvaYs(blinkMaskPad) = nan;

    vel = hypot(diff(dvaXs,1,2), diff(dvaYs,1,2)) * Fs; % dva/s, nTrial x (nSamp-1); NaN where either endpoint is blink-masked
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
        fixSD(t) = mean(std([dvaXs(t,baseIdx); dvaYs(t,baseIdx)], 0, 2, 'omitnan'), 'omitnan');

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
% hsv(8) puts a pale yellow-green (RGB 0.5,1,0) at the 90 deg slot, which
% at the 0.2-0.35 trajectory-line alpha used below nearly disappears
% against the white background -- looked like "fewer trials" for that
% direction but was a color-visibility artifact, not a real count
% difference (verified directly against trialTable counts). Fixed
% palette instead: 8 hues at consistent, moderate saturation/luminance so
% every direction stays visible at low alpha on white.
cmap = [0.89 0.10 0.11; 1.00 0.50 0.00; 0.60 0.60 0.00; 0.20 0.60 0.20; ...
        0.00 0.60 0.60; 0.12 0.47 0.71; 0.58 0.00 0.83; 0.89 0.10 0.55];
pctCorrect = 100*sum(isCorrect)/numel(isCorrect);

f = figure('Visible','off','Position',[0 0 1900 950], 'Color', 'w');
tl = tiledlayout(f, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Dotted circle of radius Stm.TargWinSz around each target: the task's
% own online hit-detection window (any landing inside this circle counts
% as "correct"). Makes it visually obvious that a saccade landing several
% dva short of dead-center is still a legitimate hit, not a labeling
% error (see PIPELINE_REPORT sec 4.4).
targXY = cat(1, S.Stm.Targ.xy_dva);
targWinSz = S.Stm.TargWinSz;
circTheta = linspace(0, 2*pi, 60);
circX = targWinSz*cos(circTheta); circY = targWinSz*sin(circTheta);

% 1) trajectory spoke plot -- ALL trials (error trials included, dimmer)
% Hit-detection circles are drawn FIRST (light grey, behind everything)
% so they read as unobtrusive background context, not a foreground element.
ax1 = nexttile(tl); hold(ax1, 'on');
for i = 1:size(targXY,1)
    plot(ax1, targXY(i,1)+circX, targXY(i,2)+circY, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
end
for d = 1:8
    trl = find(dirDeg==dirsDeg(d));
    for k = trl'
        w = 0.35 * isCorrect(k) + 0.08 * ~isCorrect(k);
        plot(ax1, es.dvaX(k,:), es.dvaY(k,:), 'Color', [cmap(d,:) w]);
    end
end
scatter(ax1, targXY(:,1), targXY(:,2), 80, 'k', 'filled');
axis(ax1, 'equal'); xlim(ax1, [-16 16]); ylim(ax1, [-16 16]);
set(ax1, 'Color', 'w');
xlabel(ax1, 'X (dva)'); ylabel(ax1, 'Y (dva)');
title(ax1, sprintf('%s run-%03d: all trials (n=%d)', S.Day, S.RunN, numel(isCorrect)));

% 1b) trajectory spoke plot -- correct trials only, with %correct
ax1b = nexttile(tl); hold(ax1b, 'on');
for i = 1:size(targXY,1)
    plot(ax1b, targXY(i,1)+circX, targXY(i,2)+circY, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
end
for d = 1:8
    trl = find(isCorrect & dirDeg==dirsDeg(d));
    for k = trl'
        plot(ax1b, es.dvaX(k,:), es.dvaY(k,:), 'Color', [cmap(d,:) 0.35]);
    end
end
scatter(ax1b, targXY(:,1), targXY(:,2), 80, 'k', 'filled');
axis(ax1b, 'equal'); xlim(ax1b, [-16 16]); ylim(ax1b, [-16 16]);
set(ax1b, 'Color', 'w');
xlabel(ax1b, 'X (dva)'); ylabel(ax1b, 'Y (dva)');
title(ax1b, sprintf('Correct trials only (%.0f%% correct, n=%d); dotted circle = hit window (TargWinSz=%.0fdva)', pctCorrect, sum(isCorrect), targWinSz));

% 2) fixation scatter (baseline window)
ax2 = nexttile(tl);
baseIdx = S.tbEye >= -0.9 & S.tbEye < -0.6;
bx = mean(es.dvaX(:,baseIdx),2,'omitnan'); by = mean(es.dvaY(:,baseIdx),2,'omitnan');
scatter(ax2, bx, by, 10, 'filled'); axis(ax2, 'equal');
xlim(ax2, [-4 4]); ylim(ax2, [-4 4]); set(ax2, 'Color', 'w');
xlabel(ax2, 'X (dva)'); ylabel(ax2, 'Y (dva)');
title(ax2, sprintf('Fixation scatter (median SD=%.2f dva)', median(es.fixSD,'omitnan')));

% 3) main sequence
ax3 = nexttile(tl);
scatter(ax3, es.amplitude, es.peakVel, 15, 'filled');
set(ax3, 'Color', 'w');
xlabel(ax3, 'Amplitude (dva)'); ylabel(ax3, 'Peak velocity (dva/s)');
title(ax3, 'Saccade main sequence');

% 4) latency histogram (ms)
ax4 = nexttile(tl);
histogram(ax4, es.onsetTime*1000, 30);
set(ax4, 'Color', 'w');
xlabel(ax4, 'Saccade onset (ms from go-cue)'); ylabel(ax4, '# trials');
title(ax4, sprintf('Latency (median %.0fms)', median(es.onsetTime,'omitnan')*1000));

% 5) pupil trace aligned to go-cue (time axis in ms, per convention)
ax5 = nexttile(tl);
pm = mean(es.pupil(isCorrect,:), 1, 'omitnan');
ps = std(es.pupil(isCorrect,:), 0, 1, 'omitnan') ./ sqrt(sum(isCorrect));
tbEyeMs = S.tbEye * 1000;
plot(ax5, tbEyeMs, pm, 'k'); hold(ax5, 'on');
fill(ax5, [tbEyeMs fliplr(tbEyeMs)], [pm+ps fliplr(pm-ps)], 'k', 'FaceAlpha', 0.2, 'EdgeColor','none');
xline(ax5, 0, 'r--');
set(ax5, 'Color', 'w');
xlabel(ax5, 'Time from go-cue (ms)'); ylabel(ax5, 'Pupil channel (V, raw)');
title(ax5, 'Pupil (raw, blink-masked)');

% 6) endpoint error by direction
ax6 = nexttile(tl);
errByDir = nan(1,8);
for d = 1:8
    errByDir(d) = mean(es.endpointErr(isCorrect & dirDeg==dirsDeg(d)), 'omitnan');
end
bar(ax6, dirsDeg, errByDir);
set(ax6, 'Color', 'w');
xlabel(ax6, 'Direction (deg)'); ylabel(ax6, 'Endpoint error (dva)');
title(ax6, 'Saccade accuracy by direction');

forceLightTheme(f);
saveas(f, fullfile(opts.figDir, sprintf('%s_run-%03d_eye_summary.png', S.Day, S.RunN)));
close(f);
end
