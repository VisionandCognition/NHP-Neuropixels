function screening = screenResponsiveChannels(sessions, opts)
% SCREENRESPONSIVECHANNELS Test every channel in every session for a
% response to any part of the saccade task (stim / delay / peri-saccadic
% / post-saccadic), relative to a pre-stimulus fixation baseline.
%
% screening = screenResponsiveChannels(sessions, opts)
%   sessions : struct array from loadAllSessions()
%   opts.windows : struct of [t0 t1] windows (s, relative to targ_start /
%                  go-cue). Defaults cover the fixed 8targ trial timing
%                  (FIXT 300 + STIMT 100 + DELAYT 400 ms before targ_start):
%       baseline [-0.9 -0.6], stim [-0.5 -0.35], delay [-0.3 -0.05],
%       peri [0 0.25], post [0.3 0.6]
%   opts.alpha       : FDR-corrected significance threshold (default 0.05)
%   opts.minEffect   : minimum |effect size| (window mean diff / baseline SD)
%                      required in addition to significance (default 0.3)
%   opts.figDir      : where to save QC figures
%
% Returns a struct array (one per session) with fields:
%   Day, RunN, CHorder, windowNames, p (nChan x nWin), fdrP, effect,
%   responsive (nChan x 1 logical, true if any window passes both criteria)

if nargin < 2, opts = struct(); end
opts = setDefault(opts, 'windows', struct( ...
    'baseline', [-0.9 -0.6], ...
    'stim',     [-0.5 -0.35], ...
    'delay',    [-0.3 -0.05], ...
    'peri',     [0    0.25], ...
    'post',     [0.3  0.6]));
opts = setDefault(opts, 'alpha', 0.05);
opts = setDefault(opts, 'minEffect', 0.3);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));

if ~isfolder(opts.figDir), mkdir(opts.figDir); end

winNames = setdiff(fieldnames(opts.windows), {'baseline'});
baseWin = opts.windows.baseline;

screening = struct([]);

for s = 1:numel(sessions)
    S = sessions(s);
    isCorrect = strcmp(S.trialTable.outcome, 'correct');
    nChan = size(S.MUA, 1);
    nWin = numel(winNames);

    baseIdx = S.tb >= baseWin(1) & S.tb < baseWin(2);
    baseMean = squeeze(mean(S.MUA(:, isCorrect, baseIdx), 3)); % nChan x nTrial
    baseSD = std(baseMean, 0, 2);

    p = nan(nChan, nWin);
    effect = nan(nChan, nWin);
    for w = 1:nWin
        win = opts.windows.(winNames{w});
        winIdx = S.tb >= win(1) & S.tb < win(2);
        winMean = squeeze(mean(S.MUA(:, isCorrect, winIdx), 3)); % nChan x nTrial
        diffs = winMean - baseMean;
        for c = 1:nChan
            d = diffs(c,:);
            d = d(~isnan(d));
            if numel(d) < 5 || all(d==d(1))
                continue
            end
            p(c,w) = signrank(d);
            effect(c,w) = mean(d) / max(baseSD(c), eps);
        end
    end

    % FDR correction per window, across channels
    fdrP = nan(size(p));
    for w = 1:nWin
        fdrP(:,w) = bhFDR(p(:,w));
    end

    responsive = any(fdrP < opts.alpha & abs(effect) > opts.minEffect, 2);

    screening(s).Day = S.Day;
    screening(s).RunN = S.RunN;
    screening(s).CHorder = S.CHorder;
    screening(s).windowNames = winNames;
    screening(s).p = p;
    screening(s).fdrP = fdrP;
    screening(s).effect = effect;
    screening(s).responsive = responsive;

    fprintf('%s run-%03d: %d/%d channels responsive\n', S.Day, S.RunN, sum(responsive), nChan);

    plotSessionOverview(S, responsive, opts, winNames, baseWin);
end

outFile = fullfile(fileparts(opts.figDir), 'screening_results.mat');
save(outFile, 'screening', 'opts');
fprintf('saved %s\n', outFile);

end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function q = bhFDR(p)
% Benjamini-Hochberg FDR correction, NaN-safe (NaNs pass through as NaN).
q = nan(size(p));
ok = ~isnan(p);
n = sum(ok);
if n == 0, return; end
[ps, order] = sort(p(ok));
ranks = (1:n)';
qs = ps .* n ./ ranks;
% enforce monotonicity from the largest p-value downward
qs = flipud(cummin(flipud(qs)));
qs = min(qs, 1);
qtmp = nan(n,1);
qtmp(order) = qs;
q(ok) = qtmp;
end

% ------------------------------------------------------------------
function plotSessionOverview(S, responsive, opts, winNames, baseWin)
isCorrect = strcmp(S.trialTable.outcome, 'correct');
meanMUA = squeeze(mean(S.MUA(:, isCorrect, :), 2)); % nChan x nTime
baseIdx = S.tb >= baseWin(1) & S.tb < baseWin(2);
base = mean(meanMUA(:, baseIdx), 2);
% per-trial baseline SD per channel, for z-scoring (raw MUA magnitude
% varies hugely across channels/depths, so a shared color scale on raw
% units hides all structure -- normalize by each channel's own baseline
% variability instead, as elsewhere in this pipeline).
baseTrials = squeeze(mean(S.MUA(:, isCorrect, baseIdx), 3));
baseSD = std(baseTrials, 0, 2);
meanMUAz = (meanMUA - base) ./ max(baseSD, eps);

f = figure('Visible', 'off', 'Position', [0 0 1400 600]);

subplot(1,3,1);
imagesc(S.tb, 1:size(meanMUAz,1), meanMUAz);
set(gca, 'CLim', [-3 3]);
colormap(gca, parula);
colorbar;
xlabel('Time from go-cue (s)'); ylabel('Channel (depth order)');
title(sprintf('%s run-%03d: mean MUA, z-scored to baseline (correct trials)', S.Day, S.RunN));
hold on;
xline(0, 'w--', 'go-cue');
for w = 1:numel(winNames)
    win = opts.windows.(winNames{w});
    xline(win(1), 'w:'); xline(win(2), 'w:');
end

subplot(1,3,2);
if any(responsive)
    imagesc(S.tb, 1:sum(responsive), meanMUAz(responsive,:));
    set(gca, 'CLim', [-3 3]);
    colorbar;
end
xlabel('Time from go-cue (s)'); ylabel('Responsive channel #');
title(sprintf('%d/%d responsive channels', sum(responsive), numel(responsive)));
xline(0, 'w--');

subplot(1,3,3);
barh(1:numel(responsive), double(responsive));
set(gca, 'YDir', 'reverse');
ylabel('Channel (depth order)'); xlabel('responsive');
title('Responsive channel map');
ylim([1 numel(responsive)]);

saveas(f, fullfile(opts.figDir, sprintf('%s_run-%03d_overview.png', S.Day, S.RunN)));
close(f);
end
