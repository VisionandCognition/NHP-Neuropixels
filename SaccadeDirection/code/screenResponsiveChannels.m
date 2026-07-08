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

    plotSessionOverview(S, screening(s), opts, winNames, baseWin);
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
function plotSessionOverview(S, sc, opts, winNames, baseWin)
% sc: this session's row of the `screening` struct array (p, fdrP,
% effect, responsive, windowNames), so the epoch-specificity panel can
% show *which* window each channel was significant in, not just a binary
% responsive flag.
responsive = sc.responsive;
isCorrect = strcmp(S.trialTable.outcome, 'correct');
nChan = size(S.MUA, 1);
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

% Channel index 1 = deepest (nearest probe tip, smallest CHdepthUm);
% channel index nChan = shallowest (nearest brain surface / probe base).
% Plot with deep at the BOTTOM (YDir normal + smallest-index-first data),
% matching anatomical convention (surface up, tip down).
chIdx = (1:nChan)';

% Trial-epoch event times relative to targ_start (go-cue, t=0), from this
% session's own Stm timing (fixed across trials by task design).
tFix   = -(S.Stm.FIXT + S.Stm.STIMT + S.Stm.DELAYT)/1000;
tStim  = -(S.Stm.STIMT + S.Stm.DELAYT)/1000;
tDelay = -(S.Stm.DELAYT)/1000;
tGo    = 0;
tEnd   = S.Stm.TARGT/1000;
evT = [tFix, tStim, tDelay, tGo, tEnd];
evLabel = {'fixation', 'stimulus', 'delay', 'go-cue', 'target end'};
% NB: pick colors that stay visible against BOTH the white figure
% background/legend AND the parula heatmap underneath (parula spans dark
% blue -> teal -> green -> yellow) -- white/pale colors disappear against
% white, and blues/teals/greens disappear against parula's own midtones
% (this is why the original cyan "delay" line was hard to see). Stick to
% colors outside parula's gamut: orange, magenta, purple, black, dark red.
evColor = {[0.9 0.45 0], [1 0 1], [0.55 0 0.85], [0 0 0], [0.65 0.1 0.1]};

anat = sessionAnatomyInfo(S.Day);

f = figure('Visible', 'off', 'Position', [0 0 1700 700], 'Color', 'w');
tl = tiledlayout(f, 1, 12, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel 1 (narrow): real electrode depth ruler, aligned to the same
% channel-index y-axis as the heatmap, so gaps/clusters from the online
% channel-selection scheme (bank/section/quarter-density picks) are
% visible directly next to the response data.
axDepth = nexttile(tl, 1, [1 2]);
if isfield(S, 'CHdepthUm')
    plot(axDepth, S.CHdepthUm/1000, chIdx, 'k.-', 'MarkerSize', 8);
    xlabel(axDepth, 'Electrode depth (mm from tip)');
else
    text(axDepth, 0.5, 0.5, 'CHdepthUm not available', 'HorizontalAlignment', 'center');
end
ylim(axDepth, [1 nChan]);
set(axDepth, 'YDir', 'normal', 'Color', 'w');
ylabel(axDepth, 'Channel # (1 = deepest / probe tip)');
title(axDepth, 'Real depth');
grid(axDepth, 'on');

% --- Panel 2 (wide): full-probe z-scored MUA heatmap
axMain = nexttile(tl, 3, [1 5]);
imagesc(axMain, S.tb, chIdx, meanMUAz);
set(axMain, 'CLim', [-3 3], 'YDir', 'normal', 'Color', 'w');
colormap(axMain, parula);
cb = colorbar(axMain); cb.Label.String = 'z (baseline SD)';
xlabel(axMain, 'Time from go-cue (s)');
ylabel(axMain, 'Channel # (1 = deepest / probe tip)');
title(axMain, sprintf('%s run-%03d: mean MUA, z-scored to baseline (correct trials, deep=bottom)', S.Day, S.RunN));
ylim(axMain, [1 nChan]);
hold(axMain, 'on');
evHandles = gobjects(1, numel(evT));
for e = 1:numel(evT)
    xline(axMain, evT(e), '-', 'Color', evColor{e}, 'LineWidth', 2.5);
    evHandles(e) = plot(axMain, nan, nan, '-', 'Color', evColor{e}, 'LineWidth', 2.5); % legend proxy
end
legend(axMain, evHandles, evLabel, 'Location', 'northoutside', 'Orientation', 'horizontal', ...
    'Box', 'off', 'TextColor', 'k', 'FontSize', 8);

% --- Panel 3 (medium): which task epoch each channel is significant in
% (sign + magnitude of the effect, gray/transparent where not
% significant), using the same channel-index y-axis as panel 2.
axEpoch = nexttile(tl, 8, [1 3]);
sigMask = sc.fdrP < opts.alpha & abs(sc.effect) > opts.minEffect;
epochEffect = sc.effect;
epochEffect(isnan(epochEffect)) = 0;
imagesc(axEpoch, 1:numel(winNames), chIdx, epochEffect, 'AlphaData', double(sigMask));
set(axEpoch, 'CLim', [-2 2], 'YDir', 'normal', 'Color', [0.85 0.85 0.85]);
colormap(axEpoch, redblue(256));
set(axEpoch, 'XTick', 1:numel(winNames), 'XTickLabel', winNames, 'XTickLabelRotation', 30);
cb2 = colorbar(axEpoch); cb2.Label.String = 'effect (baseline SD)';
ylim(axEpoch, [1 nChan]);
title(axEpoch, 'Significant epoch(s) per channel');
ylabel(axEpoch, '');

% --- Panel 4 (narrow): text annotation, planned trajectory / channel
% selection for this day (qualitative -- see sessionAnatomyInfo.m header
% for why this is not a precise per-channel depth registration).
axText = nexttile(tl, 11, [1 2]);
axis(axText, 'off');
wrapWidth = 32; % chars -- this panel is narrow, long structure names must be wrapped or they clip past the tile edge invisibly
structLines = cellfun(@(s) strjoin(wrapText(s, wrapWidth), '\n'), anat.structures, 'UniformOutput', false);
structStr = strjoin(structLines, ' ->\n');
txt = sprintf(['Burr hole: %s\nPlanned depth: %.0f mm\nChannel selection:\n%s\n\n' ...
    'Planned trajectory\n(superficial -> deep):\n%s\n\n' ...
    '%d/%d channels responsive\n(qualitative only -- no verified\nper-channel depth registration,\nsee sessionAnatomyInfo.m)'], ...
    anat.burrHole, anat.plannedDepthMm, strjoin(wrapText(anat.channelSelection, wrapWidth), '\n'), ...
    structStr, sum(responsive), nChan);
text(axText, 0, 1, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Interpreter', 'none', 'Color', 'k');

forceLightTheme(f);
saveas(f, fullfile(opts.figDir, sprintf('%s_run-%03d_overview.png', S.Day, S.RunN)));
close(f);
end

% ------------------------------------------------------------------
function cmap = redblue(n)
% Simple blue-white-red diverging colormap (no extra toolbox dependency).
half = floor(n/2);
top = [linspace(0.1,1,half)', linspace(0.1,1,half)', ones(half,1)];
bot = [ones(n-half,1), linspace(1,0.1,n-half)', linspace(1,0.1,n-half)'];
cmap = [top; bot];
end

% ------------------------------------------------------------------
function lines = wrapText(str, width)
% Greedy word-wrap to a fixed character width (no font-metrics
% dependency -- good enough for a fixed-size QC text panel).
words = strsplit(str, ' ');
lines = {};
cur = '';
for i = 1:numel(words)
    if isempty(cur)
        trial = words{i};
    else
        trial = [cur ' ' words{i}];
    end
    if numel(trial) > width && ~isempty(cur)
        lines{end+1} = cur; %#ok<AGROW>
        cur = words{i};
    else
        cur = trial;
    end
end
if ~isempty(cur), lines{end+1} = cur; end
end
