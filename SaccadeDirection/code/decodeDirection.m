function decodeResults = decodeDirection(sessions, screening, opts)
% DECODEDIRECTION Can recorded MUA actually be used to tell WHICH of the 8
% directions a trial was, in each task epoch? Stimulus position and
% saccade direction are the same label in this task (the target position
% IS the instructed saccade direction on every trial), so this is one
% 8-way classification problem, tested against the same three epochs used
% elsewhere in this pipeline (visual/planning/execution -- see
% analyzePlanningExecution.m): decoding from the visual epoch asks
% whether stimulus position is represented; decoding from planning/
% execution asks the same of the upcoming/ongoing saccade.
%
% Decoder: nearest-centroid on z-scored (per-channel, training-fold)
% responses, evaluated with leave-one-trial-out cross-validation. Chosen
% over LDA/SVM because population sizes here range up to ~300+ channels
% while trial counts are only ~85-145 -- far more features than trials,
% where a full-covariance classifier is singular/overfit. Nearest-centroid
% with per-channel (diagonal) standardization is the standard robust
% choice in this n<<p regime (equivalent to a diagonal/naive-Bayes LDA)
% and needs no hyperparameters. Chance level for 8 balanced classes is
% 1/8 = 12.5%.
%
% Two decoding questions, per session per epoch:
%   1. Individual channels: every channel decoded alone (not just
%      responsive ones, so unresponsive channels serve as an internal
%      near-chance control), to see which channels carry direction info
%      and how that lines up with the tuning categories in
%      computeDirectionTuning.m / classifyTuning.m.
%   2. Population: accuracy vs. population size, built by randomly
%      subsampling opts.popSizes channels FROM THE RESPONSIVE SET
%      (repeated opts.nReps times per size, averaged) plus one point at
%      the full responsive population -- the standard "how many channels
%      before performance saturates" curve.
%
% decodeResults = decodeDirection(sessions, screening, opts)
%   sessions  : struct array from loadAllSessions()
%   screening : struct array from screenResponsiveChannels() (for .responsive)
%   opts.baseWindow, opts.visualWindow, opts.goCueWindows : same defaults
%     as analyzePlanningExecution.m, so epochs match exactly
%   opts.popSizes : default [1 2 4 8 16 32 64 128 256]
%   opts.nReps     : random-subsample repeats per population size (default 20)
%   opts.figDir    : default <saccadeDataDir()>/figures

if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'baseWindow', [-0.9 -0.6]);
opts = setDefault(opts, 'visualWindow', [-0.5 -0.35]);
opts = setDefault(opts, 'goCueWindows', struct('planning', [0 0.1], 'execution', [0.1 0.3]));
opts = setDefault(opts, 'popSizes', [1 2 4 8 16 32 64 128 256]);
opts = setDefault(opts, 'nReps', 20);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));
if ~isfolder(opts.figDir), mkdir(opts.figDir); end

epochNames = {'visual', 'planning', 'execution'};
decodeResults = struct([]);

for s = 1:numel(sessions)
    S = sessions(s);
    isCorrect = strcmp(S.trialTable.outcome, 'correct') & ~isnan(S.trialTable.direction_deg);
    dirDeg = S.trialTable.direction_deg(isCorrect);
    dirGroup = round(dirDeg(:)'/45) + 1; % 1 x nTrial, values 1..8

    baseIdx = S.tb >= opts.baseWindow(1) & S.tb < opts.baseWindow(2);
    baseTrials = squeeze(mean(S.MUA(:, isCorrect, baseIdx), 3));

    visIdx = S.tb >= opts.visualWindow(1) & S.tb < opts.visualWindow(2);
    planIdx = S.tb >= opts.goCueWindows.planning(1) & S.tb < opts.goCueWindows.planning(2);
    execIdx = S.tb >= opts.goCueWindows.execution(1) & S.tb < opts.goCueWindows.execution(2);
    respVisual = squeeze(mean(S.MUA(:, isCorrect, visIdx), 3)) - baseTrials;
    respPlan   = squeeze(mean(S.MUA(:, isCorrect, planIdx), 3)) - baseTrials;
    respExec   = squeeze(mean(S.MUA(:, isCorrect, execIdx), 3)) - baseTrials;
    epochResp = {respVisual, respPlan, respExec};

    nChan = size(respVisual, 1);
    nTrial = numel(dirGroup);

    % --- 1. individual-channel decoding, every channel, every epoch
    chanAcc = nan(nChan, 3);
    for e = 1:3
        for c = 1:nChan
            chanAcc(c,e) = looDecodeAcc(epochResp{e}(c,:), dirGroup);
        end
    end

    % --- 2. population-size sweep, sampled from responsive channels
    respChans = find(screening(s).responsive);
    nResp = numel(respChans);
    popSizesUse = unique([opts.popSizes(opts.popSizes <= nResp), nResp]);
    popSizesUse = popSizesUse(popSizesUse >= 1);
    nP = numel(popSizesUse);
    popAcc = nan(nP, 3, opts.nReps);
    for e = 1:3
        for p = 1:nP
            P = popSizesUse(p);
            for r = 1:opts.nReps
                if P == nResp
                    chSel = respChans;
                else
                    chSel = respChans(randperm(nResp, P));
                end
                popAcc(p,e,r) = looDecodeAcc(epochResp{e}(chSel,:), dirGroup);
                if P == nResp, popAcc(p,e,:) = popAcc(p,e,r); break; end %#ok<AGROW> % deterministic -- no need to repeat
            end
        end
    end

    decodeResults(s).Day = S.Day;
    decodeResults(s).RunN = S.RunN;
    decodeResults(s).epochNames = epochNames;
    decodeResults(s).chanAcc = chanAcc;
    decodeResults(s).nChan = nChan;
    decodeResults(s).nTrial = nTrial;
    decodeResults(s).nResp = nResp;
    decodeResults(s).popSizes = popSizesUse;
    decodeResults(s).popAcc = popAcc;

    fprintf('%s run-%03d: full-population (n=%d responsive) decode accuracy -- visual %.1f%%, planning %.1f%%, execution %.1f%% (chance 12.5%%, n=%d trials)\n', ...
        S.Day, S.RunN, nResp, 100*mean(popAcc(end,1,:)), 100*mean(popAcc(end,2,:)), 100*mean(popAcc(end,3,:)), nTrial);

    plotDecodeChannels(S, chanAcc, opts);
    plotDecodePopulation(S, popSizesUse, popAcc, epochNames, opts);
end

outFile = fullfile(fileparts(opts.figDir), 'decode_results.mat');
save(outFile, 'decodeResults');
fprintf('saved %s\n', outFile);
end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function acc = looDecodeAcc(resp, dirGroup)
% Leave-one-trial-out nearest-centroid decoding accuracy.
% resp: nChan x nTrial baseline-subtracted response; dirGroup: 1 x nTrial
% labels 1..8. Uses closed-form leave-one-out updates to the per-channel
% mean/SD (standardization) and per-class centroids -- O(nChan) per fold
% instead of O(nChan*nTrial) from recomputing from scratch each time,
% since re-deriving all trial-loop decodes (up to ~380 channels x ~9
% population sizes x 20 reps x 3 epochs x 12 sessions) from scratch would
% otherwise be prohibitively slow.
nTrial = size(resp, 2);
N = nTrial;
sumAll = sum(resp, 2);
sumSqAll = sum(resp.^2, 2);
nClass = 8;
classSum = zeros(size(resp,1), nClass);
classCount = zeros(1, nClass);
for d = 1:nClass
    idx = dirGroup == d;
    classCount(d) = sum(idx);
    if classCount(d) > 0
        classSum(:,d) = sum(resp(:,idx), 2);
    end
end

correct = false(1, nTrial);
for t = 1:nTrial
    d0 = dirGroup(t);
    x = resp(:,t);
    mu = (sumAll - x) / (N-1);
    varLoo = (sumSqAll - x.^2)/(N-1) - mu.^2;
    varLoo = max(varLoo, 0);
    sd = sqrt(varLoo); sd(sd==0) = eps;

    ztest = (x - mu) ./ sd;
    dists = inf(1, nClass);
    for d = 1:nClass
        cnt = classCount(d) - (d==d0);
        if cnt <= 0, continue; end
        csum = classSum(:,d) - x*(d==d0);
        centroidZ = (csum/cnt - mu) ./ sd;
        dists(d) = sum((centroidZ - ztest).^2);
    end
    [~, predIdx] = min(dists);
    correct(t) = (predIdx == d0);
end
acc = mean(correct);
end

% ------------------------------------------------------------------
function plotDecodeChannels(S, chanAcc, opts)
% Manual axes positions (not tiledlayout) so the two scatter panels can
% each carry marginal histograms of the two axes' distributions -- a
% tiledlayout tile can't easily host a main+marginal sub-layout.
nChan = size(chanAcc, 1);
chIdx = (1:nChan)';
epochNames = {'visual', 'planning', 'execution'};
colors = {[0.2 0.4 0.8], [0.85 0.55 0.1], [0.75 0.1 0.1]};

f = figure('Visible','off','Position',[0 0 1750 620], 'Color', 'w');

ax1 = axes(f, 'Position', [0.045 0.10 0.24 0.80]); hold(ax1, 'on');
for e = 1:3
    scatter(ax1, chanAcc(:,e)*100, chIdx, 10, colors{e}, 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none');
end
xline(ax1, 12.5, 'k:');
xlabel(ax1, 'Single-channel decode accuracy (%)'); ylabel(ax1, 'Channel # (1 = deepest)');
legend(ax1, [epochNames, {'chance'}], 'Location', 'southeast', 'TextColor', 'k', 'Color', 'w');
title(ax1, 'Individual-channel decoding, by depth');
set(ax1, 'Color', 'w'); ylim(ax1, [1 nChan]);

plotScatterMarginals(f, 0.35, chanAcc(:,3)*100, chanAcc(:,2)*100, ...
    'Execution decode accuracy (%)', 'Planning decode accuracy (%)', 'Planning vs execution, per channel', [0.3 0.3 0.3]);
plotScatterMarginals(f, 0.68, chanAcc(:,1)*100, chanAcc(:,3)*100, ...
    'Visual decode accuracy (%)', 'Execution decode accuracy (%)', 'Visual vs execution, per channel', [0.3 0.3 0.3]);

st = sgtitle(f, sprintf('%s run-%03d: individual-channel direction decoding (nearest-centroid, LOO-CV, chance=12.5%%)', S.Day, S.RunN));
st.Color = 'k'; % forceLightTheme's generic text sweep doesn't reliably reach sgtitle's internal text object
forceLightTheme(f);
saveas(f, figSavePath(opts.figDir, 'decode_channels', S.Day, S.RunN));
close(f);
end

% ------------------------------------------------------------------
function plotScatterMarginals(f, x0, xdata, ydata, xlab, ylab, ttl, col)
% Scatter of (xdata,ydata) at horizontal offset x0 (normalized figure
% units), with a top marginal histogram of xdata and a right marginal
% histogram of ydata, both sharing the main axes' 0-100 range and chance
% line (12.5%) so the univariate distributions can be read directly
% alongside the joint scatter.
mainW = 0.17; mainH = 0.62; gap = 0.012; histD = 0.11;
y0 = 0.10;

axMain = axes(f, 'Position', [x0, y0, mainW, mainH]); hold(axMain, 'on');
scatter(axMain, xdata, ydata, 12, col, 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none');
xlim(axMain, [0 100]); ylim(axMain, [0 100]);
plot(axMain, [0 100], [0 100], 'k:');
xline(axMain, 12.5, 'k:'); yline(axMain, 12.5, 'k:');
xlabel(axMain, xlab); ylabel(axMain, ylab);
set(axMain, 'Color', 'w');

axTop = axes(f, 'Position', [x0, y0+mainH+gap, mainW, histD]);
histogram(axTop, xdata, 20, 'FaceColor', col, 'EdgeColor', 'none');
xlim(axTop, [0 100]);
xline(axTop, 12.5, 'k:');
set(axTop, 'Color', 'w', 'XTick', [], 'Box', 'off');
title(axTop, ttl, 'FontSize', 10, 'FontWeight', 'bold');

axRight = axes(f, 'Position', [x0+mainW+gap, y0, histD, mainH]);
histogram(axRight, ydata, 20, 'FaceColor', col, 'EdgeColor', 'none', 'Orientation', 'horizontal');
ylim(axRight, [0 100]);
yline(axRight, 12.5, 'k:');
set(axRight, 'Color', 'w', 'YTick', [], 'Box', 'off');
end

% ------------------------------------------------------------------
function plotDecodePopulation(S, popSizes, popAcc, epochNames, opts)
colors = {[0.2 0.4 0.8], [0.85 0.55 0.1], [0.75 0.1 0.1]};
f = figure('Visible','off','Position',[0 0 700 550], 'Color', 'w');
ax = axes(f); hold(ax, 'on');
hs = gobjects(1,3);
for e = 1:3
    m = squeeze(mean(popAcc(:,e,:), 3)) * 100;
    sdv = squeeze(std(popAcc(:,e,:), 0, 3)) * 100;
    fill(ax, [popSizes fliplr(popSizes)], [m'+sdv' fliplr(m'-sdv')], colors{e}, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    hs(e) = plot(ax, popSizes, m, '-o', 'Color', colors{e}, 'LineWidth', 1.8, 'MarkerFaceColor', colors{e}, 'MarkerSize', 4);
end
yline(ax, 12.5, 'k:', 'LineWidth', 1);
set(ax, 'XScale', 'log', 'Color', 'w');
xlabel(ax, 'Population size (# channels, log scale)'); ylabel(ax, 'Decode accuracy (%)');
ylim(ax, [0 100]);
legend(ax, hs, epochNames, 'Location', 'southeast', 'TextColor', 'k', 'Color', 'w');
title(ax, sprintf('%s run-%03d: population direction decoding vs. population size\n(nearest-centroid, LOO-CV, %d reps/size, chance=12.5%%)', S.Day, S.RunN, size(popAcc,3)));
forceLightTheme(f);
saveas(f, figSavePath(opts.figDir, 'decode_population', S.Day, S.RunN));
close(f);
end
