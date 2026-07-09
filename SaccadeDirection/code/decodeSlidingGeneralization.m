function slideResults = decodeSlidingGeneralization(sessions, screening, opts)
% DECODESLIDINGGENERALIZATION Fine, continuous version of
% decodeTemporalGeneralization.m's coarse 3-epoch matrix: bins the whole
% trial into fixed-width windows (opts.binWidth, default 50ms) spanning
% opts.window (default [-0.5, 0.8]s relative to go-cue), and computes the
% same train/test cross-decoding accuracy for every (train bin, test bin)
% pair -- a continuous generalization matrix showing exactly WHEN the
% population code changes, rather than just confirming (via the 3 coarse
% epochs) that it does.
%
% Method identical to decodeTemporalGeneralization.m: nearest-centroid,
% leave-one-trial-out, full responsive population; standardization
% (z-scoring) always computed from the TRAIN bin's own statistics and
% applied to the TEST bin's data, so this is a genuine test of whether a
% decision boundary trained at one moment transfers to another moment,
% not just whether both moments are independently informative.
%
% slideResults = decodeSlidingGeneralization(sessions, screening, opts)
%   sessions, screening : as in decodeDirection.m
%   opts.window   : [t0 t1] s relative to go-cue (default [-0.5 0.8])
%   opts.binWidth : bin width, s (default 0.05)
%   opts.figDir   : default <saccadeDataDir()>/figures

if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'window', [-0.5 0.8]);
opts = setDefault(opts, 'binWidth', 0.05);
opts = setDefault(opts, 'baseWindow', [-0.9 -0.6]);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));
if ~isfolder(opts.figDir), mkdir(opts.figDir); end

edges = opts.window(1):opts.binWidth:opts.window(2);
nBins = numel(edges) - 1;
binCentersMs = (edges(1:end-1) + edges(2:end))/2 * 1000;

slideResults = struct([]);

for s = 1:numel(sessions)
    S = sessions(s);
    isCorrect = strcmp(S.trialTable.outcome, 'correct') & ~isnan(S.trialTable.direction_deg);
    dirGroup = round(S.trialTable.direction_deg(isCorrect)'/45) + 1;

    baseIdx = S.tb >= opts.baseWindow(1) & S.tb < opts.baseWindow(2);
    baseTrials = squeeze(mean(S.MUA(:, isCorrect, baseIdx), 3));

    respChans = find(screening(s).responsive);
    respBins = cell(nBins, 1);
    for b = 1:nBins
        idx = S.tb >= edges(b) & S.tb < edges(b+1);
        respBins{b} = squeeze(mean(S.MUA(respChans, isCorrect, idx), 3)) - baseTrials(respChans,:);
    end

    genMat = nan(nBins, nBins);
    for iTrain = 1:nBins
        for iTest = 1:nBins
            genMat(iTrain, iTest) = crossEpochLOOAcc(respBins{iTrain}, respBins{iTest}, dirGroup);
        end
    end

    slideResults(s).Day = S.Day; %#ok<AGROW>
    slideResults(s).RunN = S.RunN;
    slideResults(s).binCentersMs = binCentersMs;
    slideResults(s).genMat = genMat;
    slideResults(s).nResp = numel(respChans);

    fprintf('%s run-%03d: sliding generalization matrix done (%d x %d bins, %d responsive channels)\n', ...
        S.Day, S.RunN, nBins, nBins, numel(respChans));

    plotSlidingMatrix(S, binCentersMs, genMat, opts);
end

outFile = fullfile(fileparts(opts.figDir), 'decode_sliding_generalization.mat');
save(outFile, 'slideResults', '-v7.3');
fprintf('saved %s\n', outFile);
end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function acc = crossEpochLOOAcc(respTrain, respTest, dirGroup)
% Identical to decodeTemporalGeneralization.m's crossEpochLOOAcc.
nTrial = size(respTrain, 2);
N = nTrial;
sumAll = sum(respTrain, 2);
sumSqAll = sum(respTrain.^2, 2);
nClass = 8;
classSum = zeros(size(respTrain,1), nClass);
classCount = zeros(1, nClass);
for d = 1:nClass
    idx = dirGroup == d;
    classCount(d) = sum(idx);
    if classCount(d) > 0
        classSum(:,d) = sum(respTrain(:,idx), 2);
    end
end

correct = false(1, nTrial);
for t = 1:nTrial
    d0 = dirGroup(t);
    xTrain = respTrain(:,t);
    xTest = respTest(:,t);
    mu = (sumAll - xTrain) / (N-1);
    varLoo = (sumSqAll - xTrain.^2)/(N-1) - mu.^2;
    varLoo = max(varLoo, 0);
    sd = sqrt(varLoo); sd(sd==0) = eps;

    ztest = (xTest - mu) ./ sd;
    dists = inf(1, nClass);
    for d = 1:nClass
        cnt = classCount(d) - (d==d0);
        if cnt <= 0, continue; end
        csum = classSum(:,d) - xTrain*(d==d0);
        centroidZ = (csum/cnt - mu) ./ sd;
        dists(d) = sum((centroidZ - ztest).^2);
    end
    [~, predIdx] = min(dists);
    correct(t) = (predIdx == d0);
end
acc = mean(correct);
end

% ------------------------------------------------------------------
function plotSlidingMatrix(S, binCentersMs, genMat, opts) %#ok<INUSD>
f = figure('Visible','off','Position',[0 0 620 560], 'Color', 'w');
ax = axes(f);
imagesc(ax, binCentersMs, binCentersMs, genMat*100);
set(ax, 'YDir', 'normal', 'Color', 'w', 'CLim', [0 max(60, max(genMat(:)*100))]);
colormap(ax, parula);
cb = colorbar(ax); cb.Label.String = 'decode accuracy (%)';
hold(ax, 'on');
xline(ax, 0, 'w:', 'LineWidth', 1); yline(ax, 0, 'w:', 'LineWidth', 1);
plot(ax, [min(binCentersMs) max(binCentersMs)], [min(binCentersMs) max(binCentersMs)], 'w:', 'LineWidth', 0.75);
xlabel(ax, 'Test time from go-cue (ms)'); ylabel(ax, 'Train time from go-cue (ms)');
title(ax, sprintf('%s run-%03d: sliding temporal generalization matrix\n(nearest-centroid, LOO-CV, chance=12.5%%, full responsive population)', S.Day, S.RunN), 'FontSize', 9);
forceLightTheme(f);
saveas(f, fullfile(opts.figDir, sprintf('%s_run-%03d_decode_slidegen.png', S.Day, S.RunN)));
close(f);
end
