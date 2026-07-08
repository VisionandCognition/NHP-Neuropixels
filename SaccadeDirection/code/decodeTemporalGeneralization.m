function genResults = decodeTemporalGeneralization(sessions, screening, opts)
% DECODETEMPORALGENERALIZATION Cross-epoch decoding: train the
% nearest-centroid decoder (decodeDirection.m) on one epoch's population
% response, test it on ANOTHER epoch's response for the same trials, for
% every (train, test) pair among {visual, planning, execution}. This
% tests whether the population code is a shared, temporally-generalizing
% code (a channel that prefers direction X during vision still prefers X
% during execution, so a visual-trained decoder still works at
% execution), or an epoch-specific code that gets reconfigured between
% epochs (population direction-discriminating structure changes, so
% cross-epoch accuracy drops toward chance even though within-epoch
% accuracy is high). The diagonal of the resulting 3x3 matrix is exactly
% decodeDirection.m's full-population, same-epoch accuracy.
%
% Method: leave-one-trial-out as in decodeDirection.m, but the per-channel
% standardization (mean/SD) and per-class centroids are always built from
% the TRAIN epoch's response (excluding the held-out trial); the held-out
% trial's TEST-epoch response is standardized using the TRAIN epoch's
% mean/SD (not its own) before computing distance to the train-epoch
% centroids -- this is what "apply the trained model to new data" means
% concretely, and is what makes the diagonal reduce exactly to
% decodeDirection.m's ordinary same-epoch decoding.
%
% Uses the full responsive population (no population-size sweep here --
% see decodeDirection.m for that); this is a first, coarse 3x3 version
% using the same fixed epoch windows as analyzePlanningExecution.m/
% decodeDirection.m. A finer, continuous sliding-window generalization
% matrix across the whole trial is a natural follow-up if this looks
% informative.
%
% genResults = decodeTemporalGeneralization(sessions, screening, opts)
%   sessions, screening : as in decodeDirection.m
%   opts.baseWindow, opts.visualWindow, opts.goCueWindows : same defaults
%     as decodeDirection.m, so epochs match exactly
%   opts.figDir : default <saccadeDataDir()>/figures

if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'baseWindow', [-0.9 -0.6]);
opts = setDefault(opts, 'visualWindow', [-0.5 -0.35]);
opts = setDefault(opts, 'goCueWindows', struct('planning', [0 0.1], 'execution', [0.1 0.3]));
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));
if ~isfolder(opts.figDir), mkdir(opts.figDir); end

epochNames = {'visual', 'planning', 'execution'};
genResults = struct([]);

for s = 1:numel(sessions)
    S = sessions(s);
    isCorrect = strcmp(S.trialTable.outcome, 'correct') & ~isnan(S.trialTable.direction_deg);
    dirDeg = S.trialTable.direction_deg(isCorrect);
    dirGroup = round(dirDeg(:)'/45) + 1;

    baseIdx = S.tb >= opts.baseWindow(1) & S.tb < opts.baseWindow(2);
    baseTrials = squeeze(mean(S.MUA(:, isCorrect, baseIdx), 3));

    visIdx = S.tb >= opts.visualWindow(1) & S.tb < opts.visualWindow(2);
    planIdx = S.tb >= opts.goCueWindows.planning(1) & S.tb < opts.goCueWindows.planning(2);
    execIdx = S.tb >= opts.goCueWindows.execution(1) & S.tb < opts.goCueWindows.execution(2);
    respVisual = squeeze(mean(S.MUA(:, isCorrect, visIdx), 3)) - baseTrials;
    respPlan   = squeeze(mean(S.MUA(:, isCorrect, planIdx), 3)) - baseTrials;
    respExec   = squeeze(mean(S.MUA(:, isCorrect, execIdx), 3)) - baseTrials;
    epochResp = {respVisual, respPlan, respExec};

    respChans = find(screening(s).responsive);
    genMat = nan(3,3);
    for iTrain = 1:3
        for iTest = 1:3
            genMat(iTrain, iTest) = crossEpochLOOAcc(epochResp{iTrain}(respChans,:), epochResp{iTest}(respChans,:), dirGroup);
        end
    end

    genResults(s).Day = S.Day;
    genResults(s).RunN = S.RunN;
    genResults(s).epochNames = epochNames;
    genResults(s).genMat = genMat; % rows = train epoch, cols = test epoch
    genResults(s).nResp = numel(respChans);

    fprintf('%s run-%03d: temporal generalization matrix (rows=train, cols=test, %%):\n', S.Day, S.RunN);
    fprintf('%12s%10s%10s%10s\n', '', epochNames{:});
    for iTrain = 1:3
        fprintf('%12s%10.1f%10.1f%10.1f\n', epochNames{iTrain}, 100*genMat(iTrain,:));
    end

    plotGenMatrix(S, genMat, epochNames, opts);
end

outFile = fullfile(fileparts(opts.figDir), 'decode_generalization.mat');
save(outFile, 'genResults');
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
% Leave-one-trial-out: standardization + centroids built from respTrain
% (excluding trial t); trial t is classified using its respTest value,
% standardized with respTrain's (t-excluded) mean/SD. Reduces exactly to
% decodeDirection.m's looDecodeAcc when respTrain and respTest are the
% same matrix. See decodeDirection.m for the closed-form leave-one-out
% update this is based on.
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

    ztest = (xTest - mu) ./ sd; % standardize the TEST-epoch point with the TRAIN epoch's own stats
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
function plotGenMatrix(S, genMat, epochNames, opts)
f = figure('Visible','off','Position',[0 0 500 460], 'Color', 'w');
ax = axes(f);
imagesc(ax, 1:3, 1:3, genMat*100);
set(ax, 'CLim', [0 max(60,max(genMat(:)*100))], 'YDir', 'reverse', 'Color', 'w');
colormap(ax, parula);
cb = colorbar(ax); cb.Label.String = 'decode accuracy (%)';
set(ax, 'XTick', 1:3, 'XTickLabel', epochNames, 'YTick', 1:3, 'YTickLabel', epochNames);
xlabel(ax, 'Test epoch'); ylabel(ax, 'Train epoch');
for i = 1:3
    for j = 1:3
        val = genMat(i,j)*100;
        txtCol = 'k'; if val > 35, txtCol = 'w'; end
        text(ax, j, i, sprintf('%.1f%%', val), 'HorizontalAlignment', 'center', 'Color', txtCol, 'FontWeight', 'bold');
    end
end
title(ax, sprintf('%s run-%03d: temporal generalization matrix\n(nearest-centroid, LOO-CV, chance=12.5%%, full responsive population)', S.Day, S.RunN), 'FontSize', 9);
forceLightTheme(f);
saveas(f, fullfile(opts.figDir, sprintf('%s_run-%03d_decode_genmatrix.png', S.Day, S.RunN)));
close(f);
end
