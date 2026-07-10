function artifactCheck = checkArtifactCorrelation(sessions, opts)
% CHECKARTIFACTCORRELATION Test whether the widespread, near-simultaneous
% MUA response at the go-cue (flagged as a possible shared movement/EMG
% artifact in PIPELINE_REPORT.md sec 6) looks like a single shared signal
% riding on every channel, or like spatially-organized neural activity.
%
% Logic: take each channel's single-trial peri-saccadic response, remove
% the direction-tuning contribution (subtract each direction-condition's
% own mean), leaving trial-to-trial fluctuation ("noise") that should be
% roughly independent across channels if responses are genuinely local
% and neural. Compute the pairwise correlation of this residual across
% all channel pairs. A shared artifact (common-mode EMG/reference/motion
% signal added to every channel) predicts uniformly high correlation
% that does NOT fall off with physical distance between channels. Locally
% organized neural activity predicts correlation that decays with
% distance (nearby channels share local circuitry/volume conduction,
% far-apart channels do not).
%
% artifactCheck = checkArtifactCorrelation(sessions, opts)
%   sessions : struct array from loadAllSessions()
%   opts.respWindow : [t0 t1] s relative to go-cue (default [0 0.3], same
%                     as computeDirectionTuning's tuning window)
%   opts.baseWindow : baseline window (default [-0.9 -0.6])
%   opts.figDir     : figure output dir

if nargin < 2, opts = struct(); end
opts = setDefault(opts, 'respWindow', [0 0.3]);
opts = setDefault(opts, 'baseWindow', [-0.9 -0.6]);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));
if ~isfolder(opts.figDir), mkdir(opts.figDir); end

artifactCheck = struct([]);

for s = 1:numel(sessions)
    S = sessions(s);
    isCorrect = strcmp(S.trialTable.outcome, 'correct') & ~isnan(S.trialTable.direction_deg);
    respIdx = S.tb >= opts.respWindow(1) & S.tb < opts.respWindow(2);
    baseIdx = S.tb >= opts.baseWindow(1) & S.tb < opts.baseWindow(2);

    respTrials = squeeze(mean(S.MUA(:, isCorrect, respIdx), 3)); % nChan x nTrial
    baseTrials = squeeze(mean(S.MUA(:, isCorrect, baseIdx), 3));
    d = respTrials - baseTrials; % nChan x nTrial, baseline-subtracted response

    dirDeg = S.trialTable.direction_deg(isCorrect);
    dirGroup = round(dirDeg/45) + 1;

    % remove each direction-condition's own mean per channel, leaving
    % trial-to-trial residual (removes genuine direction tuning so it
    % can't inflate the correlation between similarly-tuned channels)
    resid = nan(size(d));
    for g = 1:8
        idx = dirGroup == g;
        resid(:, idx) = d(:, idx) - mean(d(:, idx), 2);
    end

    R = corr(resid'); % nChan x nChan
    nChan = size(R,1);
    offDiagMask = ~eye(nChan);
    meanR = mean(R(offDiagMask));
    medianR = median(R(offDiagMask));

    depthUm = [];
    if isfield(S, 'CHdepthUm')
        depthUm = S.CHdepthUm;
    end

    % correlation vs. inter-channel physical distance
    distBins = []; corrByBin = [];
    if ~isempty(depthUm)
        distMat = abs(depthUm - depthUm');
        edges = [0 200 500 1000 2000 4000 8000 Inf];
        binIdx = discretize(distMat(offDiagMask), edges);
        rvals = R(offDiagMask);
        corrByBin = accumarray(binIdx(~isnan(binIdx)), rvals(~isnan(binIdx)), [numel(edges)-1 1], @mean, NaN);
        distBins = edges(1:end-1);
    end

    artifactCheck(s).Day = S.Day;
    artifactCheck(s).RunN = S.RunN;
    artifactCheck(s).meanOffDiagCorr = meanR;
    artifactCheck(s).medianOffDiagCorr = medianR;
    artifactCheck(s).distBinEdgesUm = distBins;
    artifactCheck(s).corrByDistBin = corrByBin;

    fprintf('%s run-%03d: mean pairwise trial-residual correlation = %.3f (median %.3f)\n', ...
        S.Day, S.RunN, meanR, medianR);

    plotArtifactCheck(S, R, depthUm, distBins, corrByBin, opts);
end

outFile = fullfile(fileparts(opts.figDir), 'artifact_correlation_check.mat');
save(outFile, 'artifactCheck');
fprintf('saved %s\n', outFile);

end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function plotArtifactCheck(S, R, depthUm, distBins, corrByBin, opts)
nChan = size(R,1);
f = figure('Visible','off','Position',[0 0 1400 550], 'Color', 'w');
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl);
imagesc(ax1, 1:nChan, 1:nChan, R);
set(ax1, 'CLim', [-1 1], 'YDir', 'normal', 'Color', 'w');
colormap(ax1, redblueLocal(256));
colorbar(ax1);
axis(ax1, 'square');
xlabel(ax1, 'Channel #'); ylabel(ax1, 'Channel #');
title(ax1, 'Trial-residual correlation matrix');

ax2 = nexttile(tl);
if ~isempty(depthUm)
    [~, order] = sort(depthUm);
    imagesc(ax2, 1:nChan, 1:nChan, R(order,order));
    set(ax2, 'CLim', [-1 1], 'YDir', 'normal', 'Color', 'w');
    colormap(ax2, redblueLocal(256));
    colorbar(ax2);
    axis(ax2, 'square');
    xlabel(ax2, 'Channel (sorted by real depth)'); ylabel(ax2, 'Channel (sorted by real depth)');
else
    axis(ax2, 'off');
    text(ax2, 0.5, 0.5, 'no depth info', 'HorizontalAlignment', 'center');
end
title(ax2, 'Same, sorted by depth');

ax3 = nexttile(tl);
if ~isempty(corrByBin)
    binCenters = distBins/1000; % mm, approx (left edges)
    plot(ax3, binCenters, corrByBin, '-o', 'LineWidth', 1.5);
    set(ax3, 'Color', 'w');
    xlabel(ax3, 'Inter-channel distance bin (mm, left edge)');
    ylabel(ax3, 'Mean trial-residual correlation');
    title(ax3, 'Correlation vs. distance');
    yline(ax3, 0, 'k:');
else
    axis(ax3, 'off');
end

title(tl, sprintf('%s run-%03d: shared-artifact check (uniform-with-distance = artifact-like; decaying = neural-like)', ...
    S.Day, S.RunN));
forceLightTheme(f);
saveas(f, figSavePath(opts.figDir, 'artifact_check', S.Day, S.RunN));
close(f);
end

% ------------------------------------------------------------------
function cmap = redblueLocal(n)
half = floor(n/2);
top = [linspace(0.1,1,half)', linspace(0.1,1,half)', ones(half,1)];
bot = [ones(n-half,1), linspace(1,0.1,n-half)', linspace(1,0.1,n-half)'];
cmap = [top; bot];
end
