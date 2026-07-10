function sharedCheck = checkSharedComponentKinematics(sessions, eyeSummary, opts)
% CHECKSHAREDCOMPONENTKINEMATICS Follow-up to checkArtifactCorrelation.m:
% that check found trial-residual correlation decays with distance in
% every session (neural-like) but plateaus at a small nonzero value even
% between far-apart channels (consistent with a smaller shared/common-mode
% component riding underneath genuine local structure). This tests
% directly whether that shared component tracks saccade KINEMATICS
% (peak velocity, amplitude) -- if it does, in a way that's the same
% regardless of which channels you look at, that's the signature of
% EMG/movement contamination rather than e.g. a shared reference issue
% (PIPELINE_REPORT.md sec 7, open item 1).
%
% Method: for each session, take the same baseline-subtracted,
% direction-demeaned trial-residual matrix used in
% checkArtifactCorrelation.m (nChan x nTrial), and extract its dominant
% shared fluctuation via the first principal component (SVD) across
% channels -- this is what "shared component" means concretely: the one
% linear combination of channels that explains the most trial-to-trial
% covariance, i.e. the closest thing to a single common-mode signal. Its
% per-trial score is then correlated (Spearman) against that trial's
% saccade peak velocity and amplitude (from analyzeEyeData.m).
%
% For sessions where this correlation is significant with BOTH peak
% velocity and amplitude (same sign -- the clean EMG/movement signature),
% this ALSO performs common-mode removal (subtracting the rank-1 PC1
% reconstruction, loading*score', from the residual) and recomputes the
% trial-residual-correlation-vs-distance relationship
% (checkArtifactCorrelation.m's own test) on the cleaned data, plotting
% original vs. common-mode-removed side by side: if the local (short
% distance) correlation survives while the long-distance plateau drops,
% that confirms the plateau specifically was the shared/kinematic
% component, and the depth-organized structure elsewhere in this pipeline
% is not an artifact of it.
%
% sharedCheck = checkSharedComponentKinematics(sessions, eyeSummary, opts)
%   sessions   : struct array from loadAllSessions()
%   eyeSummary : struct array from analyzeEyeData()
%   opts.respWindow, opts.baseWindow : same as checkArtifactCorrelation.m defaults
%   opts.figDir : default <saccadeDataDir()>/figures

if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'respWindow', [0 0.3]);
opts = setDefault(opts, 'baseWindow', [-0.9 -0.6]);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));
if ~isfolder(opts.figDir), mkdir(opts.figDir); end

sharedCheck = struct([]);

for s = 1:numel(sessions)
    S = sessions(s);
    esIdx = find(strcmp({eyeSummary.Day}, S.Day) & [eyeSummary.RunN] == S.RunN);
    isCorrect = strcmp(S.trialTable.outcome, 'correct') & ~isnan(S.trialTable.direction_deg);
    respIdx = S.tb >= opts.respWindow(1) & S.tb < opts.respWindow(2);
    baseIdx = S.tb >= opts.baseWindow(1) & S.tb < opts.baseWindow(2);

    respTrials = squeeze(mean(S.MUA(:, isCorrect, respIdx), 3)); % nChan x nTrial
    baseTrials = squeeze(mean(S.MUA(:, isCorrect, baseIdx), 3));
    d = respTrials - baseTrials;

    dirDeg = S.trialTable.direction_deg(isCorrect);
    dirGroup = round(dirDeg/45) + 1;
    resid = nan(size(d));
    for g = 1:8
        idx = dirGroup == g;
        resid(:, idx) = d(:, idx) - mean(d(:, idx), 2);
    end

    % dominant shared component: PC1 of the (channels x trials) residual,
    % via SVD on the channel-centered matrix (each channel already
    % baseline/direction demeaned above, so no further centering needed)
    [U, Sig, V] = svd(resid, 'econ');
    pc1Score = V(:,1) * Sig(1,1); % per-trial score, nTrial x 1
    pc1Loading = U(:,1); % per-channel loading, nChan x 1
    varExplained = Sig(1,1)^2 / sum(diag(Sig).^2);
    loadingSignConsistency = max(mean(pc1Loading>0), mean(pc1Loading<0)); % fraction of channels agreeing in sign

    rVel = nan; pVel = nan; rAmp = nan; pAmp = nan;
    nMatched = 0;
    es = []; correctIdxList = []; ok = [];
    if ~isempty(esIdx)
        es = eyeSummary(esIdx);
        correctIdxList = find(isCorrect);
        peakVel = es.peakVel(correctIdxList);
        amplitude = es.amplitude(correctIdxList);
        ok = ~isnan(peakVel) & ~isnan(amplitude);
        nMatched = sum(ok);
        if nMatched >= 8
            [rVel, pVel] = corr(pc1Score(ok), peakVel(ok), 'Type', 'Spearman');
            [rAmp, pAmp] = corr(pc1Score(ok), amplitude(ok), 'Type', 'Spearman');
        end
    end

    isKinematicSig = ~isnan(pVel) && ~isnan(pAmp) && pVel < 0.05 && pAmp < 0.05 && sign(rVel) == sign(rAmp);

    % Common-mode removal: subtract the rank-1 PC1 reconstruction, only
    % where the kinematic link is confirmed above (no point doing this for
    % sessions where PC1 isn't shown to be kinematics-linked in the first
    % place -- see docstring).
    distEdgesUm = [0 200 500 1000 2000 4000 8000 Inf];
    corrByDistOrig = []; corrByDistClean = [];
    if isKinematicSig && isfield(S, 'CHdepthUm')
        residClean = resid - pc1Loading * pc1Score';
        corrByDistOrig = corrVsDistance(resid, S.CHdepthUm, distEdgesUm);
        corrByDistClean = corrVsDistance(residClean, S.CHdepthUm, distEdgesUm);
    end

    sharedCheck(s).Day = S.Day;
    sharedCheck(s).RunN = S.RunN;
    sharedCheck(s).varExplainedPC1 = varExplained;
    sharedCheck(s).loadingSignConsistency = loadingSignConsistency;
    sharedCheck(s).rVel = rVel; sharedCheck(s).pVel = pVel;
    sharedCheck(s).rAmp = rAmp; sharedCheck(s).pAmp = pAmp;
    sharedCheck(s).nMatched = nMatched;
    sharedCheck(s).isKinematicSig = isKinematicSig;
    sharedCheck(s).distEdgesUm = distEdgesUm;
    sharedCheck(s).corrByDistOrig = corrByDistOrig;
    sharedCheck(s).corrByDistClean = corrByDistClean;

    fprintf('%s run-%03d: PC1 explains %.1f%% of residual variance (%.0f%% of channels agree in loading sign); PC1 vs peakVel r=%.3f (p=%.3f), vs amplitude r=%.3f (p=%.3f), n=%d%s\n', ...
        S.Day, S.RunN, 100*varExplained, 100*loadingSignConsistency, rVel, pVel, rAmp, pAmp, nMatched, ...
        ternary(isKinematicSig, ' [kinematic-linked -> common-mode removal performed]', ''));

    plotSharedCheck(S, pc1Score, pc1Loading, es, correctIdxList, ok, varExplained, isKinematicSig, ...
        distEdgesUm, corrByDistOrig, corrByDistClean, opts);
end

outFile = fullfile(fileparts(opts.figDir), 'shared_component_kinematics.mat');
save(outFile, 'sharedCheck');
fprintf('saved %s\n', outFile);
end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

% ------------------------------------------------------------------
function corrByBin = corrVsDistance(resid, depthUm, distEdgesUm)
% Mean pairwise trial-residual correlation, binned by inter-channel real
% distance -- same method as checkArtifactCorrelation.m.
R = corr(resid');
nChan = size(R,1);
offDiagMask = ~eye(nChan);
distMat = abs(depthUm - depthUm');
binIdx = discretize(distMat(offDiagMask), distEdgesUm);
rvals = R(offDiagMask);
corrByBin = accumarray(binIdx(~isnan(binIdx)), rvals(~isnan(binIdx)), [numel(distEdgesUm)-1 1], @mean, NaN);
end

% ------------------------------------------------------------------
function plotSharedCheck(S, pc1Score, pc1Loading, es, correctIdxList, ok, varExplained, ...
    isKinematicSig, distEdgesUm, corrByDistOrig, corrByDistClean, opts)
if isempty(es) || ~any(ok), return; end
peakVel = es.peakVel(correctIdxList(ok));
amplitude = es.amplitude(correctIdxList(ok));
pc1 = pc1Score(ok);

f = figure('Visible','off','Position',[0 0 1800 500], 'Color', 'w');
tl = tiledlayout(f, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl); hold(ax1, 'on');
histogram(ax1, pc1Loading, 30, 'FaceColor', [0.3 0.3 0.7]);
set(ax1, 'Color', 'w');
xlabel(ax1, 'PC1 channel loading'); ylabel(ax1, '# channels');
title(ax1, sprintf('PC1 loadings (%.1f%% var explained)', 100*varExplained));

ax2 = nexttile(tl); hold(ax2, 'on');
scatter(ax2, pc1, peakVel, 14, [0.2 0.5 0.8], 'filled');
set(ax2, 'Color', 'w');
xlabel(ax2, 'PC1 trial score (shared component)'); ylabel(ax2, 'Saccade peak velocity (dva/s)');
title(ax2, 'Shared component vs. peak velocity');

ax3 = nexttile(tl); hold(ax3, 'on');
scatter(ax3, pc1, amplitude, 14, [0.8 0.4 0.1], 'filled');
set(ax3, 'Color', 'w');
xlabel(ax3, 'PC1 trial score (shared component)'); ylabel(ax3, 'Saccade amplitude (dva)');
title(ax3, 'Shared component vs. amplitude');

% Panel 4: common-mode removal result, only for sessions confirmed
% kinematic-linked (both peakVel and amplitude correlations significant,
% same sign) -- shows whether the local (short-distance) trial-residual
% correlation survives removing PC1 while the long-distance plateau drops.
ax4 = nexttile(tl); hold(ax4, 'on');
if isKinematicSig && ~isempty(corrByDistOrig)
    binCentersMm = distEdgesUm(1:end-1)/1000;
    hOrig = plot(ax4, binCentersMm, corrByDistOrig, '-o', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.5);
    hClean = plot(ax4, binCentersMm, corrByDistClean, '-o', 'Color', [0.1 0.5 0.1], 'LineWidth', 1.5);
    yline(ax4, 0, 'k:');
    % NB: 'Location','best' can silently fail to render in offscreen
    % ('Visible','off') figures -- use a fixed location instead.
    legend(ax4, [hOrig hClean], {'original','common-mode removed'}, 'Location', 'northeast', 'TextColor', 'k', 'Color', 'w');
    xlabel(ax4, 'Inter-channel distance bin (mm, left edge)');
    ylabel(ax4, 'Mean trial-residual correlation');
    title(ax4, 'Common-mode removal: before vs after');
else
    text(ax4, 0.5, 0.5, {'Not kinematic-linked','(no common-mode removal', 'performed for this session)'}, ...
        'HorizontalAlignment', 'center', 'Units', 'normalized', 'Color', [0.5 0.5 0.5]);
    axis(ax4, 'off');
end
set(ax4, 'Color', 'w');

title(tl, sprintf('%s run-%03d: shared-component / kinematics check', S.Day, S.RunN));
forceLightTheme(f);
saveas(f, figSavePath(opts.figDir, 'shared_kinematics', S.Day, S.RunN));
close(f);
end
