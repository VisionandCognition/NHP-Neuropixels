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
% (PIPELINE_REPORT.md sec 7, open item 2).
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
    pc1Score = V(:,1) * Sig(1,1); % per-trial score, 1 x nTrial (as column here nTrial x 1)
    pc1Loading = U(:,1); % per-channel loading
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

    sharedCheck(s).Day = S.Day;
    sharedCheck(s).RunN = S.RunN;
    sharedCheck(s).varExplainedPC1 = varExplained;
    sharedCheck(s).loadingSignConsistency = loadingSignConsistency;
    sharedCheck(s).rVel = rVel; sharedCheck(s).pVel = pVel;
    sharedCheck(s).rAmp = rAmp; sharedCheck(s).pAmp = pAmp;
    sharedCheck(s).nMatched = nMatched;

    fprintf('%s run-%03d: PC1 explains %.1f%% of residual variance (%.0f%% of channels agree in loading sign); PC1 vs peakVel r=%.3f (p=%.3f), vs amplitude r=%.3f (p=%.3f), n=%d\n', ...
        S.Day, S.RunN, 100*varExplained, 100*loadingSignConsistency, rVel, pVel, rAmp, pAmp, nMatched);

    plotSharedCheck(S, pc1Score, pc1Loading, es, correctIdxList, ok, varExplained, opts); %#ok<NASGU>
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
function plotSharedCheck(S, pc1Score, pc1Loading, es, correctIdxList, ok, varExplained, opts)
if isempty(es) || ~any(ok), return; end
peakVel = es.peakVel(correctIdxList(ok));
amplitude = es.amplitude(correctIdxList(ok));
pc1 = pc1Score(ok);

f = figure('Visible','off','Position',[0 0 1400 500], 'Color', 'w');
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

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

title(tl, sprintf('%s run-%03d: shared-component / kinematics check', S.Day, S.RunN));
forceLightTheme(f);
saveas(f, fullfile(opts.figDir, sprintf('%s_run-%03d_shared_kinematics.png', S.Day, S.RunN)));
close(f);
end
