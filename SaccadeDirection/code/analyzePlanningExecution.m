function results = analyzePlanningExecution(sessions, eyeSummary, opts)
% ANALYZEPLANNINGEXECUTION Splits saccade-direction tuning into THREE
% epochs of interest instead of one blended [0,0.3s] go-cue window:
%
%   'visual'    -- the response to the target/stimulus FLASH itself,
%                  before the delay period and well before any saccade
%                  (same window as screenResponsiveChannels.m's "stim"
%                  epoch, read from this session's own Stm timing:
%                  [-(STIMT+DELAYT), -DELAYT] ms relative to go-cue).
%                  Tests whether a channel is tuned to WHERE the target
%                  appeared (retinotopic/visual), independent of anything
%                  to do with the eventual saccade.
%   'planning'  -- premotor activity, from the go-cue instruction until
%                  the eye actually starts moving.
%   'execution' -- perimovement activity, during and immediately after
%                  the saccade itself.
%
% planning/execution are computed two ways ("mode"):
%   'gocue'    -- simple fixed split relative to the go-cue (t=0):
%                 planning=[0,0.1]s, execution=[0.1,0.3]s. Easy, but the
%                 cutoff doesn't mean the same thing in every session:
%                 median saccade latency ranges 68-218ms across the 12
%                 sessions (PIPELINE_REPORT sec 5.1).
%   'saclocked' -- anchored to each TRIAL's own detected saccade onset
%                 (from analyzeEyeData.m/eyeSummary), the methodologically
%                 correct way to do this split: planning=[-0.15,-0.02]s
%                 and execution=[-0.02,+0.15]s relative to onset, not
%                 go-cue. Trials with no detected saccade are excluded
%                 from this mode (and, for a fair 3-way comparison, from
%                 the visual epoch too within this mode's output).
% The visual epoch is identical in both modes (its timing never depends
% on saccade onset), but is restricted to the same trial subset as
% planning/execution within each mode so the three epochs are directly
% comparable in N.
%
% For each epoch and each channel: responsiveness (vs. the same
% [-0.9,-0.6]s pre-stimulus baseline used elsewhere) via Wilcoxon
% signed-rank + BH-FDR + minimum effect size (matching
% screenResponsiveChannels.m), and direction tuning within that epoch for
% channels responsive in it (vector-sum + permutation test, plus a
% Kruskal-Wallis omnibus test -- classifyTuning.m -- matching
% computeDirectionTuning.m).
%
% Does NOT modify or overwrite anything from the main pipeline
% (screening_results.mat, tuning_results.mat, *_overview.png etc) --
% saves to planning_execution_results.mat and *_planexec_<mode>.png /
% *_planexec_<mode>_curves.png, each panel/title labeled with the actual
% ms window used for that epoch.
%
% results = analyzePlanningExecution(sessions, eyeSummary, opts)
%   sessions   : struct array from loadAllSessions()
%   eyeSummary : struct array from analyzeEyeData() (for saccade onsets)
%   opts.baseWindow    : default [-0.9 -0.6]
%   opts.goCueWindows  : default planning=[0 0.1], execution=[0.1 0.3]
%   opts.sacWindows    : default planning=[-0.15 -0.02], execution=[-0.02 0.15]
%   opts.alpha, opts.minEffect, opts.nPerm : same defaults as elsewhere
%   opts.figDir : default <saccadeDataDir()>/figures

if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'baseWindow', [-0.9 -0.6]);
opts = setDefault(opts, 'visualWindow', [-0.5 -0.35]); % matches screenResponsiveChannels.m's "visual" window exactly, so the same epoch shows the same ms range in every figure
opts = setDefault(opts, 'goCueWindows', struct('planning', [0 0.1], 'execution', [0.1 0.3]));
opts = setDefault(opts, 'sacWindows', struct('planning', [-0.15 -0.02], 'execution', [-0.02 0.15]));
opts = setDefault(opts, 'alpha', 0.05);
opts = setDefault(opts, 'minEffect', 0.3);
opts = setDefault(opts, 'nPerm', 1000);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));
if ~isfolder(opts.figDir), mkdir(opts.figDir); end

dirsDeg = 0:45:315;
results = struct([]);
rowIdx = 0;

for s = 1:numel(sessions)
    S = sessions(s);
    esIdx = find(strcmp({eyeSummary.Day}, S.Day) & [eyeSummary.RunN] == S.RunN);
    isCorrect = strcmp(S.trialTable.outcome, 'correct') & ~isnan(S.trialTable.direction_deg);
    dirDeg = S.trialTable.direction_deg;

    baseIdx = S.tb >= opts.baseWindow(1) & S.tb < opts.baseWindow(2);
    baseTrials = squeeze(mean(S.MUA(:, isCorrect, baseIdx), 3)); % nChan x nTrial(correct)

    % Visual epoch: the stimulus-flash window, identical to
    % screenResponsiveChannels.m's "visual" window and fixed relative to
    % go-cue regardless of mode.
    visWin = opts.visualWindow;
    visIdx = S.tb >= visWin(1) & S.tb < visWin(2);
    respVisualAll = squeeze(mean(S.MUA(:, isCorrect, visIdx), 3)) - baseTrials;

    % --- mode 1: fixed split relative to go-cue (same window for all trials)
    planIdxGC = S.tb >= opts.goCueWindows.planning(1) & S.tb < opts.goCueWindows.planning(2);
    execIdxGC = S.tb >= opts.goCueWindows.execution(1) & S.tb < opts.goCueWindows.execution(2);
    respPlanGC = squeeze(mean(S.MUA(:, isCorrect, planIdxGC), 3)) - baseTrials;
    respExecGC = squeeze(mean(S.MUA(:, isCorrect, execIdxGC), 3)) - baseTrials;

    winsGC = struct('visual', visWin, 'planning', opts.goCueWindows.planning, 'execution', opts.goCueWindows.execution);
    r = threeMode(respVisualAll, respPlanGC, respExecGC, dirDeg(isCorrect), dirsDeg, opts);
    r.Day = S.Day; r.RunN = S.RunN; r.mode = 'gocue'; r.windows = winsGC;
    rowIdx = rowIdx+1; results(rowIdx).data = r;
    plotPlanExec(S, r, opts, 'gocue');
    plotPlanExecCurves(S, r, dirsDeg, opts, 'gocue');

    % --- mode 2: saccade-onset-locked (variable per-trial center)
    if ~isempty(esIdx) && S.eyeCalib.valid
        es = eyeSummary(esIdx);
        onsetTime = es.onsetTime; % nTrial x 1, NaN where no saccade detected
        correctIdxList = find(isCorrect);
        hasOnset = ~isnan(onsetTime(correctIdxList));
        useTrials = correctIdxList(hasOnset);

        nChan = size(S.MUA,1);
        respPlanSL = nan(nChan, numel(useTrials));
        respExecSL = nan(nChan, numel(useTrials));
        for k = 1:numel(useTrials)
            t = useTrials(k);
            winP = S.tb >= (onsetTime(t)+opts.sacWindows.planning(1)) & S.tb < (onsetTime(t)+opts.sacWindows.planning(2));
            winE = S.tb >= (onsetTime(t)+opts.sacWindows.execution(1)) & S.tb < (onsetTime(t)+opts.sacWindows.execution(2));
            if any(winP), respPlanSL(:,k) = mean(S.MUA(:,t,winP),3); end
            if any(winE), respExecSL(:,k) = mean(S.MUA(:,t,winE),3); end
        end
        baseSL = baseTrials(:, hasOnset);
        respPlanSL = respPlanSL - baseSL;
        respExecSL = respExecSL - baseSL;
        respVisualSL = respVisualAll(:, hasOnset); % same trial subset, for a fair 3-way comparison
        dirDegSL = dirDeg(useTrials);

        winsSL = struct('visual', visWin, 'planning', opts.sacWindows.planning, 'execution', opts.sacWindows.execution);
        r2 = threeMode(respVisualSL, respPlanSL, respExecSL, dirDegSL, dirsDeg, opts);
        r2.Day = S.Day; r2.RunN = S.RunN; r2.mode = 'saclocked'; r2.windows = winsSL;
        r2.nTrialsUsed = numel(useTrials); r2.nTrialsTotal = numel(correctIdxList);
        rowIdx = rowIdx+1; results(rowIdx).data = r2;
        plotPlanExec(S, r2, opts, 'saclocked');
        plotPlanExecCurves(S, r2, dirsDeg, opts, 'saclocked');
        fprintf('%s run-%03d [saclocked]: %d/%d correct trials had a detected saccade\n', ...
            S.Day, S.RunN, numel(useTrials), numel(correctIdxList));
    else
        fprintf('%s run-%03d [saclocked]: skipped, no valid eye calibration/onset data\n', S.Day, S.RunN);
    end
end

outFile = fullfile(fileparts(opts.figDir), 'planning_execution_results.mat');
save(outFile, 'results');
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
q = nan(size(p));
ok = ~isnan(p);
n = sum(ok);
if n == 0, return; end
[ps, order] = sort(p(ok));
ranks = (1:n)';
qs = ps .* n ./ ranks;
qs = flipud(cummin(flipud(qs)));
qs = min(qs, 1);
qtmp = nan(n,1);
qtmp(order) = qs;
q(ok) = qtmp;
end

% ------------------------------------------------------------------
function tuning = fitTuning(resp, dirDeg, dirsDeg, nPerm)
% resp: 1 x nTrial (one channel's baseline-subtracted window response); returns struct with prefDir/resultantLen/pPerm/tuned
% NB: force row vectors -- dirDeg often arrives as a column (straight from
% a table variable) while resp is a row; combining two logicals of
% mismatched orientation via `&` silently broadcasts to an NxN matrix
% instead of an elementwise op, which then blows up indexing resp with it.
resp = resp(:)';
dirDeg = dirDeg(:)';
dirGroup = round(dirDeg/45) + 1;
meanByDir = nan(1,8);
for d = 1:8
    meanByDir(d) = mean(resp(dirGroup==d), 'omitnan');
end
w = max(meanByDir, 0);
if sum(w) > 0
    vx = sum(w.*cosd(dirsDeg))/sum(w);
    vy = sum(w.*sind(dirsDeg))/sum(w);
else
    vx = 0; vy = 0;
end
prefDir = mod(atan2d(vy,vx), 360);
resultantLen = hypot(vx,vy);

ok = ~isnan(resp);
permR = nan(nPerm,1);
for p = 1:nPerm
    permGroup = dirGroup(randperm(numel(dirGroup)));
    mBy = nan(1,8);
    for d = 1:8
        mBy(d) = mean(resp(ok & permGroup==d), 'omitnan');
    end
    wp = max(mBy,0);
    if sum(wp) > 0
        pvx = sum(wp.*cosd(dirsDeg))/sum(wp);
        pvy = sum(wp.*sind(dirsDeg))/sum(wp);
        permR(p) = hypot(pvx,pvy);
    else
        permR(p) = 0;
    end
end
pPerm = (1 + sum(permR >= resultantLen)) / (nPerm + 1);

% Kruskal-Wallis omnibus test across the 8 direction groups, alongside the
% vector-sum permutation test -- see classifyTuning.m for why both are
% needed (pPerm alone misses bimodal/complex non-unimodal patterns).
okKW = ~isnan(resp);
if sum(okKW) > 8 && numel(unique(dirGroup(okKW))) > 1
    pKW = kruskalwallis(resp(okKW)', dirGroup(okKW)', 'off');
else
    pKW = nan;
end

tuning.prefDir = prefDir;
tuning.resultantLen = resultantLen;
tuning.pPerm = pPerm;
tuning.pKW = pKW;
tuning.tuned = pPerm < 0.05;
tuning.meanByDir = meanByDir;
end

% ------------------------------------------------------------------
function [respFlag, eff, tuning] = screenAndTuneOneEpoch(resp, dirDeg, dirsDeg, opts)
% Responsiveness (signed-rank + FDR + effect size) and, for responsive
% channels, direction tuning, for ONE epoch's nChan x nTrial response
% matrix. Shared by all three epochs (visual/planning/execution) so the
% same criteria are applied uniformly.
nChan = size(resp,1);
baseSD = std(resp, 0, 2, 'omitnan');
baseSD(baseSD==0 | isnan(baseSD)) = eps;

p = nan(nChan,1); eff = nan(nChan,1);
for c = 1:nChan
    d = resp(c,:); d = d(~isnan(d));
    if numel(d) >= 5 && ~all(d==d(1))
        p(c) = signrank(d);
        eff(c) = mean(d) / baseSD(c);
    end
end
fdrP = bhFDR(p);
respFlag = fdrP < opts.alpha & abs(eff) > opts.minEffect;

emptyTuning = struct('prefDir',nan,'resultantLen',nan,'pPerm',nan,'pKW',nan,'tuned',false,'meanByDir',nan(1,8));
tuning = repmat(emptyTuning, 1, nChan);
for c = 1:nChan
    if respFlag(c)
        tuning(c) = fitTuning(resp(c,:), dirDeg, dirsDeg, opts.nPerm);
    end
end
end

% ------------------------------------------------------------------
function r = threeMode(respVisual, respPlan, respExec, dirDeg, dirsDeg, opts)
[respVisualFlag, effVisual, tuningVisual] = screenAndTuneOneEpoch(respVisual, dirDeg, dirsDeg, opts);
[respPlanFlag, effPlan, tuningPlan] = screenAndTuneOneEpoch(respPlan, dirDeg, dirsDeg, opts);
[respExecFlag, effExec, tuningExec] = screenAndTuneOneEpoch(respExec, dirDeg, dirsDeg, opts);

r = struct();
r.effVisual = effVisual; r.effPlan = effPlan; r.effExec = effExec;
r.respVisualFlag = respVisualFlag; r.respPlanFlag = respPlanFlag; r.respExecFlag = respExecFlag;
r.tuningVisual = tuningVisual; r.tuningPlan = tuningPlan; r.tuningExec = tuningExec;
end

% ------------------------------------------------------------------
function lbl = winLabel(name, winSec)
lbl = sprintf('%s [%.0f,%.0f]ms', name, winSec(1)*1000, winSec(2)*1000);
end

% ------------------------------------------------------------------
function plotPlanExec(S, r, opts, modeName) %#ok<INUSL>
nChan = numel(r.effPlan);
chIdx = (1:nChan)';

f = figure('Visible','off','Position',[0 0 1900 600], 'Color', 'w');
tl = tiledlayout(f, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: effect map, visual/planning/execution, per channel (depth order, deep=bottom)
ax1 = nexttile(tl);
M = [r.effVisual, r.effPlan, r.effExec];
M(~[r.respVisualFlag, r.respPlanFlag, r.respExecFlag]) = nan;
imagesc(ax1, 1:3, chIdx, M, 'AlphaData', double(~isnan(M)));
set(ax1, 'CLim', [-2 2], 'YDir', 'normal', 'Color', [0.85 0.85 0.85]);
colormap(ax1, redblueLocal(256));
set(ax1, 'XTick', 1:3, 'XTickLabel', {winLabel('visual',r.windows.visual), winLabel('planning',r.windows.planning), winLabel('execution',r.windows.execution)}, ...
    'XTickLabelRotation', 20, 'FontSize', 7);
colorbar(ax1);
ylabel(ax1, 'Channel # (1 = deepest)');
title(ax1, sprintf('%d/%d visual, %d/%d planning, %d/%d execution', ...
    sum(r.respVisualFlag), nChan, sum(r.respPlanFlag), nChan, sum(r.respExecFlag), nChan));

% Panel 2: tuning strength, visual vs execution (does retinotopic/visual
% tuning predict eventual motor tuning?) -- colored by category
% (classifyTuning.m: red=directional, orange=complex, gray=untuned)
ax2 = nexttile(tl); hold(ax2, 'on');
Rvis = [r.tuningVisual.resultantLen]; Rexec = [r.tuningExec.resultantLen];
[catVis, colVis] = classifyTuning([r.tuningVisual.pPerm], [r.tuningVisual.pKW]);
[catExec, colExec] = classifyTuning([r.tuningExec.pPerm], [r.tuningExec.pKW]);
isDirVis = strcmp(catVis, 'directional');
isDirExec = strcmp(catExec, 'directional');
useCol2 = colExec; useCol2(isDirVis & ~isDirExec, :) = colVis(isDirVis & ~isDirExec, :);
scatter(ax2, Rvis, Rexec, 14, useCol2, 'filled');
xlim(ax2, [0 1]); ylim(ax2, [0 1]);
plot(ax2, [0 1], [0 1], 'k:');
set(ax2, 'Color', 'w');
xlabel(ax2, 'Tuning strength, visual'); ylabel(ax2, 'Tuning strength, execution');
title(ax2, 'Tuning strength: visual vs execution');

% Panel 3: tuning strength, planning vs execution
ax3 = nexttile(tl); hold(ax3, 'on');
Rplan = [r.tuningPlan.resultantLen];
[catPlan, colPlan] = classifyTuning([r.tuningPlan.pPerm], [r.tuningPlan.pKW]);
isDirPlan = strcmp(catPlan, 'directional');
useCol3 = colExec; useCol3(isDirPlan & ~isDirExec, :) = colPlan(isDirPlan & ~isDirExec, :);
scatter(ax3, Rplan, Rexec, 14, useCol3, 'filled');
xlim(ax3, [0 1]); ylim(ax3, [0 1]);
plot(ax3, [0 1], [0 1], 'k:');
set(ax3, 'Color', 'w');
xlabel(ax3, 'Tuning strength, planning'); ylabel(ax3, 'Tuning strength, execution');
title(ax3, 'Tuning strength: planning vs execution');

% Panel 4: preferred direction, planning vs execution (only channels tuned in BOTH)
ax4 = nexttile(tl); hold(ax4, 'on');
tunedBoth = [r.tuningPlan.tuned] & [r.tuningExec.tuned];
if any(tunedBoth)
    PDplan = [r.tuningPlan.prefDir]; PDexec = [r.tuningExec.prefDir];
    scatter(ax4, PDplan(tunedBoth), PDexec(tunedBoth), 20, [0.1 0.4 0.85], 'filled');
end
xlim(ax4, [0 360]); ylim(ax4, [0 360]);
plot(ax4, [0 360], [0 360], 'k:');
set(ax4, 'Color', 'w', 'XTick', 0:90:360, 'YTick', 0:90:360);
xlabel(ax4, 'Preferred direction, planning (deg)'); ylabel(ax4, 'Preferred direction, execution (deg)');
title(ax4, sprintf('PD planning vs execution (n=%d tuned in both)', sum(tunedBoth)));

title(tl, sprintf('%s run-%03d: visual / planning / execution (%s-anchored)', S.Day, S.RunN, modeName));
forceLightTheme(f);
saveas(f, figSavePath(opts.figDir, sprintf('planexec_%s', modeName), S.Day, S.RunN));
close(f);
end

% ------------------------------------------------------------------
function plotPlanExecCurves(S, r, dirsDeg, opts, modeName) %#ok<INUSL>
% Per-channel small-multiple tuning curves, visual/planning/execution
% stacked for the SAME channel (depth order), so you can see directly
% whether/how a channel's tuning shape changes across the three epochs.
nChan = numel(r.effPlan);
respAny = r.respVisualFlag | r.respPlanFlag | r.respExecFlag;
chans = find(respAny);
nShow = numel(chans);
if nShow == 0, return; end

[catVis, colVis] = classifyTuning([r.tuningVisual.pPerm], [r.tuningVisual.pKW]);
[catPlan, colPlan] = classifyTuning([r.tuningPlan.pPerm], [r.tuningPlan.pKW]);
[catExec, colExec] = classifyTuning([r.tuningExec.pPerm], [r.tuningExec.pKW]);

f = figure('Visible','off','Position',[0 0 1600 max(700,110*ceil(nShow/6))], 'Color', 'w');
nCol = 6; nRow = ceil(nShow/nCol)*3; % 3 rows per channel (visual/planning/execution)
tl = tiledlayout(f, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nShow
    c = chans(i);
    rowBlock = floor((i-1)/nCol);
    colPos = mod(i-1,nCol);
    tileVis = rowBlock*3*nCol + colPos + 1;
    tilePlan = tileVis + nCol;
    tileExec = tilePlan + nCol;

    axV = nexttile(tl, tileVis);
    vV = r.tuningVisual(c).meanByDir;
    plot(axV, [dirsDeg dirsDeg(1)+360], [vV vV(1)], '-o', 'Color', colVis(i,:), ...
        'MarkerFaceColor', colVis(i,:), 'LineWidth', 1.2, 'MarkerSize', 3);
    set(axV, 'XTick', [], 'Color', 'w');
    title(axV, sprintf('ch%d: visual (%s)', c, catVis{i}), 'FontSize', 6, 'FontWeight', 'normal');

    axP = nexttile(tl, tilePlan);
    vP = r.tuningPlan(c).meanByDir;
    plot(axP, [dirsDeg dirsDeg(1)+360], [vP vP(1)], '-o', 'Color', colPlan(i,:), ...
        'MarkerFaceColor', colPlan(i,:), 'LineWidth', 1.2, 'MarkerSize', 3);
    set(axP, 'XTick', [], 'Color', 'w');
    title(axP, sprintf('ch%d: plan (%s)', c, catPlan{i}), 'FontSize', 6, 'FontWeight', 'normal');

    axE = nexttile(tl, tileExec);
    vE = r.tuningExec(c).meanByDir;
    plot(axE, [dirsDeg dirsDeg(1)+360], [vE vE(1)], '-o', 'Color', colExec(i,:), ...
        'MarkerFaceColor', colExec(i,:), 'LineWidth', 1.2, 'MarkerSize', 3);
    set(axE, 'XTick', [], 'Color', 'w');
    title(axE, sprintf('ch%d: exec (%s)', c, catExec{i}), 'FontSize', 6, 'FontWeight', 'normal');
end
title(tl, sprintf(['%s run-%03d: visual (top) / planning (mid) / execution (bottom) tuning shape, %s-anchored ' ...
    '-- %s | %s | %s -- red=directional, orange=complex, gray=untuned'], S.Day, S.RunN, modeName, ...
    winLabel('visual',r.windows.visual), winLabel('planning',r.windows.planning), winLabel('execution',r.windows.execution)), ...
    'FontSize', 10);
forceLightTheme(f);
saveas(f, figSavePath(opts.figDir, sprintf('planexec_%s_curves', modeName), S.Day, S.RunN));
close(f);
end

% ------------------------------------------------------------------
function cmap = redblueLocal(n)
half = floor(n/2);
top = [linspace(0.1,1,half)', linspace(0.1,1,half)', ones(half,1)];
bot = [ones(n-half,1), linspace(1,0.1,n-half)', linspace(1,0.1,n-half)'];
cmap = [top; bot];
end
