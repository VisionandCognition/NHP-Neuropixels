function results = analyzePlanningExecution(sessions, eyeSummary, opts)
% ANALYZEPLANNINGEXECUTION Split the go-cue-relative "peri-saccadic"
% window used elsewhere in this pipeline (screenResponsiveChannels.m,
% computeDirectionTuning.m) into a PLANNING (premotor) epoch and an
% EXECUTION (perimovement) epoch, two ways, per your question about
% whether this split is possible:
%
%   'gocue'    -- simple fixed split relative to the go-cue (t=0), same
%                 cutoff for every trial: planning=[0,0.1]s,
%                 execution=[0.1,0.3]s. Easy, but the cutoff doesn't mean
%                 the same thing in every session: median saccade latency
%                 ranges 68-218ms across the 12 sessions (PIPELINE_REPORT
%                 sec 5.1), so in fast sessions the animal is often
%                 already moving well before 0.1s, and in slow sessions
%                 it's still planning well past 0.1s.
%   'saclocked' -- anchored to each TRIAL's own detected saccade onset
%                 (from analyzeEyeData.m/eyeSummary), which is the
%                 methodologically correct way to do this split:
%                 planning=[-0.15,-0.02]s and execution=[-0.02,+0.15]s
%                 relative to onset, not go-cue. Trials with no detected
%                 saccade (~11% overall, more in some sessions -- see
%                 PIPELINE_REPORT sec 4.4/5.1) are excluded from this mode.
%
% For each mode and each channel: responsiveness (vs. the same
% [-0.9,-0.6]s pre-stimulus baseline used elsewhere) in the planning and
% execution windows separately (Wilcoxon signed-rank + BH-FDR + minimum
% effect size, matching screenResponsiveChannels.m), and direction tuning
% within each window separately for channels responsive in that window
% (vector-sum preferred direction/resultant length + permutation test,
% matching computeDirectionTuning.m).
%
% Does NOT modify or overwrite anything from the main pipeline
% (screening_results.mat, tuning_results.mat, *_overview.png etc020) --
% saves to planning_execution_results.mat and *_planexec_<mode>.png.
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

    % --- mode 1: fixed split relative to go-cue (same window for all trials)
    planIdxGC = S.tb >= opts.goCueWindows.planning(1) & S.tb < opts.goCueWindows.planning(2);
    execIdxGC = S.tb >= opts.goCueWindows.execution(1) & S.tb < opts.goCueWindows.execution(2);
    respPlanGC = squeeze(mean(S.MUA(:, isCorrect, planIdxGC), 3)) - baseTrials;
    respExecGC = squeeze(mean(S.MUA(:, isCorrect, execIdxGC), 3)) - baseTrials;

    r = oneMode(S, respPlanGC, respExecGC, dirDeg(isCorrect), dirsDeg, opts, 'gocue');
    r.Day = S.Day; r.RunN = S.RunN; r.mode = 'gocue';
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
        dirDegSL = dirDeg(useTrials);

        r2 = oneMode(S, respPlanSL, respExecSL, dirDegSL, dirsDeg, opts, 'saclocked');
        r2.Day = S.Day; r2.RunN = S.RunN; r2.mode = 'saclocked';
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
function r = oneMode(S, respPlan, respExec, dirDeg, dirsDeg, opts, modeName) %#ok<INUSD>
nChan = size(respPlan,1);
baseSD = ones(nChan,1); % responses already baseline-subtracted; effect size uses SD of the planning-window diffs themselves as scale
baseSD = std(respPlan, 0, 2, 'omitnan');
baseSD(baseSD==0 | isnan(baseSD)) = eps;

pPlan = nan(nChan,1); effPlan = nan(nChan,1);
pExec = nan(nChan,1); effExec = nan(nChan,1);
for c = 1:nChan
    dP = respPlan(c,:); dP = dP(~isnan(dP));
    dE = respExec(c,:); dE = dE(~isnan(dE));
    if numel(dP) >= 5 && ~all(dP==dP(1))
        pPlan(c) = signrank(dP);
        effPlan(c) = mean(dP) / baseSD(c);
    end
    if numel(dE) >= 5 && ~all(dE==dE(1))
        pExec(c) = signrank(dE);
        effExec(c) = mean(dE) / baseSD(c);
    end
end
fdrPlan = bhFDR(pPlan);
fdrExec = bhFDR(pExec);
respPlanFlag = fdrPlan < opts.alpha & abs(effPlan) > opts.minEffect;
respExecFlag = fdrExec < opts.alpha & abs(effExec) > opts.minEffect;

emptyTuning = struct('prefDir',nan,'resultantLen',nan,'pPerm',nan,'pKW',nan,'tuned',false,'meanByDir',nan(1,8));
tuningPlan = repmat(emptyTuning, 1, nChan);
tuningExec = repmat(emptyTuning, 1, nChan);
for c = 1:nChan
    if respPlanFlag(c)
        tuningPlan(c) = fitTuning(respPlan(c,:), dirDeg, dirsDeg, opts.nPerm);
    end
    if respExecFlag(c)
        tuningExec(c) = fitTuning(respExec(c,:), dirDeg, dirsDeg, opts.nPerm);
    end
end

r = struct();
r.effPlan = effPlan; r.effExec = effExec;
r.respPlanFlag = respPlanFlag; r.respExecFlag = respExecFlag;
r.tuningPlan = tuningPlan; r.tuningExec = tuningExec;
end

% ------------------------------------------------------------------
function plotPlanExec(S, r, opts, modeName) %#ok<INUSL>
nChan = numel(r.effPlan);
chIdx = (1:nChan)';

f = figure('Visible','off','Position',[0 0 1500 600], 'Color', 'w');
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: effect map, planning vs execution, per channel (depth order, deep=bottom)
ax1 = nexttile(tl);
M = [r.effPlan, r.effExec];
M(~[r.respPlanFlag, r.respExecFlag]) = nan;
imagesc(ax1, 1:2, chIdx, M, 'AlphaData', double(~isnan(M)));
set(ax1, 'CLim', [-2 2], 'YDir', 'normal', 'Color', [0.85 0.85 0.85]);
colormap(ax1, redblueLocal(256));
set(ax1, 'XTick', [1 2], 'XTickLabel', {'planning','execution'});
colorbar(ax1);
ylabel(ax1, 'Channel # (1 = deepest)');
title(ax1, sprintf('%d/%d planning, %d/%d execution', sum(r.respPlanFlag), nChan, sum(r.respExecFlag), nChan));

% Panel 2: tuning strength, planning vs execution (colored by category --
% see classifyTuning.m: red=directional, orange=complex/non-unimodal,
% gray=untuned; category shown is whichever is "most tuned" of the two
% windows for that channel, just for this scatter's color)
ax2 = nexttile(tl);
Rplan = [r.tuningPlan.resultantLen]; Rexec = [r.tuningExec.resultantLen];
[catPlan, colPlan] = classifyTuning([r.tuningPlan.pPerm], [r.tuningPlan.pKW]);
[catExec, colExec] = classifyTuning([r.tuningExec.pPerm], [r.tuningExec.pKW]);
isDirPlan = strcmp(catPlan, 'directional');
isDirExec = strcmp(catExec, 'directional');
useCol = colExec; useCol(isDirPlan & ~isDirExec, :) = colPlan(isDirPlan & ~isDirExec, :);
hold(ax2, 'on');
scatter(ax2, Rplan, Rexec, 14, useCol, 'filled');
xlim(ax2, [0 1]); ylim(ax2, [0 1]);
plot(ax2, [0 1], [0 1], 'k:');
set(ax2, 'Color', 'w');
xlabel(ax2, 'Tuning strength, planning'); ylabel(ax2, 'Tuning strength, execution');
title(ax2, 'Tuning strength: planning vs execution');

% Panel 3: preferred direction, planning vs execution (only channels tuned in BOTH)
ax3 = nexttile(tl); hold(ax3, 'on');
tunedBoth = [r.tuningPlan.tuned] & [r.tuningExec.tuned];
if any(tunedBoth)
    PDplan = [r.tuningPlan.prefDir]; PDexec = [r.tuningExec.prefDir];
    scatter(ax3, PDplan(tunedBoth), PDexec(tunedBoth), 20, [0.1 0.4 0.85], 'filled');
end
xlim(ax3, [0 360]); ylim(ax3, [0 360]);
plot(ax3, [0 360], [0 360], 'k:');
set(ax3, 'Color', 'w', 'XTick', 0:90:360, 'YTick', 0:90:360);
xlabel(ax3, 'Preferred direction, planning (deg)'); ylabel(ax3, 'Preferred direction, execution (deg)');
title(ax3, sprintf('PD planning vs execution (n=%d tuned in both)', sum(tunedBoth)));

title(tl, sprintf('%s run-%03d: planning vs execution (%s-anchored)', S.Day, S.RunN, modeName));
forceLightTheme(f);
saveas(f, fullfile(opts.figDir, sprintf('%s_run-%03d_planexec_%s.png', S.Day, S.RunN, modeName)));
close(f);
end

% ------------------------------------------------------------------
function plotPlanExecCurves(S, r, dirsDeg, opts, modeName) %#ok<INUSL>
% Per-channel small-multiple tuning curves, planning vs execution shown
% side by side for the SAME channel (depth order), so you can see
% directly whether a channel's tuning shape/direction changes between
% the planning and execution epoch, not just whether its summary
% strength/preferred-direction numbers match (see plotPlanExec panels 2-3).
nChan = numel(r.effPlan);
respAny = r.respPlanFlag | r.respExecFlag;
chans = find(respAny);
nShow = numel(chans);
if nShow == 0, return; end

[catPlan, colPlan] = classifyTuning([r.tuningPlan.pPerm], [r.tuningPlan.pKW]);
[catExec, colExec] = classifyTuning([r.tuningExec.pPerm], [r.tuningExec.pKW]);

f = figure('Visible','off','Position',[0 0 1600 max(600,80*ceil(nShow/6))], 'Color', 'w');
nCol = 6; nRow = ceil(nShow/nCol)*2; % 2 rows per channel (planning on top, execution below)
tl = tiledlayout(f, nRow, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nShow
    c = chans(i);
    rowBlock = floor((i-1)/nCol);
    colPos = mod(i-1,nCol);
    tileTop = rowBlock*2*nCol + colPos + 1;
    tileBot = tileTop + nCol;

    axP = nexttile(tl, tileTop);
    vP = r.tuningPlan(c).meanByDir;
    plot(axP, [dirsDeg dirsDeg(1)+360], [vP vP(1)], '-o', 'Color', colPlan(i,:), ...
        'MarkerFaceColor', colPlan(i,:), 'LineWidth', 1.2, 'MarkerSize', 3);
    set(axP, 'XTick', [], 'Color', 'w');
    title(axP, sprintf('ch%d: plan (%s)', c, catPlan{i}), 'FontSize', 6, 'FontWeight', 'normal');

    axE = nexttile(tl, tileBot);
    vE = r.tuningExec(c).meanByDir;
    plot(axE, [dirsDeg dirsDeg(1)+360], [vE vE(1)], '-o', 'Color', colExec(i,:), ...
        'MarkerFaceColor', colExec(i,:), 'LineWidth', 1.2, 'MarkerSize', 3);
    set(axE, 'XTick', [], 'Color', 'w');
    title(axE, sprintf('ch%d: exec (%s)', c, catExec{i}), 'FontSize', 6, 'FontWeight', 'normal');
end
title(tl, sprintf(['%s run-%03d: planning (top) vs execution (bottom) tuning shape, %s-anchored ' ...
    '-- red=directional, orange=complex, gray=untuned'], S.Day, S.RunN, modeName));
forceLightTheme(f);
saveas(f, fullfile(opts.figDir, sprintf('%s_run-%03d_planexec_%s_curves.png', S.Day, S.RunN, modeName)));
close(f);
end

% ------------------------------------------------------------------
function cmap = redblueLocal(n)
half = floor(n/2);
top = [linspace(0.1,1,half)', linspace(0.1,1,half)', ones(half,1)];
bot = [ones(n-half,1), linspace(1,0.1,n-half)', linspace(1,0.1,n-half)'];
cmap = [top; bot];
end
