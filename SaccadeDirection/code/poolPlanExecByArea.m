function areaResults = poolPlanExecByArea(sessions, opts)
% POOLPLANEXECBYAREA Redoes §6.6's planning-vs-execution recruitment
% split ("53% of execution-directional channels were planning-untuned")
% BY AREA instead of pooled across all channels: does execution-only
% recruitment concentrate in SC/motor-adjacent areas, or is it uniform
% regardless of area? Same liberal multi-label area pooling as
% poolTuningByArea.m/decodeByArea.m (extractAreaTokens.m) -- a channel
% counts toward every area token found within its ONE labelField's label.
%
% Earlier versions of this function pulled tokens from a concatenation of
% AssignedArea+L6+L5+L4+L3, which conflated genuine within-level candidate
% uncertainty ("or X or Y") with cross-level hierarchy (e.g. "SC" at the
% fine level and "colliculi (Co)" at a coarser level are the SAME
% channels, not two candidate identities). Fixed: operates on exactly one
% labelField at a time; call it once per level you want (e.g. 'L6' and
% 'L3') to compare grains.
%
% Uses analyzePlanningExecution.m's saccade-onset-locked ('saclocked')
% results by default -- the methodologically preferred alignment (§4.8) --
% and excludes 20260224/20260326/20260910/20260916 (near-zero directional channels in
% either window, same convention as summarizePlanningExecution.m).
%
% For each area token with enough execution-directional channel-instances
% (opts.minExecDirectional, to avoid noisy small-N percentages):
%   - NExecDirectional, NPlanDirectional (pooled across qualifying sessions)
%   - PctExecOnly  : % of execution-directional channels NOT also
%                    directional during planning (the "recruited during
%                    movement" fraction, §6.6's headline number, here per area)
%   - PctPlanOnly  : % of planning-directional channels NOT directional
%                    during execution
%
% areaResults = poolPlanExecByArea(sessions, opts)
%   sessions : struct array from loadAllSessions() (needs S.channelArea)
%   opts.labelField           : which S.channelArea field to pool by --
%                               'AssignedArea','L6','L5','L4','L3'
%                               (default 'L6', finest hierarchy level)
%   opts.mode                : 'saclocked' or 'gocue' (default 'saclocked')
%   opts.excludeTokens        : cellstr, non-area qualifiers to drop
%                               (default {'DEEP','HIGH'}, same as poolTuningByArea.m)
%   opts.excludeDays          : cellstr, sessions to skip (default
%                               {'20260224','20260326','20260910','20260916'})
%   opts.minExecDirectional   : minimum execution-directional
%                               channel-instances before reporting an
%                               area (default 8)
%   opts.figDir               : default <saccadeDataDir()>/figures

if nargin < 2, opts = struct(); end
opts = setDefault(opts, 'labelField', 'L6');
opts = setDefault(opts, 'mode', 'saclocked');
opts = setDefault(opts, 'excludeTokens', {'DEEP','HIGH'});
opts = setDefault(opts, 'excludeDays', {'20260224','20260326','20260910','20260916'});
opts = setDefault(opts, 'minExecDirectional', 8);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));
if ~isfolder(opts.figDir), mkdir(opts.figDir); end

load(fullfile(saccadeDataDir(), 'planning_execution_results.mat'), 'results');
allModes = cellfun(@(x) x.mode, {results.data}, 'UniformOutput', false);
rows = find(strcmp(allModes, opts.mode));

rowData = struct('Tokens', {}, 'IsDirPlan', {}, 'IsDirExec', {});
for ri = 1:numel(rows)
    r = results(rows(ri)).data;
    if ismember(r.Day, opts.excludeDays), continue; end
    s = find(strcmp({sessions.Day}, r.Day) & [sessions.RunN] == r.RunN);
    if isempty(s), continue; end
    S = sessions(s);
    nChan = numel(r.tuningPlan);
    catPlan = classifyTuning([r.tuningPlan.pPerm], [r.tuningPlan.pKW]);
    catExec = classifyTuning([r.tuningExec.pPerm], [r.tuningExec.pKW]);
    % 'directional'+'directional_complex' == old single 'directional'
    % category (classifyTuning.m 2026-07-14 split) -- ismember preserves
    % the original "has some real preferred direction" membership.
    isDirPlan = ismember(catPlan, {'directional','directional_complex'});
    isDirExec = ismember(catExec, {'directional','directional_complex'});
    for c = 1:nChan
        toks = setdiff(extractAreaTokens(S.channelArea.(opts.labelField){c}), opts.excludeTokens);
        if isempty(toks), continue; end
        i = numel(rowData) + 1;
        rowData(i).Tokens = toks; %#ok<AGROW>
        rowData(i).IsDirPlan = isDirPlan(c);
        rowData(i).IsDirExec = isDirExec(c);
    end
end

allTokens = unique([rowData.Tokens]);
areaResults = struct([]);
fprintf('%-10s %8s %8s %10s %10s\n', 'Area', 'NExecDir', 'NPlanDir', '%ExecOnly', '%PlanOnly');
for a = 1:numel(allTokens)
    area = allTokens{a};
    inPool = cellfun(@(t) any(strcmp(t, area)), {rowData.Tokens});
    idx = find(inPool);
    isDirPlan = [rowData(idx).IsDirPlan];
    isDirExec = [rowData(idx).IsDirExec];
    nExecDir = sum(isDirExec);
    if nExecDir < opts.minExecDirectional, continue; end
    nPlanDir = sum(isDirPlan);
    nExecOnly = sum(isDirExec & ~isDirPlan);
    nBoth = sum(isDirExec & isDirPlan);
    nPlanOnly = sum(isDirPlan & ~isDirExec);

    r = numel(areaResults) + 1;
    areaResults(r).Area = area; %#ok<AGROW>
    areaResults(r).NChannels = numel(idx);
    areaResults(r).NExecDirectional = nExecDir;
    areaResults(r).NPlanDirectional = nPlanDir;
    areaResults(r).NBothDirectional = nBoth;
    areaResults(r).PctExecOnly = 100 * nExecOnly / nExecDir;
    areaResults(r).PctPlanOnly = 100 * nPlanOnly / max(nPlanDir, 1);

    fprintf('%-10s %8d %8d %9.1f%% %9.1f%%\n', area, nExecDir, nPlanDir, areaResults(r).PctExecOnly, areaResults(r).PctPlanOnly);
end

[~, order] = sort([areaResults.PctExecOnly], 'descend'); % sorted to make the range across areas visible at a glance, not buried in the middle
areaResults = areaResults(order);

outFile = fullfile(saccadeDataDir(), sprintf('planexec_pooled_by_area_%s.mat', opts.labelField));
save(outFile, 'areaResults');
fprintf('saved %s\n', outFile);

plotPlanExecByArea(areaResults, opts);
end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function plotPlanExecByArea(areaResults, opts)
if isempty(areaResults), return; end
n = numel(areaResults);
f = figure('Visible','off','Position',[0 0 max(700, 70*n) 500], 'Color', 'w');
ax = axes(f); hold(ax, 'on');
x = 1:n;
hBar = gobjects(1,2);
hBar(1) = bar(ax, x-0.15, [areaResults.PctExecOnly], 0.3, 'FaceColor', [0.75 0.1 0.1]);
hBar(2) = bar(ax, x+0.15, [areaResults.PctPlanOnly], 0.3, 'FaceColor', [0.2 0.4 0.8]);
xLabels = arrayfun(@(p) sprintf('%s (execDir=%d, planDir=%d)', p.Area, p.NExecDirectional, p.NPlanDirectional), areaResults, 'UniformOutput', false);
set(ax, 'XTick', x, 'XTickLabel', xLabels, 'XTickLabelRotation', 45, 'Color', 'w', 'TickLabelInterpreter', 'none');
ylabel(ax, '%');
legend(ax, hBar, {'% execution-directional NOT tuned in planning ("execution-only")', ...
    '% planning-directional NOT tuned in execution ("planning-only")'}, 'Location', 'northeast', 'TextColor', 'k', 'Color', 'w');
title(ax, sprintf('Planning-vs-execution recruitment by area, pooled by %s (saccade-onset-locked, min %d exec-directional channels/area)', opts.labelField, opts.minExecDirectional), 'FontSize', 9, 'Interpreter', 'none');
forceLightTheme(f);
if ~isfolder(fullfile(opts.figDir,'pooled')), mkdir(fullfile(opts.figDir,'pooled')); end
outFile = fullfile(opts.figDir, 'pooled', sprintf('planexec_by_area_%s.png', opts.labelField));
saveas(f, outFile);
close(f);
fprintf('saved %s\n', outFile);
end
