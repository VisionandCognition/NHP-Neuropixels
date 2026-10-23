function summary = summarizePlanningExecution(mode, excludeDays)
% SUMMARIZEPLANNINGEXECUTION Quantifies the planning-vs-execution split
% (analyzePlanningExecution.m) across ALL 12 sessions, rather than the
% single-session "first look" in PIPELINE_REPORT.md sec 4.8. For each
% session, reports:
%   - fraction of channels 'directional' (classifyTuning.m) in execution
%     but NOT in planning ("emerges during execution")
%   - fraction 'directional' in planning but NOT execution ("emerges
%     during planning" -- the reverse, checked for completeness)
%   - for channels 'directional' in BOTH: the circular consistency of
%     preferred direction between planning and execution (mean absolute
%     angular difference; 0 = perfectly consistent, 90 = orthogonal,
%     random pairs average ~90)
%
% summary = summarizePlanningExecution(mode)
%   mode : 'gocue' or 'saclocked' (default 'saclocked', the
%          methodologically correct anchoring -- see PIPELINE_REPORT sec 4.8)
%   excludeDays : cellstr of sessions left out of the POOLED totals (still
%          listed per-session). Default: the four sessions with near-zero
%          directional channels in either window (20260224, 20260326,
%          20260910, 20260916 -- <=3 directional channels per window,
%          same convention as poolPlanExecByArea.m).

if nargin < 1 || isempty(mode), mode = 'saclocked'; end
if nargin < 2, excludeDays = {'20260224','20260326','20260910','20260916'}; end

load(fullfile(saccadeDataDir(), 'planning_execution_results.mat'), 'results');
allModes = cellfun(@(x) x.mode, {results.data}, 'UniformOutput', false);
rows = find(strcmp(allModes, mode));

summary = struct([]);
fprintf('Day\t\tnResp\tnPlanDir\tnExecDir\tnBothDir\texecOnly\tplanOnly\tmeanAngDiff(bothDir)\n');
for i = 1:numel(rows)
    r = results(rows(i)).data;
    nChan = numel(r.tuningPlan);

    % analyzePlanningExecution.m's tuning structs don't store a category
    % field directly (only pPerm/pKW) -- classifyTuning.m derives it here,
    % same as the plotting functions do.
    catP = classifyTuning([r.tuningPlan.pPerm], [r.tuningPlan.pKW]);
    catE = classifyTuning([r.tuningExec.pPerm], [r.tuningExec.pKW]);
    % 'directional'+'directional_complex' == old single 'directional'
    % category (classifyTuning.m 2026-07-14 split) -- ismember preserves
    % the original "has some real preferred direction" membership.
    isDirPlan = ismember(catP, {'directional','directional_complex'});
    isDirExec = ismember(catE, {'directional','directional_complex'});

    both = isDirPlan & isDirExec;
    execOnly = isDirExec & ~isDirPlan;
    planOnly = isDirPlan & ~isDirExec;

    angDiff = nan;
    if any(both)
        pdP = [r.tuningPlan(both).prefDir];
        pdE = [r.tuningExec(both).prefDir];
        d = abs(mod(pdP - pdE + 180, 360) - 180); % circular abs difference, 0-180
        angDiff = mean(d);
    end

    summary(i).Day = r.Day;
    summary(i).RunN = r.RunN;
    summary(i).nChan = nChan;
    summary(i).nPlanDirectional = sum(isDirPlan);
    summary(i).nExecDirectional = sum(isDirExec);
    summary(i).nBothDirectional = sum(both);
    summary(i).nExecOnly = sum(execOnly);
    summary(i).nPlanOnly = sum(planOnly);
    summary(i).meanAngDiffBoth = angDiff;

    fprintf('%s\t%d\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%.1f\n', r.Day, nChan, sum(isDirPlan), sum(isDirExec), ...
        sum(both), sum(execOnly), sum(planOnly), angDiff);
end

inPool = ~ismember({summary.Day}, excludeDays);
totalExecOnly = sum([summary(inPool).nExecOnly]);
totalPlanOnly = sum([summary(inPool).nPlanOnly]);
totalBoth = sum([summary(inPool).nBothDirectional]);
totalExecDir = sum([summary(inPool).nExecDirectional]);
totalPlanDir = sum([summary(inPool).nPlanDirectional]);
fprintf('\n=== Pooled across %d sessions (%s mode; excluded: %s) ===\n', sum(inPool), mode, strjoin(excludeDays, ' '));
fprintf('directional in execution: %d, of which also directional in planning: %d (%.0f%%), execution-only: %d (%.0f%%)\n', ...
    totalExecDir, totalBoth, 100*totalBoth/max(totalExecDir,1), totalExecOnly, 100*totalExecOnly/max(totalExecDir,1));
fprintf('directional in planning: %d, planning-only (not execution): %d (%.0f%%)\n', ...
    totalPlanDir, totalPlanOnly, 100*totalPlanOnly/max(totalPlanDir,1));
validAng = [summary(inPool).meanAngDiffBoth]; validAng = validAng(~isnan(validAng));
fprintf('mean angular difference (planning vs execution preferred direction, channels tuned in both), pooled: %.1f deg (n=%d sessions with >=1 such channel)\n', ...
    mean(validAng), numel(validAng));
end
