function areaResults = decodeByArea(sessions, screening, opts)
% DECODEBYAREA Does the SC-stands-out finding from poolTuningByArea.m
% (§6.7 -- 42.6% directional vs. 5-20% for most other areas) also show up
% in DECODING, not just tuning fraction? For each area token found by
% extractAreaTokens.m (same liberal multi-label pooling as
% poolTuningByArea.m), restricts each session's OWN decoding population
% (decodeDirection.m's method: nearest-centroid, leave-one-trial-out,
% all 3 epochs) to just that area's RESPONSIVE channels in that session,
% then pools the resulting accuracy across every session that has enough
% area-X channels to attempt it.
%
% Trial-level data can't be naively concatenated across sessions the way
% per-channel tuning summaries can (§7 item 8 in PIPELINE_REPORT.md), so
% this is necessarily "decode per session, then average the accuracy
% numbers across sessions" rather than one pooled decode -- each session
% contributes one accuracy value per area per epoch, not raw channels.
%
% SIZE-MATCHED for fairness: raw area-matching channel counts vary hugely
% (SC has ~100+ channels/session in some sessions; several other areas
% have as few as 20 total across all sessions), which would confound a
% naive "SC decodes better" comparison with "SC just has a bigger
% population." Every session's area-X decode is capped at
% opts.capSize channels (randomly subsampled, opts.nReps repeats,
% averaged) if it has more, or uses all available (no subsampling needed)
% if it has fewer.
%
% areaResults = decodeByArea(sessions, screening, opts)
%   sessions, screening : as in decodeDirection.m (sessions must come
%     from loadAllSessions() so S.channelArea is populated)
%   opts.baseWindow, opts.visualWindow, opts.goCueWindows : same defaults
%     as decodeDirection.m, so epochs match exactly
%   opts.excludeTokens        : cellstr, non-area qualifier tokens to
%                               drop (default {'DEEP','HIGH'}, same as
%                               poolTuningByArea.m)
%   opts.minChannelsPerSession : minimum area-X responsive channels a
%                               session must have to attempt a decode
%                               there (default 15)
%   opts.minSessions          : an area must be decodable in at least
%                               this many sessions to be reported (default 2)
%   opts.capSize              : population size cap for fair cross-area
%                               comparison (default 20)
%   opts.nReps                : random-subsample repeats when a session
%                               has more than capSize channels (default 20)
%   opts.figDir               : default <saccadeDataDir()>/figures

if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'baseWindow', [-0.9 -0.6]);
opts = setDefault(opts, 'visualWindow', [-0.5 -0.35]);
opts = setDefault(opts, 'goCueWindows', struct('planning', [0 0.1], 'execution', [0.1 0.3]));
opts = setDefault(opts, 'excludeTokens', {'DEEP','HIGH'});
opts = setDefault(opts, 'minChannelsPerSession', 15);
opts = setDefault(opts, 'minSessions', 2);
opts = setDefault(opts, 'capSize', 20);
opts = setDefault(opts, 'nReps', 20);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));
if ~isfolder(opts.figDir), mkdir(opts.figDir); end

epochNames = {'visual','planning','execution'};

% --- per-session: epoch response matrices + area-token set per channel
sessData = struct([]);
for s = 1:numel(sessions)
    S = sessions(s);
    isCorrect = strcmp(S.trialTable.outcome, 'correct') & ~isnan(S.trialTable.direction_deg);
    dirGroup = round(S.trialTable.direction_deg(isCorrect)'/45) + 1;

    baseIdx = S.tb >= opts.baseWindow(1) & S.tb < opts.baseWindow(2);
    baseTrials = squeeze(mean(S.MUA(:, isCorrect, baseIdx), 3));
    visIdx = S.tb >= opts.visualWindow(1) & S.tb < opts.visualWindow(2);
    planIdx = S.tb >= opts.goCueWindows.planning(1) & S.tb < opts.goCueWindows.planning(2);
    execIdx = S.tb >= opts.goCueWindows.execution(1) & S.tb < opts.goCueWindows.execution(2);
    epochResp = { ...
        squeeze(mean(S.MUA(:, isCorrect, visIdx), 3)) - baseTrials, ...
        squeeze(mean(S.MUA(:, isCorrect, planIdx), 3)) - baseTrials, ...
        squeeze(mean(S.MUA(:, isCorrect, execIdx), 3)) - baseTrials};

    respChans = find(screening(s).responsive);
    nChan = numel(screening(s).responsive);
    chanTokens = cell(nChan,1);
    for c = respChans'
        combined = strjoin({S.channelArea.AssignedArea{c}, S.channelArea.L6{c}, ...
            S.channelArea.L5{c}, S.channelArea.L4{c}, S.channelArea.L3{c}}, ' ');
        chanTokens{c} = setdiff(extractAreaTokens(combined), opts.excludeTokens);
    end

    sessData(s).Day = S.Day; %#ok<AGROW>
    sessData(s).RunN = S.RunN;
    sessData(s).epochResp = epochResp;
    sessData(s).dirGroup = dirGroup;
    sessData(s).respChans = respChans;
    sessData(s).chanTokens = chanTokens;
end

allTokens = {};
for s = 1:numel(sessData)
    allTokens = [allTokens, sessData(s).chanTokens{sessData(s).respChans}]; %#ok<AGROW>
end
allTokens = unique(allTokens);

fprintf('%-10s %7s %8s %8s %8s\n', 'Area', 'NSess', 'visual%', 'planning%', 'execution%');
areaResults = struct([]);
for a = 1:numel(allTokens)
    area = allTokens{a};
    accBySession = nan(numel(sessData), 3);
    nChanUsed = nan(numel(sessData), 1);
    for s = 1:numel(sessData)
        sd = sessData(s);
        matchChans = sd.respChans(cellfun(@(t) any(strcmp(t, area)), sd.chanTokens(sd.respChans)));
        if numel(matchChans) < opts.minChannelsPerSession, continue; end
        nChanUsed(s) = min(numel(matchChans), opts.capSize);
        for e = 1:3
            if numel(matchChans) <= opts.capSize
                accBySession(s,e) = looDecodeAcc(sd.epochResp{e}(matchChans,:), sd.dirGroup);
            else
                accs = nan(opts.nReps,1);
                for rep = 1:opts.nReps
                    sel = matchChans(randperm(numel(matchChans), opts.capSize));
                    accs(rep) = looDecodeAcc(sd.epochResp{e}(sel,:), sd.dirGroup);
                end
                accBySession(s,e) = mean(accs);
            end
        end
    end
    validSess = find(~isnan(accBySession(:,1)));
    if numel(validSess) < opts.minSessions, continue; end

    r = numel(areaResults) + 1;
    areaResults(r).Area = area; %#ok<AGROW>
    areaResults(r).NSessions = numel(validSess);
    areaResults(r).MeanNChanUsed = mean(nChanUsed(validSess));
    for e = 1:3
        areaResults(r).(sprintf('%s_mean', epochNames{e})) = mean(accBySession(validSess,e));
        areaResults(r).(sprintf('%s_sem', epochNames{e})) = std(accBySession(validSess,e)) / sqrt(numel(validSess));
    end

    fprintf('%-10s %7d %7.1f%% %7.1f%% %7.1f%%\n', area, numel(validSess), ...
        100*areaResults(r).visual_mean, 100*areaResults(r).planning_mean, 100*areaResults(r).execution_mean);
end

[~, order] = sort([areaResults.NSessions], 'descend');
areaResults = areaResults(order);

outFile = fullfile(saccadeDataDir(), 'decode_by_area.mat');
save(outFile, 'areaResults');
fprintf('saved %s\n', outFile);

plotDecodeByArea(areaResults, opts);
end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function acc = looDecodeAcc(resp, dirGroup)
% Identical method to decodeDirection.m's looDecodeAcc (duplicated
% locally, matching this codebase's existing convention).
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
function plotDecodeByArea(areaResults, opts)
if isempty(areaResults), return; end
n = numel(areaResults);
colors = {[0.2 0.4 0.8], [0.85 0.55 0.1], [0.75 0.1 0.1]};
epochNames = {'visual','planning','execution'};

f = figure('Visible','off','Position',[0 0 max(700, 90*n) 500], 'Color', 'w');
ax = axes(f); hold(ax, 'on');
x = 1:n;
hBar = gobjects(1,3);
for e = 1:3
    m = [areaResults.(sprintf('%s_mean', epochNames{e}))] * 100;
    se = [areaResults.(sprintf('%s_sem', epochNames{e}))] * 100;
    xOff = x + (e-2)*0.25;
    hBar(e) = bar(ax, xOff, m, 0.25, 'FaceColor', colors{e});
    errorbar(ax, xOff, m, se, 'k.', 'LineWidth', 1, 'CapSize', 3);
end
yline(ax, 12.5, 'k:', 'LineWidth', 1);
xLabels = arrayfun(@(p) sprintf('%s (n=%d sess, ~%.0f chan)', p.Area, p.NSessions, p.MeanNChanUsed), areaResults, 'UniformOutput', false);
set(ax, 'XTick', x, 'XTickLabel', xLabels, 'XTickLabelRotation', 45, 'Color', 'w', 'TickLabelInterpreter', 'none');
ylabel(ax, 'Decode accuracy (%)');
% NB: legend(ax, epochNames) alone grabs whatever 3 graphics objects it
% finds first in the axes' children (a mix of bar AND errorbar handles,
% since each loop iteration creates both) -- explicit handles needed so
% the legend entries actually match the bar colors.
legend(ax, hBar, epochNames, 'Location', 'northeast', 'TextColor', 'k', 'Color', 'w');
title(ax, sprintf('Direction decoding by area (population size-matched at %d channels, chance=12.5%%)', opts.capSize), 'FontSize', 9, 'Interpreter', 'none');
forceLightTheme(f);
if ~isfolder(fullfile(opts.figDir,'pooled')), mkdir(fullfile(opts.figDir,'pooled')); end
saveas(f, fullfile(opts.figDir, 'pooled', 'decode_by_area.png'));
close(f);
fprintf('saved %s\n', fullfile(opts.figDir, 'pooled', 'decode_by_area.png'));
end
