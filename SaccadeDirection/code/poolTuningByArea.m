function poolResults = poolTuningByArea(sessions, screening, opts)
% POOLTUNINGBYAREA Pools responsiveness/tuning stats across ALL 12
% sessions BY ANATOMICAL AREA (using S.channelArea, §4.12/§4.13) instead
% of by session -- the thing §6.3 explicitly avoided doing before
% ("sessions are independent penetrations... not seated to the same
% depth"), now possible because channels carry an actual area label
% rather than just an opaque depth index.
%
% LIBERAL POOLING: many hand-assigned labels genuinely span multiple
% candidate areas (e.g. "Posterior Commisure (or Habenula or CM?)"), so a
% channel counts toward EVERY area token found in ANY of its label
% fields (AssignedArea, L6, L5, L4, L3), not just its single primary
% label -- see extractAreaTokens.m. This means a channel can appear in
% more than one area's pool; that is intentional, not a bug.
%
% For each area token found anywhere in the dataset (except
% opts.excludeTokens -- a couple of tokens that turned out to be
% position qualifiers, not areas, e.g. "(deep)"/"(high up)" describing
% which end of a compound cluster label, not a candidate structure):
%   - NChannels, NSessions represented
%   - %responsive, %directional, %complex (of screenResponsiveChannels.m /
%     classifyTuning.m's categories)
%   - pooled dominant preferred direction + concentration (resultant-
%     length-weighted circular mean across all directional channels in
%     the pool, regardless of session)
%   - mean tuning strength (resultant length)
%
% poolResults = poolTuningByArea(sessions, screening, opts)
%   sessions, screening : as elsewhere (sessions must come from
%     loadAllSessions.m so S.channelArea is populated)
%   opts.excludeTokens : cellstr of tokens to drop (default {'DEEP','HIGH'})
%   opts.minChannels   : only report/plot areas with at least this many
%                        pooled channels (default 20)
%   opts.minDirForRose : only include an area in the preferred-direction
%                        rose-plot grid if it has at least this many
%                        directional channels (default 15)
%   opts.figDir        : default <saccadeDataDir()>/figures

if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'excludeTokens', {'DEEP','HIGH'});
opts = setDefault(opts, 'minChannels', 20);
opts = setDefault(opts, 'minDirForRose', 15);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));
if ~isfolder(opts.figDir), mkdir(opts.figDir); end

load(fullfile(saccadeDataDir(), 'tuning_results.mat'), 'tuning');

% --- build one row per (session, channel): area tokens + tuning outcome.
% nameCounts: token -> Map(full-name-phrase -> occurrence count), so the
% most frequently-seen phrasing for each token can be reported as its
% human-readable full name (labels aren't perfectly consistent across
% sessions, e.g. capitalization, so this picks the most common one seen
% rather than just the first).
rows = struct('Tokens', {}, 'Category', {}, 'PrefDir', {}, 'ResultantLen', {});
nameCounts = containers.Map();
for s = 1:numel(sessions)
    S = sessions(s);
    sc = screening(s);
    nChan = numel(sc.responsive);
    tRows = find(strcmp({tuning.Day}, S.Day) & [tuning.RunN] == S.RunN);
    tByChan = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for r = tRows, tByChan(tuning(r).channel) = tuning(r); end %#ok<NASGU>

    for c = 1:nChan
        combined = strjoin({S.channelArea.AssignedArea{c}, S.channelArea.L6{c}, ...
            S.channelArea.L5{c}, S.channelArea.L4{c}, S.channelArea.L3{c}}, ' ');
        [rawToks, rawNames] = extractAreaTokens(combined);
        keep = ~ismember(rawToks, opts.excludeTokens);
        toks = rawToks(keep);
        if isempty(toks), continue; end % no area info for this channel at all -- excluded from every pool

        for k = 1:numel(toks)
            nm = rawNames{strcmp(rawToks, toks{k})};
            if isempty(nm), continue; end
            if isKey(nameCounts, toks{k})
                m = nameCounts(toks{k});
            else
                m = containers.Map('KeyType', 'char', 'ValueType', 'double');
            end
            if isKey(m, nm), m(nm) = m(nm) + 1; else, m(nm) = 1; end %#ok<NASGU>
            nameCounts(toks{k}) = m;
        end

        i = numel(rows) + 1;
        rows(i).Tokens = toks;
        if sc.responsive(c) && isKey(tByChan, c)
            tr = tByChan(c);
            rows(i).Category = tr.category;
            rows(i).PrefDir = tr.prefDir;
            rows(i).ResultantLen = tr.resultantLen;
        else
            rows(i).Category = 'not_responsive';
            rows(i).PrefDir = nan;
            rows(i).ResultantLen = 0;
        end
    end
end

allTokens = unique([rows.Tokens]);
poolResults = struct([]);
for a = 1:numel(allTokens)
    area = allTokens{a};
    inPool = cellfun(@(t) any(strcmp(t, area)), {rows.Tokens});
    idx = find(inPool);
    if isempty(idx), continue; end

    cats = {rows(idx).Category};
    isDir = strcmp(cats, 'directional');
    isCx = strcmp(cats, 'complex');
    isResp = ~strcmp(cats, 'not_responsive');
    lens = [rows(idx).ResultantLen];
    prefDirs = [rows(idx).PrefDir];

    domPrefDir = nan; conc = nan;
    if any(isDir)
        w = lens(isDir); pd = prefDirs(isDir);
        cx = sum(w.*cosd(pd)); cy = sum(w.*sind(pd));
        domPrefDir = mod(atan2d(cy, cx), 360);
        conc = hypot(cx, cy) / max(sum(w), eps);
    end

    r = numel(poolResults) + 1;
    poolResults(r).Area = area; %#ok<AGROW>
    poolResults(r).FullName = mostCommonName(nameCounts, area);
    poolResults(r).NChannels = numel(idx);
    poolResults(r).PctResponsive = 100*mean(isResp);
    poolResults(r).PctDirectional = 100*mean(isDir);
    poolResults(r).PctComplex = 100*mean(isCx);
    poolResults(r).NDirectional = sum(isDir);
    poolResults(r).DominantPrefDir = domPrefDir;
    poolResults(r).Concentration = conc;
    poolResults(r).MeanTuningStrength = mean(lens);
    poolResults(r).DirectionalPrefDirs = prefDirs(isDir); % raw per-channel preferred directions, for the distribution plot -- NOT just the collapsed vector mean
end

[~, order] = sort([poolResults.NChannels], 'descend');
poolResults = poolResults(order);

fprintf('%-10s %-32s %7s %8s %8s %8s %11s %7s %7s\n', 'Area', 'FullName', 'NChan', '%Resp', '%Dir', '%Cx', 'domPrefDir', 'conc', 'nDir');
for a = 1:numel(poolResults)
    p = poolResults(a);
    fprintf('%-10s %-32s %7d %7.1f%% %7.1f%% %7.1f%% %10.0fdeg %7.2f %7d\n', ...
        p.Area, p.FullName, p.NChannels, p.PctResponsive, p.PctDirectional, p.PctComplex, p.DominantPrefDir, p.Concentration, p.NDirectional);
end

outFile = fullfile(saccadeDataDir(), 'tuning_pooled_by_area.mat');
save(outFile, 'poolResults');
fprintf('saved %s\n', outFile);

plotPooledByArea(poolResults, opts);
end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function plotAreaLegend(show, opts)
% Plain text table, Area (abbreviation) -> FullName -> NChannels, sorted
% to match the bar chart. Kept as a SEPARATE figure rather than folded
% into the bar chart's x-tick labels -- ~40 areas' full descriptive names
% would not fit legibly as rotated tick labels.
n = numel(show);
f = figure('Visible','off','Position',[0 0 700 max(300,24*n+60)], 'Color', 'w');
ax = axes(f, 'Position', [0.02 0.02 0.96 0.96]);
axis(ax, 'off');
lines = cell(n+1,1);
lines{1} = sprintf('%-10s %-45s %s', 'Area', 'Full name (most common phrasing)', 'NChan');
for i = 1:n
    lines{i+1} = sprintf('%-10s %-45s %d', show(i).Area, show(i).FullName, show(i).NChannels);
end
text(ax, 0, 1, strjoin(lines, newline), 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontName', 'FixedWidth', 'FontSize', 9, 'Interpreter', 'none', 'Color', 'k');
forceLightTheme(f);
saveas(f, fullfile(opts.figDir, 'pooled', 'pooled_tuning_by_area_legend.png'));
close(f);
fprintf('saved %s\n', fullfile(opts.figDir, 'pooled', 'pooled_tuning_by_area_legend.png'));
end

% ------------------------------------------------------------------
function nm = mostCommonName(nameCounts, token)
% Picks the most frequently-seen full-name phrase for a token across the
% whole dataset (nameCounts: token -> Map(phrase -> count), built while
% scanning every channel's labels). Falls back to the bare token itself
% if no descriptive phrase was ever captured for it (e.g. a token that
% only ever appeared via the simple abbreviation-only regex, with no
% clean preceding phrase found).
nm = token;
if ~isKey(nameCounts, token), return; end
m = nameCounts(token);
if m.Count == 0, return; end
phrases = keys(m); counts = cell2mat(values(m));
[~, best] = max(counts);
nm = phrases{best};
end

% ------------------------------------------------------------------
function plotPooledByArea(poolResults, opts)
show = poolResults([poolResults.NChannels] >= opts.minChannels);
if isempty(show), return; end
[~, order] = sort([show.NChannels], 'descend');
show = show(order);
nArea = numel(show);

f = figure('Visible','off','Position',[0 0 max(700, 60*nArea) 500], 'Color', 'w');
ax = axes(f); hold(ax, 'on');
x = 1:nArea;
bar(ax, x-0.25, [show.PctResponsive], 0.25, 'FaceColor', [0.5 0.5 0.5]);
bar(ax, x, [show.PctDirectional], 0.25, 'FaceColor', [0.75 0.1 0.1]);
bar(ax, x+0.25, [show.PctComplex], 0.25, 'FaceColor', [0.85 0.55 0.1]);
xLabels = arrayfun(@(p) sprintf('%s (n=%d)', p.Area, p.NChannels), show, 'UniformOutput', false);
set(ax, 'XTick', x, 'XTickLabel', xLabels, 'XTickLabelRotation', 45, 'Color', 'w', 'TickLabelInterpreter', 'none');
ylabel(ax, '%');
legend(ax, {'% responsive', '% directional', '% complex'}, 'Location', 'northeast', 'TextColor', 'k', 'Color', 'w');
title(ax, sprintf('Pooled by area across all sessions (liberal pooling, min %d channels/area, chan can count toward multiple areas)', opts.minChannels), 'FontSize', 9);
forceLightTheme(f);
if ~isfolder(fullfile(opts.figDir,'pooled')), mkdir(fullfile(opts.figDir,'pooled')); end
saveas(f, fullfile(opts.figDir, 'pooled', 'pooled_tuning_by_area_summary.png'));
close(f);
fprintf('saved %s\n', fullfile(opts.figDir, 'pooled', 'pooled_tuning_by_area_summary.png'));

plotAreaLegend(show, opts);

% Preferred-direction DISTRIBUTION (polar histogram) per area, only for
% areas with enough directional channels -- NOT just the collapsed
% vector-mean line this used to show. There's no a priori reason a given
% area should have one preferred direction at all: for a retinotopically
% organized structure (e.g. SC), preferred direction should track WHERE
% in the map a given penetration happened to land, not converge to a
% single value across sessions -- and any area's pooled distribution
% here is itself shaped by which part of that area each penetration
% happened to sample (§6.7 caveat), not just intrinsic properties of the
% area. The actual per-channel distribution is shown directly (8 bins,
% matching the task's 8 discrete directions); the resultant vector mean
% (dotted line + concentration, as before) is kept as a reference
% overlay, not the headline result.
% Manual polaraxes positioning (not tiledlayout/nexttile -- nexttile
% returns a plain Cartesian axes, incompatible with polarhistogram).
showRose = show([show.NDirectional] >= opts.minDirForRose);
if isempty(showRose), return; end
nRose = numel(showRose);
nCol = min(6, nRose); nRow = ceil(nRose/nCol);
f2 = figure('Visible','off','Position',[0 0 220*nCol 220*nRow+60], 'Color', 'w');
margin = 0.03; topMargin = 0.12; % reserve space so the (now 2-line) sgtitle doesn't overlap the top row
tileW = (1-margin)/nCol; tileH = (1-margin-topMargin)/nRow;
edges = deg2rad(0:45:360); % 8 equal bins spanning the full circle (polarhistogram requires edges within a 2*pi span)
for a = 1:nRose
    p = showRose(a);
    row = ceil(a/nCol); col = a - (row-1)*nCol;
    pos = [margin + (col-1)*tileW, 1 - topMargin - row*tileH, tileW*0.85, tileH*0.8];
    ax = polaraxes(f2, 'Units', 'normalized', 'Position', pos);
    hold(ax, 'on');
    polarhistogram(ax, deg2rad(p.DirectionalPrefDirs), edges, 'FaceColor', [0.2 0.4 0.8], 'FaceAlpha', 0.7);
    maxCount = max(1, max(histcounts(deg2rad(p.DirectionalPrefDirs), edges)));
    polarplot(ax, deg2rad([p.DominantPrefDir p.DominantPrefDir]), [0 p.Concentration*maxCount], 'r:', 'LineWidth', 1.5);
    title(ax, sprintf('%s (n=%d, %d dir)', p.Area, p.NChannels, p.NDirectional), 'FontSize', 8, 'Interpreter', 'none');
end
st = sgtitle(f2, {'Preferred-direction DISTRIBUTION by area (8 bins, not a single value) -- red dotted line = resultant vector mean, for reference only', ...
    'Distribution shape reflects which part of each area a penetration happened to sample, not necessarily an intrinsic property of the area (see PIPELINE REPORT sec 6.7)'}, 'FontSize', 9, 'Interpreter', 'none');
st.Color = 'k'; % forceLightTheme's generic text sweep doesn't reliably reach sgtitle's internal text object
forceLightTheme(f2);
saveas(f2, fullfile(opts.figDir, 'pooled', 'pooled_tuning_by_area_rose.png'));
close(f2);
fprintf('saved %s\n', fullfile(opts.figDir, 'pooled', 'pooled_tuning_by_area_rose.png'));
end
