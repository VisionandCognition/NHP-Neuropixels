function tuning = computeDirectionTuning(sessions, screening, opts)
% COMPUTEDIRECTIONTUNING For every channel flagged responsive by
% screenResponsiveChannels, compute an 8-direction saccade tuning curve,
% preferred direction, tuning strength, and significance (permutation
% test on the circular resultant vector length).
%
% tuning = computeDirectionTuning(sessions, screening, opts)
%   sessions  : struct array from loadAllSessions()
%   screening : struct array from screenResponsiveChannels()
%   opts.respWindow : [t0 t1] s relative to targ_start used for the
%                     tuning response measure (default [0 0.3], the
%                     peri-saccadic window)
%   opts.baseWindow : baseline window (default [-0.9 -0.6])
%   opts.nPerm      : permutations for significance test (default 1000)
%   opts.figDir     : figure output dir

if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'respWindow', [0 0.3]);
opts = setDefault(opts, 'baseWindow', [-0.9 -0.6]);
opts = setDefault(opts, 'nPerm', 1000);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures'));

if ~isfolder(opts.figDir), mkdir(opts.figDir); end

dirsDeg = 0:45:315;
dirsRad = deg2rad(dirsDeg);

tuning = struct([]);
rowIdx = 0;

for s = 1:numel(sessions)
    S = sessions(s);
    sc = screening(s);
    assert(strcmp(S.Day, sc.Day) && S.RunN == sc.RunN, 'session/screening mismatch');

    isCorrect = strcmp(S.trialTable.outcome, 'correct') & ~isnan(S.trialTable.direction_deg);
    respIdx = S.tb >= opts.respWindow(1) & S.tb < opts.respWindow(2);
    baseIdx = S.tb >= opts.baseWindow(1) & S.tb < opts.baseWindow(2);

    respTrials = squeeze(mean(S.MUA(:, isCorrect, respIdx), 3)); % nChan x nTrial
    baseTrials = squeeze(mean(S.MUA(:, isCorrect, baseIdx), 3));
    dResp = respTrials - baseTrials; % baseline-subtracted response, nChan x nTrial

    dirDeg = S.trialTable.direction_deg(isCorrect);
    dirGroup = round(dirDeg/45) + 1; % 1..8

    respChans = find(sc.responsive)';
    for c = respChans
        r = dResp(c,:); % 1 x nTrial

        meanByDir = nan(1,8); semByDir = nan(1,8); nByDir = nan(1,8);
        for d = 1:8
            v = r(dirGroup==d);
            meanByDir(d) = mean(v);
            semByDir(d) = std(v)/sqrt(numel(v));
            nByDir(d) = numel(v);
        end

        % vector-sum preferred direction & resultant length (rectified weights)
        w = max(meanByDir, 0);
        if sum(w) > 0
            vx = sum(w.*cosd(dirsDeg))/sum(w);
            vy = sum(w.*sind(dirsDeg))/sum(w);
        else
            vx = 0; vy = 0;
        end
        prefDir = mod(atan2d(vy,vx), 360);
        resultantLen = hypot(vx,vy);

        % permutation test: shuffle direction labels, recompute resultant length
        permR = nan(opts.nPerm,1);
        for p = 1:opts.nPerm
            permGroup = dirGroup(randperm(numel(dirGroup)));
            mBy = nan(1,8);
            for d = 1:8
                mBy(d) = mean(r(permGroup==d));
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
        pPerm = (1 + sum(permR >= resultantLen)) / (opts.nPerm + 1);

        % also a simple Kruskal-Wallis across direction groups
        pKW = kruskalwallis(r', dirGroup', 'off');

        rowIdx = rowIdx + 1;
        tuning(rowIdx).Day = S.Day;
        tuning(rowIdx).RunN = S.RunN;
        tuning(rowIdx).channel = c;
        tuning(rowIdx).depthIdx = c; % already depth-ordered (S.MUA rows)
        tuning(rowIdx).meanByDir = meanByDir;
        tuning(rowIdx).semByDir = semByDir;
        tuning(rowIdx).nByDir = nByDir;
        tuning(rowIdx).prefDir = prefDir;
        tuning(rowIdx).resultantLen = resultantLen;
        tuning(rowIdx).pPerm = pPerm;
        tuning(rowIdx).pKW = pKW;
        tuning(rowIdx).tuned = pPerm < 0.05;
        % category distinguishes a clean single-preferred-direction
        % response ('directional', what pPerm alone tests) from a channel
        % with real but non-unimodal direction structure ('complex' --
        % e.g. bimodal/opposite-direction responses, which cancel out in
        % the vector-sum test but show up in the omnibus Kruskal-Wallis
        % test) from no direction effect at all ('untuned'). See
        % classifyTuning.m header for the full reasoning.
        tuning(rowIdx).category = classifyTuning(pPerm, pKW);
    end

    plotSessionTuning(S, sc, tuning, dirsDeg, opts);
    plotPopulationSummary(S, tuning, dirsDeg, opts);
end

outFile = fullfile(fileparts(opts.figDir), 'tuning_results.mat');
save(outFile, 'tuning', 'opts');
fprintf('saved %s (%d channel x session tuning rows, %d significantly tuned)\n', ...
    outFile, numel(tuning), sum([tuning.tuned]));

end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function plotSessionTuning(S, sc, tuning, dirsDeg, opts)
rows = find(strcmp({tuning.Day}, S.Day) & [tuning.RunN] == S.RunN);
if isempty(rows), return; end

% --- tuning_examples.png: ALL responsive channels (not just tuned, not
% capped at an arbitrary count -- that was the source of "why is the
% count different per run": it was min(nTuned,12), which just tracks
% however many happened to be tuned that day). Sorted by depth (channel
% index), matching the overview figure's channel ordering.
nShow = numel(rows);
[~, depthOrder] = sort([tuning(rows).depthIdx]);
rowsByDepth = rows(depthOrder);

f = figure('Visible','off','Position',[0 0 2200 2000], 'Color', 'w');
nCol = ceil(sqrt(nShow*1.3)); nRow = ceil(nShow/nCol);
tl = tiledlayout(f, nRow, nCol, 'TileSpacing', 'none', 'Padding', 'compact');
nDirectional = sum(strcmp({tuning(rows).category}, 'directional'));
nComplex = sum(strcmp({tuning(rows).category}, 'complex'));
for i = 1:nShow
    r = rowsByDepth(i);
    ax = nexttile(tl); hold(ax, 'on');
    v = tuning(r).meanByDir;
    normFactor = max(range(v), eps);
    vn = (v - min(v)) / normFactor; % per-channel min-max, just for compact display
    % SEM scaled by the SAME normalization factor -- not just cosmetic:
    % a curve can look sharply "peaked" once stretched to fill [0,1]
    % regardless of whether that peak is bigger than the trial-to-trial
    % noise or not (min-max normalization erases absolute effect size).
    % Plotting error bars in the same normalized units directly shows
    % whether a visually-obvious-looking peak is actually distinguishable
    % from noise -- e.g. a channel can fail the permutation test at
    % p=0.06 with a real-looking peak simply because its error bars at
    % the peak and trough overlap once you look at them (verified
    % directly against several borderline "untuned" channels before
    % making this change).
    semN = tuning(r).semByDir / normFactor;
    [~, col] = classifyTuning(tuning(r).pPerm, tuning(r).pKW);
    errorbar(ax, [dirsDeg dirsDeg(1)+360], [vn vn(1)], [semN semN(1)], '-o', ...
        'Color', col, 'MarkerFaceColor', col, 'LineWidth', 1, 'MarkerSize', 2, 'CapSize', 2);
    ylim(ax, [-0.5 1.5]); % wider than [-0.1 1.1] -- SEM whiskers can extend past the normalized curve's own min/max
    xticks(ax, []); yticks(ax, []);
    if strcmp(tuning(r).category, 'untuned')
        set(ax, 'XColor', [0.85 0.85 0.85], 'YColor', [0.85 0.85 0.85], 'LineWidth', 0.5, 'Box', 'on');
    else
        set(ax, 'XColor', col, 'YColor', col, 'LineWidth', 1.5, 'Box', 'on');
    end
    title(ax, sprintf('%d', tuning(r).channel), 'FontSize', 6, 'FontWeight', 'normal');
end
title(tl, sprintf(['%s run-%03d: ALL responsive channels (n=%d) -- ' ...
    'red = directional (n=%d), orange = complex/non-unimodal (n=%d), gray = untuned'], ...
    S.Day, S.RunN, nShow, nDirectional, nComplex));
forceLightTheme(f);
saveas(f, figSavePath(opts.figDir, 'tuning_examples', S.Day, S.RunN));
close(f);

% --- tuning_matrix.png: all responsive channels, sorted by preferred
% direction, with a *bounded* per-channel normalization (min-max to
% [0,1]) instead of dividing by the row max -- the old version divided
% raw (signed) baseline-subtracted responses by their own max, so a
% weakly-modulated/noisy channel with a tiny positive max but a much
% larger negative dip produced huge out-of-range ratios, blowing out the
% shared color axis and washing out every other (well-behaved) row.
% Min-max keeps every row in [0,1] regardless of sign/magnitude, so the
% shared colorbar is never dominated by a handful of outlier channels.
% A tuning-strength strip and a tuned/not marker are shown alongside so
% amplitude information isn't lost in the normalization.
if numel(rows) >= 2
    prefs = [tuning(rows).prefDir];
    [~, order] = sort(prefs);
    allRows = rows(order);
    M = cat(1, tuning(allRows).meanByDir);
    Mnorm = (M - min(M,[],2)) ./ max(range(M,2), eps);
    strength = [tuning(allRows).resultantLen]';
    [~, catColor] = classifyTuning([tuning(allRows).pPerm], [tuning(allRows).pKW]);

    f2 = figure('Visible','off','Position',[0 0 1000 700], 'Color', 'w');
    tl2 = tiledlayout(f2, 1, 10, 'TileSpacing', 'compact', 'Padding', 'compact');

    axM = nexttile(tl2, 1, [1 6]);
    imagesc(axM, dirsDeg, 1:numel(allRows), Mnorm);
    set(axM, 'CLim', [0 1], 'Color', 'w');
    xlabel(axM, 'Saccade direction (deg)'); ylabel(axM, 'Channel (sorted by preferred direction)');
    title(axM, sprintf('%s run-%03d: tuning matrix (per-channel min-max)', S.Day, S.RunN));
    colorbar(axM);

    axS = nexttile(tl2, 7, [1 2]);
    imagesc(axS, 1, 1:numel(allRows), strength);
    set(axS, 'Color', 'w', 'XTick', []);
    title(axS, 'strength', 'FontSize', 8);
    cbS = colorbar(axS); cbS.Label.String = 'resultant length';

    axT = nexttile(tl2, 9, [1 2]);
    image(axT, 1, 1:numel(allRows), reshape(catColor, [numel(allRows) 1 3]));
    set(axT, 'Color', 'w', 'XTick', []);
    title(axT, {'category', '(red=dir,orange=cx)'}, 'FontSize', 7);

    forceLightTheme(f2);
    saveas(f2, figSavePath(opts.figDir, 'tuning_matrix', S.Day, S.RunN));
    close(f2);
end
end

% ------------------------------------------------------------------
function plotPopulationSummary(S, tuning, dirsDeg, opts)
% Per-SESSION summary (not pooled across sessions): each session is an
% independent penetration, often a different burr hole/target structure
% (see sessionAnatomyInfo.m), and the probe is not seated to the same
% depth across sessions -- pooling channel depth index across sessions
% (as an earlier version of this figure did) implied a false
% cross-session comparability. Real depth (um from tip) is used here
% instead of abstract channel index, meaningful within this one session.
rows = find(strcmp({tuning.Day}, S.Day) & [tuning.RunN] == S.RunN);
if isempty(rows), return; end
tuned = [tuning(rows).tuned];
anat = sessionAnatomyInfo(S.Day);

depthUm = nan(1, numel(rows));
if isfield(S, 'CHdepthUm')
    depthUm = S.CHdepthUm([tuning(rows).depthIdx])';
end

f = figure('Visible','off','Position',[0 0 1400 500], 'Color', 'w');
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = polaraxes(tl);
ax1.Layout.Tile = 1;
if any(tuned)
    polarhistogram(ax1, deg2rad([tuning(rows(tuned)).prefDir]), 16);
end
title(ax1, sprintf('Preferred direction (n=%d/%d tuned)', sum(tuned), numel(rows)));

ax2 = nexttile(tl);
histogram(ax2, [tuning(rows).resultantLen], 20);
set(ax2, 'Color', 'w');
xlabel(ax2, 'Resultant vector length (tuning strength)'); ylabel(ax2, '# channels');
title(ax2, 'Tuning strength, all responsive channels');

ax3 = nexttile(tl);
scatter(ax3, depthUm/1000, [tuning(rows).prefDir], 15, [tuning(rows).resultantLen], 'filled');
set(ax3, 'Color', 'w');
c = colorbar(ax3); c.Label.String = 'tuning strength';
xlabel(ax3, 'Real electrode depth (mm from tip)'); ylabel(ax3, 'Preferred direction (deg)');
title(ax3, 'Preferred direction vs. real depth (this session only)');

title(tl, sprintf('%s run-%03d population summary -- burr hole %s (%s)', ...
    S.Day, S.RunN, anat.burrHole, strjoin(anat.structures, ' -> ')), 'Interpreter', 'none');
forceLightTheme(f);
saveas(f, figSavePath(opts.figDir, 'population_summary', S.Day, S.RunN));
close(f);
end
