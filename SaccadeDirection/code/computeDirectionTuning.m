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
    end

    plotSessionTuning(S, sc, tuning, dirsDeg, opts);
end

outFile = fullfile(fileparts(opts.figDir), 'tuning_results.mat');
save(outFile, 'tuning', 'opts');
fprintf('saved %s (%d channel x session tuning rows, %d significantly tuned)\n', ...
    outFile, numel(tuning), sum([tuning.tuned]));

plotPopulationSummary(tuning, dirsDeg, opts);

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

tunedRows = rows([tuning(rows).tuned]);
nShow = min(numel(tunedRows), 12);
if nShow == 0, return; end

f = figure('Visible','off','Position',[0 0 1400 900]);
nCol = ceil(sqrt(nShow)); nRow = ceil(nShow/nCol);
for i = 1:nShow
    r = tunedRows(i);
    subplot(nRow, nCol, i);
    theta = deg2rad([dirsDeg dirsDeg(1)]);
    rho = [tuning(r).meanByDir, tuning(r).meanByDir(1)];
    polarplot(theta, max(rho,0), '-o', 'LineWidth', 1.5);
    title(sprintf('ch%d PD=%.0f R=%.2f p=%.3f', tuning(r).channel, ...
        tuning(r).prefDir, tuning(r).resultantLen, tuning(r).pPerm), 'FontSize', 8);
end
sgtitle(sprintf('%s run-%03d: direction-tuned channels', S.Day, S.RunN));
saveas(f, fullfile(opts.figDir, sprintf('%s_run-%03d_tuning_examples.png', S.Day, S.RunN)));
close(f);

% tuning matrix heatmap sorted by preferred direction (all tuned channels this session)
if numel(tunedRows) >= 2
    allRows = rows;
    prefs = [tuning(allRows).prefDir];
    [~, order] = sort(prefs);
    allRows = allRows(order);
    M = cat(1, tuning(allRows).meanByDir);
    Mnorm = M ./ max(M, [], 2);

    f2 = figure('Visible','off','Position',[0 0 800 600]);
    imagesc(dirsDeg, 1:numel(allRows), Mnorm);
    xlabel('Saccade direction (deg)'); ylabel('Channel (sorted by preferred direction)');
    title(sprintf('%s run-%03d: tuning matrix (normalized)', S.Day, S.RunN));
    colorbar;
    saveas(f2, fullfile(opts.figDir, sprintf('%s_run-%03d_tuning_matrix.png', S.Day, S.RunN)));
    close(f2);
end
end

% ------------------------------------------------------------------
function plotPopulationSummary(tuning, dirsDeg, opts)
if isempty(tuning), return; end
tuned = [tuning.tuned];

f = figure('Visible','off','Position',[0 0 1200 500]);

subplot(1,3,1);
polarhistogram(deg2rad([tuning(tuned).prefDir]), 16);
title(sprintf('Preferred direction (n=%d tuned channels)', sum(tuned)));

subplot(1,3,2);
histogram([tuning.resultantLen], 20);
xlabel('Resultant vector length (tuning strength)'); ylabel('# channels');
title('Tuning strength, all screened channels');

subplot(1,3,3);
scatter([tuning.depthIdx], [tuning.prefDir], 15, [tuning.resultantLen], 'filled');
c = colorbar; c.Label.String = 'tuning strength';
xlabel('Channel depth index'); ylabel('Preferred direction (deg)');
title('Preferred direction vs. probe depth');

sgtitle(sprintf('Population summary across %d sessions', numel(unique({tuning.Day}))));
saveas(f, fullfile(opts.figDir, 'population_tuning_summary.png'));
close(f);
end
