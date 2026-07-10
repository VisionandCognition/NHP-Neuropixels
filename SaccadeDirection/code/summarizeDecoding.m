function summary = summarizeDecoding()
% SUMMARIZEDECODING Quantifies decodeDirection.m's results across ALL 12
% sessions: per-session full-population decode accuracy per epoch, and a
% pooled population-size-vs-accuracy curve (averaged across sessions at
% each common population size, standard error across sessions), plus a
% pooled summary figure.
%
% summary = summarizeDecoding()

load(fullfile(saccadeDataDir(), 'decode_results.mat'), 'decodeResults');
epochNames = decodeResults(1).epochNames;
nEpoch = numel(epochNames);

summary = struct([]);
fprintf('Day\t\tnResp\tnTrial\tvisual%%\tplanning%%\texecution%%\n');
for i = 1:numel(decodeResults)
    r = decodeResults(i);
    fullAcc = squeeze(mean(r.popAcc(end,:,:), 3)) * 100; % 1 x nEpoch
    summary(i).Day = r.Day;
    summary(i).RunN = r.RunN;
    summary(i).nResp = r.nResp;
    summary(i).nTrial = r.nTrial;
    summary(i).fullPopAcc = fullAcc;
    fprintf('%s\t%d\t%d\t%.1f\t%.1f\t\t%.1f\n', r.Day, r.nResp, r.nTrial, fullAcc(1), fullAcc(2), fullAcc(3));
end

allFull = cat(1, summary.fullPopAcc); % nSession x nEpoch
fprintf('\n=== Pooled across %d sessions (full responsive population, chance=12.5%%) ===\n', numel(summary));
for e = 1:nEpoch
    fprintf('%s: mean %.1f%%, SEM %.1f%%, range [%.1f, %.1f]%%, %d/%d sessions above chance\n', ...
        epochNames{e}, mean(allFull(:,e)), std(allFull(:,e))/sqrt(size(allFull,1)), ...
        min(allFull(:,e)), max(allFull(:,e)), sum(allFull(:,e) > 12.5), size(allFull,1));
end

% pooled population-size curve: average accuracy across sessions at each
% common population size (only where that session's own sweep reached
% that size -- smaller-nResp sessions drop out of the curve at high N).
% Restricted to the standard geometric sizes (decodeDirection.m's default
% opts.popSizes) -- each session's *_decode_population.png ALSO includes
% one extra point at its own idiosyncratic full-responsive-population
% size (e.g. 156, 271, 361...), which is not comparable across sessions
% and would otherwise make most of the pooled curve's tail reflect just
% one single session (not a real pooled average). Each session's own
% full-population number is reported separately above instead.
standardPopSizes = [1 2 4 8 16 32 64 128 256];
allPopSizes = intersect(standardPopSizes, unique([decodeResults.popSizes]));
pooledMean = nan(numel(allPopSizes), nEpoch);
pooledSEM = nan(numel(allPopSizes), nEpoch);
nSessAtSize = nan(numel(allPopSizes), 1);
for p = 1:numel(allPopSizes)
    P = allPopSizes(p);
    vals = nan(numel(decodeResults), nEpoch);
    for i = 1:numel(decodeResults)
        idx = find(decodeResults(i).popSizes == P, 1);
        if ~isempty(idx)
            vals(i,:) = squeeze(mean(decodeResults(i).popAcc(idx,:,:), 3)) * 100;
        end
    end
    nSessAtSize(p) = sum(~isnan(vals(:,1)));
    pooledMean(p,:) = mean(vals, 1, 'omitnan');
    pooledSEM(p,:) = std(vals, 0, 1, 'omitnan') ./ sqrt(max(nSessAtSize(p),1));
end

plotPooled(allPopSizes, pooledMean, pooledSEM, nSessAtSize, epochNames);
end

% ------------------------------------------------------------------
function plotPooled(popSizes, pooledMean, pooledSEM, nSessAtSize, epochNames)
colors = {[0.2 0.4 0.8], [0.85 0.55 0.1], [0.75 0.1 0.1]};
f = figure('Visible','off','Position',[0 0 750 580], 'Color', 'w');
ax = axes(f); hold(ax, 'on');
hs = gobjects(1,3);
for e = 1:3
    m = pooledMean(:,e)'; se = pooledSEM(:,e)';
    ok = ~isnan(m);
    fill(ax, [popSizes(ok) fliplr(popSizes(ok))], [m(ok)+se(ok) fliplr(m(ok)-se(ok))], colors{e}, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    hs(e) = plot(ax, popSizes(ok), m(ok), '-o', 'Color', colors{e}, 'LineWidth', 1.8, 'MarkerFaceColor', colors{e}, 'MarkerSize', 4);
end
yline(ax, 12.5, 'k:', 'LineWidth', 1);
set(ax, 'XScale', 'log', 'Color', 'w');
xlabel(ax, 'Population size (# channels, log scale)'); ylabel(ax, 'Decode accuracy (%, mean +/- SEM across sessions)');
ylim(ax, [0 60]);
legend(ax, hs, epochNames, 'Location', 'northwest', 'TextColor', 'k', 'Color', 'w');
title(ax, sprintf('Pooled across 12 sessions: direction decoding vs. population size\n(nearest-centroid, LOO-CV, chance=12.5%%; N sessions contributing per size: %s)', ...
    strjoin(arrayfun(@(n) num2str(n), nSessAtSize, 'UniformOutput', false), ',')), 'FontSize', 9);
forceLightTheme(f);
figDir = fullfile(saccadeDataDir(), 'figures');
if ~isfolder(fullfile(figDir,'pooled')), mkdir(fullfile(figDir,'pooled')); end
saveas(f, fullfile(figDir, 'pooled', 'pooled_decode_population.png'));
close(f);
fprintf('saved %s\n', fullfile(figDir, 'pooled', 'pooled_decode_population.png'));
end
