function plotDirectionTimecourses(sessions, screening, eyeSummary, opts)
% PLOTDIRECTIONTIMECOURSES Per-channel mean MUA timecourse split by the 8
% saccade directions (correct trials only), two alignments:
%   - go-cue-locked: the trial's own time base (S.tb), as used everywhere
%     else in this pipeline.
%   - saccade-onset-locked: each trial's MUA re-centered on its own
%     detected saccade onset (from analyzeEyeData.m), via linear
%     interpolation onto a common relative time axis -- necessary because
%     onset time varies trial to trial, so a simple index shift would mix
%     up different absolute sample times across trials.
%
% One figure PER RESPONSIVE CHANNEL per session per alignment (this is
% deliberately a lot of files -- ~3200 responsive channels across all 12
% sessions, x2 alignments), so figures are organized into subfolders:
%   <saccadeDataDir()>/figures/timecourses/<Day>_run-<RunN>/gocue/ch<NNN>.png
%   <saccadeDataDir()>/figures/timecourses/<Day>_run-<RunN>/saclocked/ch<NNN>.png
% <NNN> is the channel INDEX (1 = deepest / probe tip, matching every
% other figure in this pipeline), not a physical electrode ID.
%
% Each direction's mean trace is smoothed (zero-phase Gaussian, default
% 20ms) before plotting -- purely a display transform applied to the
% already-averaged trace, not to the underlying MUA/S.MUA or any other
% analysis in this pipeline (tuning/decoding/screening all use windowed
% averages and are unaffected). Convolution commutes with averaging, so
% smoothing the final per-direction mean is equivalent to smoothing every
% trial first, at a fraction of the cost.
%
% plotDirectionTimecourses(sessions, screening, eyeSummary, opts)
%   sessions   : struct array from loadAllSessions()
%   screening  : struct array from screenResponsiveChannels() (for .responsive)
%   eyeSummary : struct array from analyzeEyeData() (for saccade onsets)
%   opts.gocueWindow : display window, s, relative to go-cue (default [-0.9 1.0])
%   opts.sacWindow   : display window, s, relative to saccade onset (default [-0.5 0.8])
%   opts.smoothMs    : Gaussian smoothing window, ms, applied to each
%                      direction's mean trace before plotting (default 50
%                      -- 20ms tried first and was barely visible against
%                      this 200Hz-envelope MUA's high-frequency content)
%   opts.figDir      : default <saccadeDataDir()>/figures/timecourses

if nargin < 4, opts = struct(); end
opts = setDefault(opts, 'gocueWindow', [-0.9 1.0]);
opts = setDefault(opts, 'sacWindow', [-0.5 0.8]);
opts = setDefault(opts, 'smoothMs', 50);
opts = setDefault(opts, 'figDir', fullfile(saccadeDataDir(), 'figures', 'timecourses'));
if ~isfolder(opts.figDir), mkdir(opts.figDir); end

dirsDeg = 0:45:315;
dirLabels = arrayfun(@(x) sprintf('%d deg', x), dirsDeg, 'UniformOutput', false);
% Same fixed 8-direction palette used in analyzeEyeData.m's trajectory
% plots (hsv(8) puts a nearly-invisible pale yellow-green at 90 deg) --
% reused here so direction colors are consistent across every figure type
% in this pipeline.
cmap = [0.89 0.10 0.11; 1.00 0.50 0.00; 0.60 0.60 0.00; 0.20 0.60 0.20; ...
        0.00 0.60 0.60; 0.12 0.47 0.71; 0.58 0.00 0.83; 0.89 0.10 0.55];

for s = 1:numel(sessions)
    S = sessions(s);
    isCorrect = strcmp(S.trialTable.outcome, 'correct') & ~isnan(S.trialTable.direction_deg);
    dirGroup = round(S.trialTable.direction_deg/45) + 1;
    correctIdx = find(isCorrect);
    respChans = find(screening(s).responsive);
    nChan = numel(respChans);

    dayDir = fullfile(opts.figDir, sprintf('%s_run-%03d', S.Day, S.RunN));
    gocueDir = fullfile(dayDir, 'gocue');
    sacDir = fullfile(dayDir, 'saclocked');
    if ~isfolder(gocueDir), mkdir(gocueDir); end

    % --- go-cue-locked window (shared time base, no interpolation needed)
    gcIdx = S.tb >= opts.gocueWindow(1) & S.tb <= opts.gocueWindow(2);
    tbGC = S.tb(gcIdx) * 1000;
    dt = median(diff(S.tb)); % shared sampling interval, used for both alignments' smoothing window
    smoothWin = max(1, round(opts.smoothMs/1000/dt));

    % --- saccade-onset-locked: interpolate each usable correct trial onto
    % a common relative time axis, once per session (not per channel).
    esIdx = find(strcmp({eyeSummary.Day}, S.Day) & [eyeSummary.RunN] == S.RunN);
    haveOnsets = ~isempty(esIdx) && S.eyeCalib.valid;
    useTrials = []; snippet = []; relTbMs = [];
    if haveOnsets
        if ~isfolder(sacDir), mkdir(sacDir); end
        es = eyeSummary(esIdx);
        onsetTime = es.onsetTime;
        relTb = opts.sacWindow(1):dt:opts.sacWindow(2);
        relTbMs = relTb * 1000;
        useTrials = correctIdx(~isnan(onsetTime(correctIdx)));
        snippet = nan(size(S.MUA,1), numel(relTb), numel(useTrials));
        for k = 1:numel(useTrials)
            t = useTrials(k);
            absT = onsetTime(t) + relTb;
            X = squeeze(S.MUA(:,t,:))'; % nTime x nChan
            snippet(:,:,k) = interp1(S.tb, X, absT, 'linear')'; % nChan x nRelT
        end
    end

    % Reuse ONE figure/axes per alignment across all channels in this
    % session (cla + redraw, not figure()/close() per channel) --
    % per-channel figure creation/destruction dominated runtime in initial
    % testing (~2.6s/channel-pair, which would be several hours across all
    % 12 sessions' ~3200 responsive channels). forceLightTheme still runs
    % every save since new Text/Line objects are created each iteration.
    f1 = figure('Visible','off','Position',[0 0 560 400], 'Color', 'w');
    ax1 = axes(f1);
    if haveOnsets
        f2 = figure('Visible','off','Position',[0 0 560 400], 'Color', 'w');
        ax2 = axes(f2);
    end

    for ci = 1:nChan
        c = respChans(ci);
        depthStr = '';
        if isfield(S, 'CHdepthUm'), depthStr = sprintf(' (%.2fmm)', S.CHdepthUm(c)/1000); end

        cla(ax1); hold(ax1, 'on');
        for d = 1:8
            idx = correctIdx(dirGroup(correctIdx) == d);
            if isempty(idx), continue; end
            tr = squeeze(mean(S.MUA(c, idx, gcIdx), 2, 'omitnan'));
            tr = smoothdata(tr, 'gaussian', smoothWin);
            plot(ax1, tbGC, tr, 'Color', cmap(d,:), 'LineWidth', 1.3);
        end
        xline(ax1, 0, 'k:');
        xlabel(ax1, 'Time from go-cue (ms)'); ylabel(ax1, 'MUA (a.u.)');
        title(ax1, sprintf('%s run-%03d ch%d%s: go-cue-locked', S.Day, S.RunN, c, depthStr), 'FontSize', 9);
        legend(ax1, dirLabels, 'Location', 'eastoutside', 'FontSize', 7, 'TextColor', 'k', 'Color', 'w');
        set(ax1, 'Color', 'w');
        forceLightTheme(f1);
        saveas(f1, fullfile(gocueDir, sprintf('ch%03d.png', c)));

        if haveOnsets
            cla(ax2); hold(ax2, 'on');
            for d = 1:8
                idxMask = dirGroup(useTrials) == d;
                if ~any(idxMask), continue; end
                tr = squeeze(mean(snippet(c,:,idxMask), 3, 'omitnan'));
                tr = smoothdata(tr, 'gaussian', smoothWin);
                plot(ax2, relTbMs, tr, 'Color', cmap(d,:), 'LineWidth', 1.3);
            end
            xline(ax2, 0, 'k:');
            xlabel(ax2, 'Time from saccade onset (ms)'); ylabel(ax2, 'MUA (a.u.)');
            title(ax2, sprintf('%s run-%03d ch%d%s: saccade-onset-locked', S.Day, S.RunN, c, depthStr), 'FontSize', 9);
            legend(ax2, dirLabels, 'Location', 'eastoutside', 'FontSize', 7, 'TextColor', 'k', 'Color', 'w');
            set(ax2, 'Color', 'w');
            forceLightTheme(f2);
            saveas(f2, fullfile(sacDir, sprintf('ch%03d.png', c)));
        end
    end
    close(f1);
    if haveOnsets, close(f2); end

    if haveOnsets
        fprintf('%s run-%03d: saved %d go-cue-locked + %d saccade-onset-locked timecourse figures\n', S.Day, S.RunN, nChan, nChan);
    else
        fprintf('%s run-%03d: saved %d go-cue-locked timecourse figures (saccade-onset-locked skipped: no valid eye calibration)\n', S.Day, S.RunN, nChan);
    end
end
end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end
