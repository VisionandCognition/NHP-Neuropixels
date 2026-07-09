function generateAnatomyTemplate(sessions, screening, opts)
% GENERATEANATOMYTEMPLATE Builds a fillable template for manually
% assigning anatomical area names to electrode depth ranges, by
% combining: (1) real per-channel depth (S.CHdepthUm), (2) this
% pipeline's own recorded activity (responsiveness, tuning category,
% preferred direction, tuning strength), and (3) the planned
% burr-hole/trajectory info already in sessionAnatomyInfo.m. This does
% NOT assign area names itself -- it proposes candidate channel clusters
% and gives you the data to compare against an atlas, with blank columns
% to fill in.
%
% Two CSVs are written per call (all 12 sessions each), to
% <REPO_ROOT>/SaccadeDirection/anatomy/:
%
%   anatomy_template_clusters.csv -- the document to fill out. One row
%     per (session, contiguous cluster of channels whose recorded
%     activity looks similar): channel/depth range, % responsive, %
%     directional/complex, dominant preferred direction (circular mean +
%     concentration, only meaningful if PctDirectional is reasonably
%     high), mean tuning strength, planned burr-hole/structure list for
%     context, and blank AssignedArea/Confidence/Notes columns.
%
%   anatomy_template_channels.csv -- full per-channel reference (every
%     channel, not just responsive ones), with a ClusterID column linking
%     back to the cluster table -- use this if you want to check/override
%     a cluster boundary at finer grain, or spot-check individual
%     channels near a boundary.
%
% CLUSTER DETECTION METHOD (transparent, no exotic toolbox dependency):
% bin channels into fixed-size bins (opts.binSize, default 20 channels)
% and summarize each bin's responsiveness rate and resultant-length-
% weighted circular mean preferred direction, then greedily merge a bin
% into the cluster to its left if both are "similar enough"
% (opts.respMergeThresh responsiveness-rate difference AND
% opts.dirMergeThreshDeg circular preferred-direction difference).
% Binning first (rather than only flagging sharp jumps) matters here: at
% least one session (20260312) shows preferred direction drifting
% SMOOTHLY and continuously with depth (no sharp jump anywhere, checked
% directly), which a pure jump-detector collapses into one meaningless
% whole-session "cluster" -- binning keeps a continuous gradient visible
% as a sequence of clusters instead of hiding it. This is a heuristic
% proposal, not a certified segmentation -- it exists to save you from
% segmenting ~3200 channels by eye, not to replace your judgment;
% adjust/merge/split clusters freely when filling in the template.
%
% generateAnatomyTemplate(sessions, screening, opts)
%   sessions  : struct array from loadAllSessions()
%   screening : struct array from screenResponsiveChannels()
%   opts.binSize            : base bin size, channels (default 20)
%   opts.respMergeThresh    : max responsiveness-rate difference to merge
%                             adjacent bins into one cluster (default 0.15)
%   opts.dirMergeThreshDeg  : max circular preferred-direction difference
%                             to merge adjacent bins (default 30)
%   opts.minDirCountInBin   : min directional channels in a bin before its
%                             preferred-direction estimate is trusted at
%                             all (default 3)
%   opts.outDir             : default <SaccadeDirection repo root>/anatomy

if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'binSize', 20);
opts = setDefault(opts, 'respMergeThresh', 0.15);
opts = setDefault(opts, 'dirMergeThreshDeg', 30);
opts = setDefault(opts, 'minDirCountInBin', 3);
opts = setDefault(opts, 'outDir', fullfile(fileparts(fileparts(mfilename('fullpath'))), 'anatomy'));
if ~isfolder(opts.outDir), mkdir(opts.outDir); end

load(fullfile(saccadeDataDir(), 'tuning_results.mat'), 'tuning');

clusterRows = {};
channelRows = {};

for s = 1:numel(sessions)
    S = sessions(s);
    sc = screening(s);
    assert(strcmp(S.Day, sc.Day) && S.RunN == sc.RunN, 'session/screening mismatch');
    anat = sessionAnatomyInfo(S.Day);
    nChan = numel(sc.responsive);

    tRows = find(strcmp({tuning.Day}, S.Day) & [tuning.RunN] == S.RunN);
    tByChan = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for r = tRows
        tByChan(tuning(r).channel) = tuning(r);
    end

    depthMm = nan(nChan,1);
    if isfield(S, 'CHdepthUm'), depthMm = S.CHdepthUm(:)/1000; end
    responsive = sc.responsive(:);
    category = repmat({'untested'}, nChan, 1);
    prefDir = nan(nChan,1);
    resultantLen = zeros(nChan,1);
    for c = 1:nChan
        if responsive(c) && isKey(tByChan, c)
            tr = tByChan(c);
            category{c} = tr.category;
            resultantLen(c) = tr.resultantLen;
            if strcmp(tr.category, 'directional')
                prefDir(c) = tr.prefDir;
            end
        elseif responsive(c)
            category{c} = 'responsive_no_tuning_row'; % shouldn't normally happen
        else
            category{c} = 'not_responsive';
        end
    end

    % --- per-channel reference rows
    for c = 1:nChan
        channelRows(end+1,:) = {S.Day, S.RunN, c, depthMm(c), responsive(c), category{c}, ...
            prefDir(c), resultantLen(c)}; %#ok<AGROW>
    end

    % --- cluster detection along channel index (already depth-ordered, 1=deepest)
    boundaries = detectClusterBoundaries(responsive, prefDir, resultantLen, opts);

    structStr = strjoin(anat.structures, ' -> ');
    for k = 1:numel(boundaries)-1
        idx = boundaries(k):(boundaries(k+1)-1);
        pctResp = 100*mean(responsive(idx));
        pctDir = 100*mean(strcmp(category(idx), 'directional'));
        pctCx = 100*mean(strcmp(category(idx), 'complex'));
        isDir = strcmp(category(idx), 'directional');
        if any(isDir)
            w = resultantLen(idx(isDir));
            pd = prefDir(idx(isDir));
            cx = sum(w.*cosd(pd)); cy = sum(w.*sind(pd));
            domPrefDir = mod(atan2d(cy,cx), 360);
            concentration = hypot(cx,cy)/max(sum(w), eps);
        else
            domPrefDir = nan; concentration = nan;
        end
        meanStrength = mean(resultantLen(idx));
        depthRange = depthMm(idx);
        depthRange = depthRange(~isnan(depthRange));
        if isempty(depthRange), dLo = nan; dHi = nan; else, dLo = min(depthRange); dHi = max(depthRange); end

        clusterRows(end+1,:) = {S.Day, S.RunN, anat.burrHole, anat.channelSelection, k, ...
            idx(1), idx(end), dLo, dHi, numel(idx), pctResp, pctDir, pctCx, ...
            domPrefDir, concentration, meanStrength, structStr, '', '', ''}; %#ok<AGROW>
    end

    fprintf('%s run-%03d: %d clusters detected (%d channels)\n', S.Day, S.RunN, numel(boundaries)-1, nChan);
end

clusterT = cell2table(clusterRows, 'VariableNames', ...
    {'Day','RunN','BurrHole','ChannelSelection','ClusterID','ChanIdxStart','ChanIdxEnd', ...
     'DepthMM_start','DepthMM_end','NChan','PctResponsive','PctDirectional','PctComplex', ...
     'DominantPrefDir_deg','PrefDirConcentration_0to1','MeanTuningStrength', ...
     'PlannedStructures_superficialToDeep','AssignedArea','Confidence','Notes'});
channelT = cell2table(channelRows, 'VariableNames', ...
    {'Day','RunN','ChannelIdx','DepthMM','Responsive','Category','PrefDir_deg','TuningStrength'});

clusterFile = fullfile(opts.outDir, 'anatomy_template_clusters.csv');
channelFile = fullfile(opts.outDir, 'anatomy_template_channels.csv');
writetable(clusterT, clusterFile);
writetable(channelT, channelFile);
fprintf('\nsaved %s (%d rows -- the document to fill out)\n', clusterFile, height(clusterT));
fprintf('saved %s (%d rows -- per-channel reference)\n', channelFile, height(channelT));
end

% ------------------------------------------------------------------
function opts = setDefault(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end

% ------------------------------------------------------------------
function boundaries = detectClusterBoundaries(responsive, prefDir, resultantLen, opts)
n = numel(responsive);
binSize = opts.binSize;
nBins = ceil(n/binSize);
binStart = nan(nBins,1); binResp = nan(nBins,1); binPrefDir = nan(nBins,1); binHasDir = false(nBins,1);
for b = 1:nBins
    idx = ((b-1)*binSize+1) : min(b*binSize, n);
    binStart(b) = idx(1);
    binResp(b) = mean(responsive(idx));
    isDir = ~isnan(prefDir(idx));
    if sum(isDir) >= opts.minDirCountInBin
        w = resultantLen(idx(isDir)); pd = prefDir(idx(isDir));
        cx = sum(w.*cosd(pd)); cy = sum(w.*sind(pd));
        if hypot(cx,cy) > eps
            binPrefDir(b) = mod(atan2d(cy,cx), 360);
            binHasDir(b) = true;
        end
    end
end

% Sequential greedy merge: a bin joins the CURRENT running cluster if
% it's similar to that cluster's running (circular-mean-updated)
% summary; otherwise it starts a new cluster. This keeps a continuous
% gradient (many bins, each only slightly different from the last) as a
% sequence of clusters, while collapsing genuinely uniform stretches.
boundaries = binStart(1);
runResp = binResp(1); runPrefDir = binPrefDir(1); runHasDir = binHasDir(1); runN = 1;
for b = 2:nBins
    similar = abs(binResp(b) - runResp) < opts.respMergeThresh;
    if similar && runHasDir && binHasDir(b)
        dAngle = abs(mod(binPrefDir(b) - runPrefDir + 180, 360) - 180);
        similar = dAngle < opts.dirMergeThreshDeg;
    end
    if similar
        runResp = (runResp*runN + binResp(b)) / (runN+1);
        if binHasDir(b)
            if runHasDir
                runPrefDir = mod(atan2d(sind(runPrefDir)*runN + sind(binPrefDir(b)), cosd(runPrefDir)*runN + cosd(binPrefDir(b))), 360);
            else
                runPrefDir = binPrefDir(b);
            end
            runHasDir = true;
        end
        runN = runN + 1;
    else
        boundaries(end+1) = binStart(b); %#ok<AGROW>
        runResp = binResp(b); runPrefDir = binPrefDir(b); runHasDir = binHasDir(b); runN = 1;
    end
end
boundaries(end+1) = n+1;
end
