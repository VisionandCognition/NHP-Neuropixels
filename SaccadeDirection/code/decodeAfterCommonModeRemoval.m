function results = decodeAfterCommonModeRemoval(sessions, screening, opts)
% DECODEAFTERCOMMONMODEREMOVAL Closes the open item from PIPELINE_REPORT
% sec 7 item 6: the 4 sessions confirmed kinematics-linked (§6.5:
% 20260310/20260320/20260323/20260324) are also among the higher-decoding
% sessions (§5.3). Is that because the shared/kinematic component is
% itself somewhat decodable (inflating accuracy), or because these
% sessions simply have more real signal? This decodes BEFORE and AFTER
% removing each epoch's own dominant shared component (PC1 of the
% direction-demeaned trial residual, same method as
% checkSharedComponentKinematics.m/§6.5), at the full responsive
% population, for all 3 epochs, for just these 4 sessions.
%
% Common-mode removal is done SEPARATELY per epoch (not once using a
% single combined window) since a kinematic/EMG-like artifact plausibly
% has a different profile in the visual epoch (well before any movement)
% than in execution (during the movement itself) -- reusing one
% epoch's PC1 for all three would conflate them. The PC1 is identified
% from the direction-DEMEANED residual (so it can't just be "real
% tuning"), but subtracted from the ORIGINAL (non-demeaned) response
% before decoding, so genuine direction-related mean differences are
% preserved and only the shared per-trial fluctuation is removed.
%
% results = decodeAfterCommonModeRemoval(sessions, screening, opts)
%   sessions, screening : as in decodeDirection.m (only the 4
%     kinematic-linked sessions are processed; others are skipped)
%   opts.baseWindow, opts.visualWindow, opts.goCueWindows : same defaults
%     as decodeDirection.m, so epochs match exactly

if nargin < 3, opts = struct(); end
opts = setDefault(opts, 'baseWindow', [-0.9 -0.6]);
opts = setDefault(opts, 'visualWindow', [-0.5 -0.35]);
opts = setDefault(opts, 'goCueWindows', struct('planning', [0 0.1], 'execution', [0.1 0.3]));

kinematicLinkedDays = {'20260310','20260320','20260323','20260324'};
epochNames = {'visual','planning','execution'};
results = struct([]);

fprintf('%-10s %-10s %8s %8s %8s\n', 'Day', 'Epoch', 'Original', 'Cleaned', 'PC1var%');
for s = 1:numel(sessions)
    S = sessions(s);
    if ~ismember(S.Day, kinematicLinkedDays), continue; end

    isCorrect = strcmp(S.trialTable.outcome, 'correct') & ~isnan(S.trialTable.direction_deg);
    dirDeg = S.trialTable.direction_deg(isCorrect);
    dirGroup = round(dirDeg(:)'/45) + 1;

    baseIdx = S.tb >= opts.baseWindow(1) & S.tb < opts.baseWindow(2);
    baseTrials = squeeze(mean(S.MUA(:, isCorrect, baseIdx), 3));
    visIdx = S.tb >= opts.visualWindow(1) & S.tb < opts.visualWindow(2);
    planIdx = S.tb >= opts.goCueWindows.planning(1) & S.tb < opts.goCueWindows.planning(2);
    execIdx = S.tb >= opts.goCueWindows.execution(1) & S.tb < opts.goCueWindows.execution(2);
    respVisual = squeeze(mean(S.MUA(:, isCorrect, visIdx), 3)) - baseTrials;
    respPlan   = squeeze(mean(S.MUA(:, isCorrect, planIdx), 3)) - baseTrials;
    respExec   = squeeze(mean(S.MUA(:, isCorrect, execIdx), 3)) - baseTrials;
    epochResp = {respVisual, respPlan, respExec};

    respChans = find(screening(s).responsive);

    rowIdx = numel(results) + 1;
    results(rowIdx).Day = S.Day; %#ok<AGROW>
    results(rowIdx).RunN = S.RunN;
    for e = 1:3
        resp = epochResp{e}(respChans, :);
        [cleaned, varExplained] = removeCommonMode(resp, dirGroup);

        accOrig = looDecodeAcc(resp, dirGroup);
        accClean = looDecodeAcc(cleaned, dirGroup);

        results(rowIdx).(sprintf('%s_orig', epochNames{e})) = accOrig;
        results(rowIdx).(sprintf('%s_cleaned', epochNames{e})) = accClean;
        results(rowIdx).(sprintf('%s_pc1var', epochNames{e})) = varExplained;

        fprintf('%-10s %-10s %7.1f%% %7.1f%% %7.1f%%\n', S.Day, epochNames{e}, 100*accOrig, 100*accClean, 100*varExplained);
    end
end

outFile = fullfile(saccadeDataDir(), 'decode_common_mode_removed.mat');
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
function [cleaned, varExplained] = removeCommonMode(resp, dirGroup)
% Direction-demean (so PC1 can't just be "real tuning"), extract PC1 of
% the residual, then subtract the PC1 reconstruction from the ORIGINAL
% (non-demeaned) response -- same method as checkSharedComponentKinematics.m.
resid = nan(size(resp));
for d = 1:8
    idx = dirGroup == d;
    resid(:, idx) = resp(:, idx) - mean(resp(:, idx), 2);
end
[U, Sig, V] = svd(resid, 'econ');
loading = U(:,1);
score = V(:,1) * Sig(1,1);
varExplained = Sig(1,1)^2 / sum(diag(Sig).^2);
cleaned = resp - loading * score';
end

% ------------------------------------------------------------------
function acc = looDecodeAcc(resp, dirGroup)
% Leave-one-trial-out nearest-centroid decoding accuracy -- identical
% method to decodeDirection.m's looDecodeAcc (duplicated locally rather
% than exported, matching this codebase's existing convention of small
% per-file helper duplication, e.g. bhFDR/redblueLocal).
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
