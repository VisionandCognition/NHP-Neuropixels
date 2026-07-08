% PATCHEYEBLINKANDCALIB One-off patch (2026-07-08) for the 12 existing
% saccade-direction _extracted.mat files, predating the discovery that
% blinks corrupt the raw eye X/Y position channels (not just pupil) --
% see extractSaccadeSession.m's eyeBlinkMask comment for the full story.
% Adds S.eyeBlinkMask and recomputes S.eyeCalib excluding blink samples
% from the landing-window fit. Does NOT touch MUA/AP data, so this is
% cheap (no need to redo the expensive raw-AP extraction).

fld.script = '/media/DOCUMENTS/DOCUMENTS/EPHYS_ANALYSIS/NHP-Neuropixels';
addpath(genpath(fullfile(fld.script,'OpenEphys')));
addpath(fullfile(fld.script,'SaccadeDirection','code'));

dataDir = saccadeDataDir();
files = dir(fullfile(dataDir, '*_extracted.mat'));
files = files(~contains({files.name}, 'summary'));

for i = 1:numel(files)
    fn = fullfile(files(i).folder, files(i).name);
    S = load(fn, 'Day', 'RunN', 'eyeRaw', 'tbEye', 'trialTable', 'Stm', 'eyeCalib');

    railV = -5.00; railTol = 0.5;
    eyeBlinkMask = squeeze(abs(S.eyeRaw(1,:,:) - railV) < railTol | abs(S.eyeRaw(2,:,:) - railV) < railTol); %#ok<NASGU>

    targXY = cat(1, S.Stm.Targ.xy_dva);
    landWin = S.tbEye > 0.20 & S.tbEye < 0.35;
    isCorrect = strcmp(S.trialTable.outcome, 'correct') & ~isnan(S.trialTable.targ_idx);
    trls = find(isCorrect);

    landX = nan(numel(trls),1); landY = nan(numel(trls),1);
    tgtX = nan(numel(trls),1); tgtY = nan(numel(trls),1);
    blinkMask = eyeBlinkMask;
    for k = 1:numel(trls)
        t = trls(k);
        blink = blinkMask(t, landWin);
        x = squeeze(S.eyeRaw(1,t,landWin)); x(blink) = nan;
        y = squeeze(S.eyeRaw(2,t,landWin)); y(blink) = nan;
        landX(k) = mean(x, 'omitnan');
        landY(k) = mean(y, 'omitnan');
        tgtX(k) = targXY(S.trialTable.targ_idx(t),1);
        tgtY(k) = targXY(S.trialTable.targ_idx(t),2);
    end
    ok = ~isnan(landX) & ~isnan(landY);

    eyeCalib = struct();
    if sum(ok) < 8
        eyeCalib.valid = false;
        eyeCalib.note = 'insufficient correct trials for calibration';
    else
        A = [landX(ok), landY(ok), ones(sum(ok),1)];
        Mx = A \ tgtX(ok);
        My = A \ tgtY(ok);
        predX = A*Mx; predY = A*My;
        r2X = 1 - sum((tgtX(ok)-predX).^2)/sum((tgtX(ok)-mean(tgtX(ok))).^2);
        r2Y = 1 - sum((tgtY(ok)-predY).^2)/sum((tgtY(ok)-mean(tgtY(ok))).^2);

        eyeCalib.valid = true;
        eyeCalib.M = [Mx, My];
        eyeCalib.r2X = r2X; eyeCalib.r2Y = r2Y;
        eyeCalib.nCalibTrials = sum(ok);
        eyeCalib.landWin = [0.20 0.35];
    end

    fprintf('%s run-%03d: old R2x=%.3f R2y=%.3f (n=%d) -> new R2x=%.3f R2y=%.3f (n=%d), %d/%d calib trials had a blink in the landing window\n', ...
        S.Day, S.RunN, S.eyeCalib.r2X, S.eyeCalib.r2Y, S.eyeCalib.nCalibTrials, ...
        eyeCalib.r2X, eyeCalib.r2Y, eyeCalib.nCalibTrials, ...
        sum(any(blinkMask(trls,landWin),2)), numel(trls));

    save(fn, 'eyeBlinkMask', 'eyeCalib', '-append');
end
fprintf('done\n');
