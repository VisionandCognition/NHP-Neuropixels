% PATCHEXISTINGEXTRACTIONS One-off patch for the 11 saccade-direction
% _extracted.mat files generated before extractSaccadeSession.m stored a
% verified S.LogFile field. Adds S.LogFile (the real behavioral log path,
% from findSaccadeRuns) to each existing file in-place, since the
% S.Par.LogFn/S.Par.LogPath already inside those files were found to be
% stale (see extractSaccadeSession.m comment). Also rebuilds
% extraction_summary.mat from the actual files on disk, since the
% original summary was stale and wrongly marked 20260317/20260318 as
% FAILED even though both have valid extracted data.

fld.script = '/media/DOCUMENTS/DOCUMENTS/EPHYS_ANALYSIS/NHP-Neuropixels';
addpath(genpath(fullfile(fld.script,'OpenEphys')));
addpath(fullfile(fld.script,'SaccadeDirection','code'));

runs = findSaccadeRuns();
dataDir = saccadeDataDir();

results = struct('Day',{},'RunN',{},'status',{},'message',{});
for i = 1:numel(runs)
    r = runs(i);
    results(i).Day = r.Day;
    results(i).RunN = r.RunN;
    outFile = fullfile(dataDir, sprintf('%s_run-%03d_extracted.mat', r.Day, r.RunN));
    if isfile(outFile)
        vars = who('-file', outFile);
        if ~ismember('LogFile', vars)
            S = load(outFile, 'trialTable');
            LogFile = r.LogFile;
            save(outFile, 'LogFile', '-append');
            fprintf('patched %s: added LogFile = %s\n', outFile, LogFile);
        else
            S = load(outFile, 'trialTable');
            fprintf('%s already has LogFile, skipped patch\n', outFile);
        end
        results(i).status = 'ok';
        results(i).message = sprintf('%d trials (existing extraction)', height(S.trialTable));
    else
        try
            dayPath = fileparts(r.RunPath);
            rec = matchRunRecording(dayPath, r.Day, r.RunTimeStr); %#ok<NASGU>
            results(i).status = 'FAILED';
            results(i).message = 'recording matched but no extracted output file present on disk';
        catch ME
            results(i).status = 'FAILED';
            results(i).message = ME.message;
        end
    end
end

fprintf('\n=== Rebuilt extraction summary ===\n');
for i = 1:numel(results)
    fprintf('%s run-%03d: %s -- %s\n', results(i).Day, results(i).RunN, ...
        results(i).status, results(i).message);
end

save(fullfile(dataDir, 'extraction_summary.mat'), 'results');
fprintf('\nSaved corrected extraction_summary.mat\n');
