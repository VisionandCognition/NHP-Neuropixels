% RUNALLEXTRACTION Batch-extract MUA + eye data for every identified
% 8-target saccade-direction run of m032, saving one .mat file per run
% under saccadeDataDir() (on the server, not the local repo). Sessions where no OpenEphys recording
% covers the behavioral run (clock mismatch / recording not yet started)
% are skipped and logged.
%
% Set skipExisting = true in the workspace before calling to only extract
% runs that don't have an *_extracted.mat yet (e.g. after new recordings
% were added); default false re-extracts everything.

fld.script = '/media/DOCUMENTS/DOCUMENTS/EPHYS_ANALYSIS/NHP-Neuropixels';
addpath(genpath(fullfile(fld.script,'OpenEphys')));
addpath(fullfile(fld.script,'SaccadeDirection','code'));

runs = findSaccadeRuns();
if ~exist('skipExisting','var'), skipExisting = false; end

results = struct('Day',{},'RunN',{},'status',{},'message',{});
for i = 1:numel(runs)
    r = runs(i);
    outFile = fullfile(saccadeDataDir(), sprintf('%s_run-%03d_extracted.mat', r.Day, r.RunN));
    if skipExisting && isfile(outFile)
        results(i).Day = r.Day;
        results(i).RunN = r.RunN;
        results(i).status = 'ok';
        results(i).message = 'already extracted (skipped)';
        continue
    end
    try
        tic;
        S = extractSaccadeSession(r, struct('save', true));
        et = toc;
        results(i).Day = r.Day;
        results(i).RunN = r.RunN;
        results(i).status = 'ok';
        results(i).message = sprintf('%d trials, %.1fs', height(S.trialTable), et);
    catch ME
        results(i).Day = r.Day;
        results(i).RunN = r.RunN;
        results(i).status = 'FAILED';
        results(i).message = ME.message;
        fprintf(2, 'FAILED %s run-%03d: %s\n', r.Day, r.RunN, ME.message);
    end
end

fprintf('\n=== Extraction summary ===\n');
for i = 1:numel(results)
    fprintf('%s run-%03d: %s -- %s\n', results(i).Day, results(i).RunN, ...
        results(i).status, results(i).message);
end

save(fullfile(saccadeDataDir(), 'extraction_summary.mat'), 'results');
