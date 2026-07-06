% RUNALLEXTRACTION Batch-extract MUA + eye data for every identified
% 8-target saccade-direction run of m032, saving one .mat file per run
% under saccadeDataDir() (on the server, not the local repo). Sessions where no OpenEphys recording
% covers the behavioral run (clock mismatch / recording not yet started)
% are skipped and logged.

fld.script = '/media/DOCUMENTS/DOCUMENTS/EPHYS_ANALYSIS/NHP-Neuropixels';
addpath(genpath(fullfile(fld.script,'OpenEphys')));
addpath(fullfile(fld.script,'SaccadeDirection','code'));

runs = findSaccadeRuns();

results = struct('Day',{},'RunN',{},'status',{},'message',{});
for i = 1:numel(runs)
    r = runs(i);
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
