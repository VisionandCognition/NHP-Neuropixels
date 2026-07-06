function runs = findSaccadeRuns(dataRoot, monkey)
% FINDSACCADERUNS Locate all 8-target saccade-direction runs for a subject.
%
% A run is identified as a saccade-direction run by the presence of
% stim_ct_classic_8targ.m in its run-*/ folder (confirmed convention).
%
% runs = findSaccadeRuns(dataRoot, monkey)
%   dataRoot : path to Data_collection folder
%   monkey   : subject id, e.g. 'm032'
%
% Returns a struct array with fields:
%   Day, RunFolder, RunPath, RunN, RunTimeStr, LogFile, SessionJsonFile

if nargin < 1 || isempty(dataRoot)
    dataRoot = '/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_collection';
end
if nargin < 2 || isempty(monkey)
    monkey = 'm032';
end

subjDir = fullfile(dataRoot, monkey);
dayDirs = dir(subjDir);
dayDirs = dayDirs([dayDirs.isdir] & ~startsWith({dayDirs.name},'.'));

runs = struct('Day',{},'RunFolder',{},'RunPath',{},'RunN',{}, ...
    'RunTimeStr',{},'LogFile',{},'SessionJsonFile',{});

for d = 1:numel(dayDirs)
    day = dayDirs(d).name;
    dayPath = fullfile(subjDir, day);
    runDirs = dir(fullfile(dayPath, 'run-*'));
    runDirs = runDirs([runDirs.isdir]);
    for r = 1:numel(runDirs)
        runFolder = runDirs(r).name;
        runPath = fullfile(dayPath, runFolder);
        marker = dir(fullfile(runPath, 'stim_ct_classic_8targ.m'));
        if isempty(marker)
            continue
        end

        tok = regexp(runFolder, '^run-(\d+)_T-(\d+)$', 'tokens', 'once');
        if isempty(tok)
            warning('findSaccadeRuns:badRunFolder', ...
                'Could not parse run folder name %s', runFolder);
            continue
        end
        runN = str2double(tok{1});
        runTimeStr = tok{2};

        logFile = dir(fullfile(runPath, 'sub-*.mat'));
        sessJson = dir(fullfile(runPath, 'sub-*_session.json'));
        if isempty(logFile)
            warning('findSaccadeRuns:noLogFile', ...
                'No log .mat file found in %s', runPath);
            continue
        end

        i = numel(runs) + 1;
        runs(i).Day = day;
        runs(i).RunFolder = runFolder;
        runs(i).RunPath = runPath;
        runs(i).RunN = runN;
        runs(i).RunTimeStr = runTimeStr;
        runs(i).LogFile = fullfile(logFile(1).folder, logFile(1).name);
        if ~isempty(sessJson)
            runs(i).SessionJsonFile = fullfile(sessJson(1).folder, sessJson(1).name);
        else
            runs(i).SessionJsonFile = '';
        end
    end
end

fprintf('findSaccadeRuns: found %d saccade-direction (8targ) runs for %s\n', ...
    numel(runs), monkey);

end
