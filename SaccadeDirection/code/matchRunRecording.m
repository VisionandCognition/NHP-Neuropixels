function rec = matchRunRecording(dayPath, day, runTimeStr)
% MATCHRUNRECORDING Find the OpenEphys experiment/recording folder that
% corresponds to a given behavioral run, by matching the run's start time
% (encoded in its folder name, local Europe/Amsterdam time) against each
% recording's "Software Time" in sync_messages.txt (ms since epoch, UTC).
%
% The recording chosen is the one whose start time is CLOSEST in absolute
% value to the run's start time. Operator workflow is not consistent
% across sessions: usually OE recording is started a couple of seconds
% before the behavioral run script (matching folder timestamps to within
% ~2-10s), but on some days recording was started a few seconds AFTER the
% run folder was created instead (confirmed directly: 20260303 run-003's
% matching recording starts 6s *after* its folder timestamp). A match is
% only accepted within maxOffsetSec, so a day where the true recording
% hasn't started yet (e.g. 20260224 run-001, nearest recording is 72
% minutes away) correctly reports no match rather than a bogus one.
%
% rec = matchRunRecording(dayPath, day, runTimeStr)
%   dayPath    : full path to the day folder (.../m032/20260326)
%   day        : 'YYYYMMDD' string
%   runTimeStr : 'HHMMSS' string parsed from the run folder name
%
% Returns a struct with fields:
%   ExpN, RecN, RecordingPath, StructureOebinFile, SettingsFile,
%   APContinuousDir, APEventsTTLDir, NIDAQContinuousDir, NIDAQEventsTTLDir,
%   APInfo, NIDAQInfo (parsed continuous-stream info from structure.oebin)

maxOffsetSec = 90;

runDT = datetime([day runTimeStr], 'InputFormat', 'yyyyMMddHHmmss', ...
    'TimeZone', 'Europe/Amsterdam');

expDirs = dir(fullfile(dayPath, 'experiment*'));
expDirs = expDirs([expDirs.isdir]);

bestDT = NaT('TimeZone', 'Europe/Amsterdam');
bestAbsDiff = Inf;
best = struct('ExpN', [], 'RecN', [], 'RecordingPath', '');

for e = 1:numel(expDirs)
    expTok = regexp(expDirs(e).name, '^experiment(\d+)$', 'tokens', 'once');
    if isempty(expTok), continue; end
    expN = str2double(expTok{1});

    recDirs = dir(fullfile(dayPath, expDirs(e).name, 'recording*'));
    recDirs = recDirs([recDirs.isdir]);
    for rr = 1:numel(recDirs)
        recTok = regexp(recDirs(rr).name, '^recording(\d+)$', 'tokens', 'once');
        if isempty(recTok), continue; end
        recN = str2double(recTok{1});

        recPath = fullfile(dayPath, expDirs(e).name, recDirs(rr).name);
        syncFile = fullfile(recPath, 'sync_messages.txt');
        if ~isfile(syncFile), continue; end

        recDT = readRecordingStartTime(syncFile);
        if isnat(recDT), continue; end

        absDiff = abs(seconds(recDT - runDT));
        if absDiff < bestAbsDiff
            bestAbsDiff = absDiff;
            bestDT = recDT;
            best.ExpN = expN;
            best.RecN = recN;
            best.RecordingPath = recPath;
        end
    end
end

if isempty(best.RecordingPath) || bestAbsDiff > maxOffsetSec
    error('matchRunRecording:noMatch', ...
        'No recording starts within %ds of run time %s on %s (closest: %.0fs)', ...
        maxOffsetSec, runTimeStr, day, bestAbsDiff);
end

rec = best;
rec.RecordingStartDT = bestDT; % absolute start time of the matched recording
% (Europe/Amsterdam), needed to localize a run within a recording that
% spans multiple behavioral runs (see extractSaccadeSession.m)
rec.StructureOebinFile = fullfile(rec.RecordingPath, 'structure.oebin');
rec.SettingsFile = fullfile(dayPath, 'settings.xml');
if ~isfile(rec.SettingsFile)
    % settings.xml can also be experiment-specific in some sessions
    altSettings = fullfile(dayPath, sprintf('settings_%d.xml', rec.ExpN));
    if isfile(altSettings)
        rec.SettingsFile = altSettings;
    end
end

info = jsondecode(fileread(rec.StructureOebinFile));
cont = info.continuous;
apIdx = find(cellfun(@(x) contains(x,'ProbeA-AP'), {cont.folder_name}), 1);
niIdx = find(cellfun(@(x) contains(x,'NI-DAQmx'), {cont.folder_name}), 1);
if isempty(apIdx)
    error('matchRunRecording:noAP', 'No ProbeA-AP stream found in %s', rec.StructureOebinFile);
end
if isempty(niIdx)
    error('matchRunRecording:noNIDAQ', 'No NI-DAQmx stream found in %s', rec.StructureOebinFile);
end

rec.APInfo = cont(apIdx);
rec.NIDAQInfo = cont(niIdx);

rec.APContinuousDir = fullfile(rec.RecordingPath, 'continuous', strtrim(cont(apIdx).folder_name));
rec.NIDAQContinuousDir = fullfile(rec.RecordingPath, 'continuous', strtrim(cont(niIdx).folder_name));
rec.APEventsTTLDir = fullfile(rec.RecordingPath, 'events', strtrim(cont(apIdx).folder_name), 'TTL');
rec.NIDAQEventsTTLDir = fullfile(rec.RecordingPath, 'events', strtrim(cont(niIdx).folder_name), 'TTL');

end

% ---------------------------------------------------------------------
function dt = readRecordingStartTime(syncFile)
dt = NaT('TimeZone', 'Europe/Amsterdam');
txt = fileread(syncFile);
tok = regexp(txt, 'Software Time \(milliseconds since midnight Jan 1st 1970 UTC\):\s*(\d+)', ...
    'tokens', 'once');
if isempty(tok)
    return
end
epochMs = str2double(tok{1});
dt = datetime(epochMs/1000, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
dt.TimeZone = 'Europe/Amsterdam';
end
