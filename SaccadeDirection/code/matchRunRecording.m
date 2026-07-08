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
% genuinely hasn't started yet correctly reports no match rather than a
% bogus one.
%
% Recording folders are searched RECURSIVELY under dayPath (not just
% dayPath/experiment*/recording*), because at least one day (20260224)
% stores its real recordings under a separate, differently-named
% "Record Node" folder sitting alongside dayPath/experiment1 rather than
% inside it (dayPath/m032-2026-02-24_12-59-28/Record Node 110/experiment1/
% recording4). settings.xml is looked up next to wherever the matched
% experimentN folder actually lives, not assumed to be at dayPath root.
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

% Recursive search for recording* folders under dayPath. Not all days use
% the flat dayPath/experimentN/recordingM layout -- e.g. 20260224 has its
% real matching recording under dayPath/m032-2026-02-24_12-59-28/Record
% Node 110/experiment1/recording4, a separate/older-style OpenEphys
% "Record Node" folder sitting *alongside* dayPath/experiment1 rather than
% inside it. A non-recursive dir(fullfile(dayPath,'experiment*')) misses
% this entirely and silently searches only the (wrong/irrelevant) flat
% experiment1 folder, which is how 20260224/run-001 was originally
% mis-reported as having no matching recording at all: the true match was
% just never looked at.
recDirs = dir(fullfile(dayPath, '**', 'recording*'));
recDirs = recDirs([recDirs.isdir]);

bestDT = NaT('TimeZone', 'Europe/Amsterdam');
bestAbsDiff = Inf;
best = struct('ExpN', [], 'RecN', [], 'RecordingPath', '', 'BaseDir', '');

for rr = 1:numel(recDirs)
    recTok = regexp(recDirs(rr).name, '^recording(\d+)$', 'tokens', 'once');
    if isempty(recTok), continue; end
    recN = str2double(recTok{1});

    [expDir, expName] = fileparts(recDirs(rr).folder);
    expTok = regexp(expName, '^experiment(\d+)$', 'tokens', 'once');
    if isempty(expTok), continue; end
    expN = str2double(expTok{1});

    recPath = fullfile(recDirs(rr).folder, recDirs(rr).name);
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
        best.BaseDir = expDir; % folder containing experimentN, where its settings.xml lives
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
rec.SettingsFile = fullfile(rec.BaseDir, 'settings.xml');
if ~isfile(rec.SettingsFile)
    % settings.xml can also be experiment-specific in some sessions
    altSettings = fullfile(rec.BaseDir, sprintf('settings_%d.xml', rec.ExpN));
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
