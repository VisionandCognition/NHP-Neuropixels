function U = loadSortedUnits(monkey, day, expN, varargin)
% LOADSORTEDUNITS Load sorted units for one (monkey, day, experiment).
%
% U = loadSortedUnits('m032', '20260312', 1)
% U = loadSortedUnits(..., 'goodOnly', true)          keep only units flagged good
% U = loadSortedUnits(..., 'node', 'm032-.../Record Node 110')   non-flat layouts (e.g. 20260224)
%
% Spike times are in seconds on the OpenEphys timeline -- the same clock as the
% NI-DAQ TTL timestamps used in extractSaccadeSession.m, so trial event times
% (ttl.timestamps) can be compared directly with U.times{k}.
%
% Fields:
%   U.times     {nUnits x 1} spike times (s)
%   U.units     table (unit_id, firing_rate, snr, ..., x_um, y_um, peak_channel_y_um, good, ...)
%   U.segments  struct: one entry per concatenated OE recording (recording, n_samples,
%               first_timestamp_s, ...). 
%   U.fs        AP sampling rate
%   U.dir       folder loaded from
p = inputParser;
addParameter(p, 'goodOnly', false);
addParameter(p, 'node', '.');
parse(p, varargin{:});

if strcmp(p.Results.node, '.')
    name = sprintf('exp%d', expN);
else
    slug = strrep(strrep(strrep(p.Results.node, '\', '/'), '/', '__'), ' ', '_');
    name = sprintf('%s__exp%d', slug, expN);
end
d = fullfile(spikesortDataDir(), monkey, day, name);
if ~isfile(fullfile(d, 'DONE.json'))
    error('loadSortedUnits:notDone', 'No finished spike sort in %s', d);
end

S = load(fullfile(d, 'spikes.mat'));
T = readtable(fullfile(d, 'units.csv'));

U.times = S.times;
U.units = T;
U.segments = S.segments;
U.fs = S.fs;
U.dir = d;

if p.Results.goodOnly
    keep = T.good == 1;
    U.times = U.times(keep);
    U.units = T(keep, :);
end
end
