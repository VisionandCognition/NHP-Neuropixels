function anat = getChannelArea(day, runN, channelIdx)
% GETCHANNELAREA Looks up the tentative anatomical area for one or more
% channels of a given session, from the hand-filled
% anatomy/anatomy_template_channels.csv (see generateAnatomyTemplate.m
% for how that file was built, and PIPELINE_REPORT.md sec 4.12).
%
% These are TENTATIVE, hand-assigned labels compared against an atlas by
% eye (per-cluster, then propagated to every channel in that cluster) --
% not a verified histological registration. Treat accordingly.
%
% anat = getChannelArea(day, runN, channelIdx)
%   day        : 'YYYYMMDD' string (or numeric, matches either)
%   runN       : run number (numeric)
%   channelIdx : vector of channel indices (1 = deepest / probe tip,
%                matching every other per-channel array in this
%                pipeline, e.g. S.CHdepthUm). Default: every channel
%                present for that session in the anatomy table.
%
% Returns a scalar struct, each field a 1 x numel(channelIdx) array
% aligned to channelIdx:
%   .ClusterID     (double, NaN if no match)
%   .AssignedArea  (cellstr, '' if no match) -- the finest/tentative label
%   .L6 .L5 .L4 .L3 (cellstr, '' if no match) -- hierarchical parcellation,
%                    L6 finest to L3 coarsest (see anatomy_template_clusters_hierarchical.csv)
%
% Usage in this pipeline: loadAllSessions.m attaches this automatically
% as S.channelArea for every session, so any analysis can use e.g.
% S.channelArea.L4{c} (coarse-ish label) or S.channelArea.AssignedArea{c}
% (finest tentative label) for channel c without calling this directly.

persistent T
if isempty(T)
    anatFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'anatomy', 'anatomy_template_channels.csv');
    T = readtable(anatFile, 'TextType', 'string');
    T.Day = string(T.Day); % readtable may infer Day as numeric (all-digit strings) -- force string so lookups are robust regardless
end

dayStr = string(day);
rows = T.Day == dayStr & T.RunN == runN;
Tsub = T(rows, :);

if nargin < 3 || isempty(channelIdx)
    channelIdx = sort(Tsub.ChannelIdx)';
end
channelIdx = channelIdx(:)';
n = numel(channelIdx);

anat = struct();
anat.ClusterID = nan(1, n);
anat.AssignedArea = repmat({''}, 1, n);
anat.L6 = repmat({''}, 1, n);
anat.L5 = repmat({''}, 1, n);
anat.L4 = repmat({''}, 1, n);
anat.L3 = repmat({''}, 1, n);

if isempty(Tsub)
    warning('getChannelArea:noSession', 'No anatomy entries found for day %s run %d -- returning empty labels', dayStr, runN);
    return
end

[found, loc] = ismember(channelIdx, Tsub.ChannelIdx');
if any(~found)
    warning('getChannelArea:missingChannels', '%d/%d requested channels have no anatomy entry for day %s run %d', ...
        sum(~found), n, dayStr, runN);
end
anat.ClusterID(found) = Tsub.ClusterID(loc(found));
anat.AssignedArea(found) = cellstr(Tsub.AssignedArea(loc(found)));
anat.L6(found) = cellstr(Tsub.L6(loc(found)));
anat.L5(found) = cellstr(Tsub.L5(loc(found)));
anat.L4(found) = cellstr(Tsub.L4(loc(found)));
anat.L3(found) = cellstr(Tsub.L3(loc(found)));
end
