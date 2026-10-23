function anat = getChannelArea(day, runN, channelIdx, chDepthUm)
% GETCHANNELAREA Looks up the tentative anatomical area for one or more
% channels of a given session, from the hand-filled
% anatomy/anatomy_template_channels.csv (see generateAnatomyTemplate.m
% for how that file was built, and PIPELINE_REPORT.md sec 4.12), plus any
% anatomy/chamber<N>/anatomy_template_channels.csv (later chambers, same
% label columns; stacked automatically).
%
% These are TENTATIVE, hand-assigned labels compared against an atlas by
% eye (per-cluster, then propagated to every channel in that cluster) --
% not a verified histological registration. Treat accordingly.
%
% If chDepthUm is supplied (S.CHdepthUm, i.e. each channel's real
% distance from the probe TIP in um), this also anchors an ESTIMATED
% absolute depth below the cortical/brain surface per channel, using
% sessionAnatomyInfo.m's estimatedInsertionDepthMm for that day (§4.13):
%   DepthBelowSurfaceMm(c) = estimatedInsertionDepthMm - chDepthUm(c)/1000
% Channel 1 (the tip) went deepest, so it gets the full insertion depth;
% a NEGATIVE value means that channel is estimated to sit ABOVE the
% brain surface (in the guide tube / dura / CSF, not brain tissue) at the
% depth the probe was actually driven to that day.
%
% anat = getChannelArea(day, runN, channelIdx, chDepthUm)
%   day        : 'YYYYMMDD' string (or numeric, matches either)
%   runN       : run number (numeric)
%   channelIdx : vector of channel indices (1 = deepest / probe tip,
%                matching every other per-channel array in this
%                pipeline, e.g. S.CHdepthUm). Default: every channel
%                present for that session in the anatomy table.
%   chDepthUm  : optional, S.CHdepthUm for this session (same length/
%                order as the full channel set, indexed by channelIdx) --
%                if omitted, DepthBelowSurfaceMm/EstimatedInsertionDepthMm/
%                DepthSource are returned as NaN/''.
%
% Returns a scalar struct, each field a 1 x numel(channelIdx) array
% aligned to channelIdx:
%   .ClusterID     (double, NaN if no match)
%   .AssignedArea  (cellstr, '' if no match) -- the finest/tentative label
%   .L6 .L5 .L4 .L3 (cellstr, '' if no match) -- hierarchical parcellation,
%                    L6 finest to L3 coarsest (see anatomy_template_clusters_hierarchical.csv)
%   .DepthBelowSurfaceMm (double, NaN unless chDepthUm supplied)
%   .EstimatedInsertionDepthMm (double, same value repeated -- for reference)
%   .DepthSource (cellstr, 'explicit' or 'planned', see sessionAnatomyInfo.m)
%
% Usage in this pipeline: loadAllSessions.m attaches this automatically
% as S.channelArea for every session, so any analysis can use e.g.
% S.channelArea.L4{c} (coarse-ish label) or S.channelArea.AssignedArea{c}
% (finest tentative label) for channel c without calling this directly.

persistent T
if isempty(T)
    % Chamber-2 file (original, hand-edited) plus any per-chamber files in
    % anatomy/chamber<N>/ -- same label columns, so just stack them. Only
    % the columns used below are kept (chamber files carry extra ones).
    anatDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'anatomy');
    anatFiles = {fullfile(anatDir, 'anatomy_template_channels.csv')};
    sub = dir(fullfile(anatDir, 'chamber*', 'anatomy_template_channels.csv'));
    for i = 1:numel(sub), anatFiles{end+1} = fullfile(sub(i).folder, sub(i).name); end %#ok<AGROW>
    keep = {'Day','RunN','ChannelIdx','ClusterID','AssignedArea','L6','L5','L4','L3'};
    parts = cell(size(anatFiles));
    for i = 1:numel(anatFiles)
        Ti = readtable(anatFiles{i}, 'TextType', 'string');
        Ti.Day = string(Ti.Day); % readtable may infer Day as numeric (all-digit strings) -- force string so lookups are robust regardless
        parts{i} = Ti(:, keep);
    end
    T = vertcat(parts{:});
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

anat.DepthBelowSurfaceMm = nan(1, n);
anat.EstimatedInsertionDepthMm = nan(1, n);
anat.DepthSource = repmat({''}, 1, n);
if nargin >= 4 && ~isempty(chDepthUm)
    info = sessionAnatomyInfo(char(dayStr));
    chDepthUm = chDepthUm(:)';
    valid = channelIdx >= 1 & channelIdx <= numel(chDepthUm);
    anat.DepthBelowSurfaceMm(valid) = info.estimatedInsertionDepthMm - chDepthUm(channelIdx(valid))/1000;
    anat.EstimatedInsertionDepthMm(valid) = info.estimatedInsertionDepthMm;
    anat.DepthSource(valid) = {info.depthSource};
end
end
