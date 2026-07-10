function sessions = loadAllSessions(dataDir)
% LOADALLSESSIONS Load every *_extracted.mat produced by
% extractSaccadeSession.m / runAllExtraction.m into a struct array.
%
% Also attaches S.channelArea (see getChannelArea.m) -- tentative,
% hand-assigned anatomical area per channel, from
% anatomy/anatomy_template_channels.csv (PIPELINE_REPORT.md sec 4.12).
% S.channelArea.AssignedArea{c} is the finest tentative label for channel
% c; S.channelArea.L6/.L5/.L4/.L3{c} are the hierarchical parcellation
% levels (L6 finest, L3 coarsest); S.channelArea.DepthBelowSurfaceMm(c) is
% an ESTIMATED absolute depth below the brain surface (§4.13), anchored
% from S.CHdepthUm and that day's estimated probe insertion depth.

if nargin < 1 || isempty(dataDir)
    dataDir = saccadeDataDir();
end

files = dir(fullfile(dataDir, '*_run-*_extracted.mat'));
sessions = struct([]);
for i = 1:numel(files)
    fp = fullfile(files(i).folder, files(i).name);
    fprintf('loading %s\n', files(i).name);
    S = load(fp);
    nChan = size(S.MUA, 1);
    chDepthUm = []; if isfield(S, 'CHdepthUm'), chDepthUm = S.CHdepthUm; end
    S.channelArea = getChannelArea(S.Day, S.RunN, 1:nChan, chDepthUm);
    if isempty(sessions)
        sessions = S;
    else
        sessions(end+1) = S; %#ok<AGROW>
    end
end
fprintf('loaded %d sessions\n', numel(sessions));

end
