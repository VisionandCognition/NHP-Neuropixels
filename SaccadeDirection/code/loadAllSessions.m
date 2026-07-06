function sessions = loadAllSessions(dataDir)
% LOADALLSESSIONS Load every *_extracted.mat produced by
% extractSaccadeSession.m / runAllExtraction.m into a struct array.

if nargin < 1 || isempty(dataDir)
    dataDir = saccadeDataDir();
end

files = dir(fullfile(dataDir, '*_run-*_extracted.mat'));
sessions = struct([]);
for i = 1:numel(files)
    fp = fullfile(files(i).folder, files(i).name);
    fprintf('loading %s\n', files(i).name);
    S = load(fp);
    if isempty(sessions)
        sessions = S;
    else
        sessions(end+1) = S; %#ok<AGROW>
    end
end
fprintf('loaded %d sessions\n', numel(sessions));

end
