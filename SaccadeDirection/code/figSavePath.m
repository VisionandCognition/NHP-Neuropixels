function fp = figSavePath(figDir, typeName, day, runN)
% FIGSAVEPATH Returns <figDir>/<typeName>/<day>_run-<runN>.png, creating
% the <typeName> subfolder if needed. Centralizes the "organize by figure
% TYPE, not by session" convention used by every per-session plotting
% function in this pipeline -- figures/<type>/<day>_run-<runN>.png -- so
% the same plot across all 12 sessions is easy to find and compare side
% by side. Cross-session summary figures (no single day/run) go in
% figures/pooled/ instead, and the per-channel direction timecourses keep
% their own separate session-first hierarchy (figures/timecourses/, see
% plotDirectionTimecourses.m -- too many files for either convention to
% help, organized differently on purpose).
%
% fp = figSavePath(figDir, typeName, day, runN)

typeDir = fullfile(figDir, typeName);
if ~isfolder(typeDir), mkdir(typeDir); end
fp = fullfile(typeDir, sprintf('%s_run-%03d.png', day, runN));
end
