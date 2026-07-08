% RUNFULLPIPELINE End-to-end saccade-direction analysis for m032:
%   1. extract MUA + eye data for every 8-target saccade run (runAllExtraction)
%   2. screen every channel for task responsiveness (screenResponsiveChannels)
%   3. compute saccade-direction tuning for responsive channels (computeDirectionTuning)
%   4. calibrate eye data and characterize saccade dynamics / fixation / pupil (analyzeEyeData)
%
% Run this after code has been reviewed; step 1 takes several minutes per
% session (raw AP filtering across 384 channels). All outputs land in
% saccadeDataDir() (per-session *_extracted.mat, screening_results.mat,
% tuning_results.mat, eye_summary.mat) and <saccadeDataDir()>/figures/ --
% moved off the local repo onto the server; see saccadeDataDir.m.

fld.script = '/media/DOCUMENTS/DOCUMENTS/EPHYS_ANALYSIS/NHP-Neuropixels';
addpath(genpath(fullfile(fld.script,'OpenEphys')));
addpath(fullfile(fld.script,'SaccadeDirection','code'));

runAllExtraction;

sessions = loadAllSessions();
screening = screenResponsiveChannels(sessions);
tuning = computeDirectionTuning(sessions, screening);
eyeSummary = analyzeEyeData(sessions);
artifactCheck = checkArtifactCorrelation(sessions);
sharedCheck = checkSharedComponentKinematics(sessions, eyeSummary);
planExecResults = analyzePlanningExecution(sessions, eyeSummary);
planExecSummary = summarizePlanningExecution('saclocked'); %#ok<NASGU>

fprintf('\n=== Pipeline complete ===\n');
fprintf('%d sessions, %d channel-session rows tuned (of %d screened responsive)\n', ...
    numel(sessions), sum([tuning.tuned]), numel(tuning));
