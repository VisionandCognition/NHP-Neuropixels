function info = sessionAnatomyInfo(day)
% SESSIONANATOMYINFO Burr hole, online channel-selection scheme, and
% planned trajectory (target structures, superficial to deep) for a given
% m032 recording day.
%
% info = sessionAnatomyInfo(day)
%   day : 'YYYYMMDD' string
%
% Returns a struct with fields: burrHole (string, may end '?' if the
% source log itself was uncertain), channelSelection (string, the online
% electrode-bank/section selection used that day), plannedDepthMm
% (approximate planned total insertion depth for that burr hole, from
% Chamber2Holes.docx), structures (cell array, superficial to deep),
% estimatedInsertionDepthMm (best available estimate of how far the probe
% tip actually went into the brain that day -- an explicit per-session
% value where the elab log encodes one directly in the channel-selection
% string (e.g. "Ch2H2@53mm" -> 53mm; only 20260324/20260325 have this),
% otherwise the burr hole's Chamber2Holes PLANNED depth, since actual
% insertion depth tracks the plan closely in the sessions where both are
% known -- see PIPELINE_REPORT.md sec 4.12/4.13), depthSource ('explicit'
% or 'planned', which of the above applied for that day).
%
% SOURCE AND LIMITATIONS (read before using this for anything beyond a
% qualitative "roughly what region" annotation):
%   - Day -> burr hole and channel-selection strings are transcribed
%     directly from Kwibus_elab.pdf (the lab notebook), one entry per
%     recording day. 20260224 is marked with '?' in the source log itself
%     (burr hole 5 assumed but not certain) -- reflected here.
%   - structures/plannedDepthMm are the PLANNED trajectory per burr hole
%     from Chamber2Holes.docx, i.e. what the MRI-based planning predicted
%     that hole should pass through at full planned depth -- not a
%     verified, session-specific, per-channel anatomical registration.
%   - There is NOT enough clean, consistently-recorded numeric data in
%     Kwibus_elab.pdf to convert this into an absolute "mm below cortical
%     surface" per channel for most sessions: the log has fields for
%     "guide tube final depth", "probe tip enter guide tube depth", and
%     "signal onset depth" (the latter being the actual dura/brain-surface
%     crossing reference needed for this), but "signal onset depth" is
%     blank/unfilled for most of these 12 days, and only 2 of 12
%     (20260324, 20260325) have an explicit final recording depth encoded
%     in the channel-selection string itself (e.g. "Ch2H2@53mm"). Treat
%     this function's output as "which structures this burr hole's
%     trajectory generally passes through," not as a per-channel lookup.

persistent T
if isempty(T)
    % day -> {burrHole, channelSelection}, from Kwibus_elab.pdf
    T = containers.Map();
    T('20260224') = struct('burrHole', '5?', 'channelSelection', 'Ch2H5quarter');
    T('20260226') = struct('burrHole', '4',  'channelSelection', 'Ch2H4quarter');
    T('20260303') = struct('burrHole', '5',  'channelSelection', 'Ch2h5quarterD');
    T('20260310') = struct('burrHole', '5',  'channelSelection', '4bankeven');
    T('20260312') = struct('burrHole', '1',  'channelSelection', 'Single column');
    T('20260317') = struct('burrHole', '2',  'channelSelection', 'Ch2H2sections3quarterD');
    T('20260318') = struct('burrHole', '0',  'channelSelection', 'Ch2H0sections2');
    T('20260320') = struct('burrHole', '4',  'channelSelection', 'Ch2H4quarterD');
    T('20260323') = struct('burrHole', '5',  'channelSelection', 'Single column');
    T('20260324') = struct('burrHole', '2',  'channelSelection', 'Ch2H2@53mm');
    T('20260325') = struct('burrHole', '1',  'channelSelection', 'Ch2H1@53');
    T('20260326') = struct('burrHole', '2',  'channelSelection', '4bankeven');
end

% burr hole -> planned trajectory, from Chamber2Holes.docx (superficial to deep)
holeInfo = containers.Map();
holeInfo('0') = struct('plannedDepthMm', 27, 'structures', {{'S1 BA1/2, BA3', 'VLP, VPI, ZI'}});
holeInfo('1') = struct('plannedDepthMm', 28, 'structures', {{'S1 BA1/2, PCC BA23', 'SC'}});
holeInfo('2') = struct('plannedDepthMm', 27, 'structures', {{'S1 BA1/2 or PE BA5', 'VIP', 'Inferior/Anterior pulvinar (Medial pulvinar, Suprageniculate, MGN)'}});
holeInfo('4') = struct('plannedDepthMm', 25, 'structures', {{'S1 BA1/2, BA3, (M1), PCC BA23', 'Thalamus MD, CL, CM'}});
holeInfo('5') = struct('plannedDepthMm', 25, 'structures', {{'S1 BA1/2, BA3, PCC BA23, retrosplenial BA29', 'Habenula, Posterior commissure or CM'}});

if ~isKey(T, day)
    info = struct('burrHole', '?', 'channelSelection', 'unknown', ...
        'plannedDepthMm', nan, 'structures', {{}});
    return
end

entry = T(day);
holeKey = regexprep(entry.burrHole, '\?', ''); % strip uncertainty marker for the lookup
hi = holeInfo(holeKey);

% explicit per-session final depth, where the elab channel-selection
% string encodes one directly (e.g. "Ch2H2@53mm", "Ch2H1@53") -- more
% precise than the generic per-hole plan where available.
explicitDepthMm = containers.Map();
explicitDepthMm('20260324') = 53;
explicitDepthMm('20260325') = 53;

info = struct();
info.burrHole = entry.burrHole;
info.channelSelection = entry.channelSelection;
info.plannedDepthMm = hi.plannedDepthMm;
info.structures = hi.structures;
if isKey(explicitDepthMm, day)
    info.estimatedInsertionDepthMm = explicitDepthMm(day);
    info.depthSource = 'explicit';
else
    info.estimatedInsertionDepthMm = hi.plannedDepthMm;
    info.depthSource = 'planned';
end

end
