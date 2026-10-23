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
% or 'planned', which of the above applied for that day), chamber (2 or 3;
% '' for unknown days).
%
% CHAMBER 3 (prefrontal, days 20260819-20260918; burr holes from
% Chamber3Holes.docx -- hole numbers REUSE chamber-2's numbers, so always
% key on chamber + hole). For these days the chamber-2 convention of
% "estimatedInsertionDepthMm = planned depth" does NOT hold: the elab
% 'probe final depth' minus 'signal onset depth' (in-tissue depth) is
% often far shallower than the Chamber3Holes brain depth (e.g. 20260918:
% ~4 mm in tissue vs 19 mm planned). So for chamber 3:
%   estimatedInsertionDepthMm = probe final depth - signal onset depth
%     (depthSource 'elab_final_minus_onset'; onsetUncertain=true where the
%     log has a '?'), or NaN (depthSource 'unknown') where onset is blank.
%   plannedDepthMm is still returned (the Chamber3Holes brain depth) for
%   comparison, plus probeFinalDepthMm / signalOnsetDepthMm raw fields.
% 20260819 was added to Kwibus_elab.pdf in the 2026-10-23 re-export (hole 1,
% '4 bank', onset 32?, final 46; guide-tube final depth blank).
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

persistent T T3 holeInfo holeInfo3
if isempty(T)
    % CHAMBER 2: day -> {burrHole, channelSelection}, from Kwibus_elab.pdf
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

    % CHAMBER 3: day -> burr hole, channel selection, elab probe final depth
    % and signal onset depth (mm; NaN = blank/"don't know" in the log;
    % onsetUncertain = log entry carries a '?'). From Kwibus_elab.pdf.
    T3 = containers.Map();
    T3('20260819') = c3('1', '4 bank',        46,   32,   true);
    T3('20260821') = c3('5', 'single column', 46,   37.5, false); % log: 37-38
    T3('20260825') = c3('4', 'single column', 46.5, 38,   true);
    T3('20260903') = c3('2', 'single column', 38,   31,   true);
    T3('20260904') = c3('2', '4bankeven',     46,   nan,  false);
    T3('20260909') = c3('1', 'single column', 48,   32,   true);
    T3('20260910') = c3('4', 'single column', 43,   35,   false);
    T3('20260915') = c3('2', 'Bank A',        35.5, nan,  false); % final "35(.5?)"
    T3('20260916') = c3('4', 'single column', 44,   nan,  false); % no signal at any depth
    T3('20260918') = c3('2', 'Bank A',        38.5, 34.5, true);

    % burr hole -> planned trajectory, from Chamber2Holes.docx (superficial to deep)
    holeInfo = containers.Map();
    holeInfo('0') = struct('plannedDepthMm', 27, 'structures', {{'S1 BA1/2, BA3', 'VLP, VPI, ZI'}});
    holeInfo('1') = struct('plannedDepthMm', 28, 'structures', {{'S1 BA1/2, PCC BA23', 'SC'}});
    holeInfo('2') = struct('plannedDepthMm', 27, 'structures', {{'S1 BA1/2 or PE BA5', 'VIP', 'Inferior/Anterior pulvinar (Medial pulvinar, Suprageniculate, MGN)'}});
    holeInfo('4') = struct('plannedDepthMm', 25, 'structures', {{'S1 BA1/2, BA3, (M1), PCC BA23', 'Thalamus MD, CL, CM'}});
    holeInfo('5') = struct('plannedDepthMm', 25, 'structures', {{'S1 BA1/2, BA3, PCC BA23, retrosplenial BA29', 'Habenula, Posterior commissure or CM'}});

    % burr hole -> planned trajectory, from Chamber3Holes.docx ("Brain N mm"
    % = planned depth within brain). Cortical-only trajectories.
    holeInfo3 = containers.Map();
    holeInfo3('1') = struct('plannedDepthMm', 16,   'structures', {{'OFC, VLPFC, DLPFC (lateral BA13, ventral BA46, dorsal BA46)'}});
    holeInfo3('2') = struct('plannedDepthMm', 19,   'structures', {{'OFC, VLPFC, DLPFC, FEF dorsal (medial BA13, BA46, dorsal 8A)'}});
    holeInfo3('4') = struct('plannedDepthMm', 13,   'structures', {{'OFC, VLPFC, DLPFC (medial BA12, ventral BA46, dorsal BA46)'}});
    holeInfo3('5') = struct('plannedDepthMm', 10.8, 'structures', {{'?, VLPFC (medial BA12, (BA45a), ventral BA46)'}});
end

if isKey(T3, day)
    e = T3(day);
    info = struct();
    info.chamber = 3;
    info.burrHole = e.burrHole;
    info.channelSelection = e.channelSelection;
    if isKey(holeInfo3, e.burrHole)
        hi = holeInfo3(e.burrHole);
        info.plannedDepthMm = hi.plannedDepthMm;
        info.structures = hi.structures;
    else
        info.plannedDepthMm = nan;
        info.structures = {};
    end
    info.probeFinalDepthMm = e.probeFinalDepthMm;
    info.signalOnsetDepthMm = e.signalOnsetDepthMm;
    info.onsetUncertain = e.onsetUncertain;
    if ~isnan(e.probeFinalDepthMm) && ~isnan(e.signalOnsetDepthMm)
        info.estimatedInsertionDepthMm = e.probeFinalDepthMm - e.signalOnsetDepthMm;
        info.depthSource = 'elab_final_minus_onset';
    else
        info.estimatedInsertionDepthMm = nan;
        info.depthSource = 'unknown';
    end
    return
end

if ~isKey(T, day)
    info = struct('chamber', '', 'burrHole', '?', 'channelSelection', 'unknown', ...
        'plannedDepthMm', nan, 'structures', {{}}, ...
        'estimatedInsertionDepthMm', nan, 'depthSource', 'unknown');
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
info.chamber = 2;
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

function e = c3(hole, sel, finalDepth, onset, uncertain)
e = struct('burrHole', hole, 'channelSelection', sel, 'probeFinalDepthMm', finalDepth, ...
    'signalOnsetDepthMm', onset, 'onsetUncertain', uncertain);
end
