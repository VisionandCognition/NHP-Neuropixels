function [lineMap, lineCounts] = matchTTLtoLogEvents(lines, eventTypeCounts, verbose)
% MATCHTTLTOLOGEVENTS Match TTL digital line numbers to behavioral event
% types purely by rising-edge count, since the hardcoded bit numbers in
% par_default.m (StimB=1, TargetB=2, RewardB=3, CorrectB=7, ...) do not
% match what is actually wired on this rig's NI-DAQ digital port (verified
% empirically per-session instead of trusted from that file).
%
% [lineMap, lineCounts] = matchTTLtoLogEvents(lines, eventTypeCounts, verbose)
%   lines           : signed TTL line codes from states.npy (+n = rising
%                     edge on line n, -n = falling edge)
%   eventTypeCounts : containers.Map from event-type name -> count in Log.events
%   verbose         : print the matches (default true)
%
% Returns:
%   lineMap    : containers.Map from event-type name -> line number, only
%                for event types whose count uniquely matches exactly one
%                TTL line's rising-edge count.
%   lineCounts : containers.Map from line number (as string) -> rising count

if nargin < 3
    verbose = true;
end

posLines = lines(lines > 0);
uLines = unique(posLines);
counts = arrayfun(@(L) sum(posLines == L), uLines);

lineCounts = containers.Map('KeyType','char','ValueType','double');
for i = 1:numel(uLines)
    lineCounts(num2str(uLines(i))) = counts(i);
end

lineMap = containers.Map('KeyType','char','ValueType','double');
types = keys(eventTypeCounts);
for i = 1:numel(types)
    ty = types{i};
    n = eventTypeCounts(ty);
    matchIdx = find(counts == n);
    if numel(matchIdx) == 1
        lineMap(ty) = uLines(matchIdx);
        if verbose
            fprintf('  TTL match: line %d <-> ''%s'' (count=%d)\n', uLines(matchIdx), ty, n);
        end
    elseif verbose && numel(matchIdx) > 1
        fprintf('  TTL match: ''%s'' (count=%d) ambiguous with lines %s -- skipped\n', ...
            ty, n, mat2str(uLines(matchIdx)));
    end
end

end
