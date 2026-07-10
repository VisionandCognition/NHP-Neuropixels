function [tokens, names] = extractAreaTokens(labelStr)
% EXTRACTAREATOKENS Pulls out canonical area abbreviation tokens (the
% parenthetical short codes, e.g. "(SC)", "(Pul)", "(MThal)") from a
% free-text anatomy label string, for LIBERAL area pooling: many
% AssignedArea/L6-L3 labels in this pipeline are hand-typed with genuine
% uncertainty between candidates (e.g. "Posterior Commisure (or Habenula
% or CM?)", "SC or mesencephalic reticular formation"), so a channel
% should count toward EVERY area it might plausibly be, not just its
% single primary label.
%
% Using the parenthetical abbreviation as the canonical key (rather than
% the full free-text name before it) sidesteps inconsistent capitalization
% across entries for the same area (e.g. "(MThal)" vs "(Mthal)",
% "(GThal)" vs "(Gthal)", "(PCgG)" vs "(PcgG)" all appear in this dataset
% for what is clearly the same area) -- tokens are upper-cased so these
% collapse together automatically.
%
% [tokens, names] = extractAreaTokens(labelStr)
%   labelStr : a label string, or several concatenated together (e.g.
%              AssignedArea + L6 + L5 + L4 + L3 for one channel, so a
%              token appearing at ANY hierarchy level still counts)
%
% Returns:
%   tokens : cellstr of unique uppercase tokens found, e.g. {'SC','MRT'}.
%            Empty {} if none found (e.g. an empty or purely-qualifier
%            label). Extraction here is intentionally permissive and
%            UNCHANGED from before names was added, so pooling membership
%            never regresses based on the (best-effort, secondary) name
%            capture below.
%   names  : cellstr, same length as tokens -- the descriptive phrase
%            found immediately before that token's parenthetical in THIS
%            string (e.g. "colliculi" for token "CO"), for human-readable
%            reporting only, NOT used for pooling. '' if no clean
%            preceding phrase was found for that occurrence (harmless --
%            the caller aggregates names across many channels/sessions,
%            so a missed capture here just means one fewer example of
%            that token's full name, not a lost pooling membership).

if isempty(labelStr)
    tokens = {};
    names = {};
    return
end
matches = regexp(labelStr, '\(([A-Za-z_][A-Za-z0-9_\-]{0,11})\)', 'tokens');
tokens = cellfun(@(c) upper(c{1}), matches, 'UniformOutput', false);
tokens = unique(tokens);

names = repmat({''}, size(tokens));
% Preceding-phrase capture is bounded to at most 4 words -- callers
% concatenate several label fields together with spaces (AssignedArea +
% L6 + L5 + L4 + L3), and an unbounded quantifier here would happily span
% BACK ACROSS those field boundaries (e.g. "SC" (AssignedArea) + "SC"
% (L6) + "superior colliculus (Co)" (L5) all in the same string
% otherwise captures as "SC SC superior colliculus" for token CO).
nameMatches = regexp(labelStr, '([A-Za-z][A-Za-z0-9_]*(?:[ _][A-Za-z0-9_]+){0,2})\s*\(([A-Za-z_][A-Za-z0-9_\-]{0,11})\)', 'tokens');
for i = 1:numel(nameMatches)
    tok = upper(nameMatches{i}{2});
    nm = regexprep(nameMatches{i}{1}, '^(or|and)\s+', '', 'ignorecase');
    nm = strtrim(strrep(nm, '_', ' '));
    % Defensive cleanup: collapse immediately-repeated words
    % (case-insensitive), e.g. "SC SC superior colliculus" ->
    % "SC superior colliculus" -- can happen when concatenated label
    % fields repeat the same bare abbreviation right before a different
    % field's parenthetical version of it.
    words = strsplit(nm, ' ');
    keep = true(size(words));
    for wi = 2:numel(words)
        if strcmpi(words{wi}, words{wi-1}), keep(wi) = false; end
    end
    nm = strjoin(words(keep), ' ');
    idx = find(strcmp(tokens, tok), 1);
    if ~isempty(idx) && isempty(names{idx}) && ~isempty(nm)
        names{idx} = nm;
    end
end
end
