function [category, color] = classifyTuning(pPerm, pKW, alpha)
% CLASSIFYTUNING Categorize a channel's direction-response pattern using
% two different tests that are sensitive to different SHAPES of tuning,
% not just "tuned vs not":
%
%   - pPerm : permutation test on the vector-sum resultant length (used
%     as the sole "tuned" criterion before 2026-07-08). This is only
%     GUARANTEED to fail for a perfectly symmetric bimodal pattern (two
%     opposite directions, EXACTLY equal magnitude) -- any asymmetry in
%     peak height or angular separation leaves a nonzero residual vector,
%     which can easily be large enough to pass this test even though the
%     underlying curve visibly has two lobes, not one. So "pPerm
%     significant" does NOT by itself mean "single clean peak."
%   - pKW : Kruskal-Wallis test across the 8 direction groups (omnibus:
%     does response differ across directions AT ALL, regardless of
%     shape). Catches bimodal/multimodal/complex patterns, including many
%     that ALSO pass pPerm (see above) -- checked directly against
%     2026-07-14's tuning_results.mat: 69% of pPerm-significant channels
%     also pass pKW, and 49% have a secondary local peak at least 90 deg
%     from the primary one and at least half its height, i.e. a
%     genuinely two-lobed shape is the COMMON case among pPerm-significant
%     channels, not a rare edge case.
%
% Before 2026-07-14 these were combined with pPerm taking priority
% whenever both were significant, giving only 3 categories and silently
% discarding the pKW result for every channel where a real dominant peak
% coexisted with real secondary structure -- readers were told
% prefDir/resultantLen are "meaningful" for every 'directional' channel,
% which wasn't true for the ~69% majority with a significant secondary
% lobe. Fixed: the two tests are now crossed directly (2x2), giving FOUR
% categories:
%   'directional'         : pPerm significant, pKW NOT significant --
%                            clean single-peaked response; prefDir/
%                            resultantLen are a reliable summary.
%   'directional_complex' : BOTH significant -- a real dominant peak
%                            exists (prefDir/resultantLen still computed
%                            and usable AS A DOMINANT-DIRECTION ESTIMATE)
%                            but real secondary structure also exists
%                            (a second, weaker-but-real peak elsewhere) --
%                            treat prefDir as "the stronger of (at least)
%                            two real preferences," not a clean summary.
%   'complex'              : pKW significant, pPerm NOT -- real
%                            direction-dependent structure with no
%                            reliable single dominant peak; prefDir/
%                            resultantLen are NOT meaningful for these,
%                            look at meanByDir directly instead.
%   'untuned'              : neither test significant.
%
% Membership-preserving for existing counts: 'directional' (old) ==
% 'directional' (new) UNION 'directional_complex' (new) by construction,
% so any code that needs the OLD "has some real preferred direction, of
% either purity" semantics should test
% ismember(category, {'directional','directional_complex'}) rather than
% strcmp(category,'directional') alone.
%
% [category, color] = classifyTuning(pPerm, pKW, alpha)
%   pPerm, pKW : p-values (scalar or array, elementwise)
%   alpha      : significance threshold (default 0.05)
%   color      : Nx3 RGB for direct use in plotting

if nargin < 3, alpha = 0.05; end

dirSig = pPerm < alpha;
kwSig = pKW < alpha;

n = numel(pPerm);
category = repmat({'untuned'}, n, 1);
category(kwSig & ~dirSig) = {'complex'};
category(dirSig & ~kwSig) = {'directional'};
category(dirSig & kwSig) = {'directional_complex'};

colDirectional        = [0.85 0.10 0.10]; % red
colDirectionalComplex = [0.60 0.10 0.75]; % purple -- dominant peak + real secondary structure
colComplex            = [0.90 0.55 0.00]; % orange
colUntuned            = [0.70 0.70 0.70]; % gray

color = repmat(colUntuned, n, 1);
color(kwSig & ~dirSig, :) = repmat(colComplex, sum(kwSig & ~dirSig), 1);
color(dirSig & ~kwSig, :) = repmat(colDirectional, sum(dirSig & ~kwSig), 1);
color(dirSig & kwSig, :) = repmat(colDirectionalComplex, sum(dirSig & kwSig), 1);

if n == 1
    category = category{1};
end
end
