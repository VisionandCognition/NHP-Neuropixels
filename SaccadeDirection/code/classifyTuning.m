function [category, color] = classifyTuning(pPerm, pKW, alpha)
% CLASSIFYTUNING Categorize a channel's direction-response pattern using
% two different tests that are sensitive to different SHAPES of tuning,
% not just "tuned vs not":
%
%   - pPerm : permutation test on the vector-sum resultant length (used
%     as the sole "tuned" criterion before 2026-07-08). This is only
%     sensitive to a UNIMODAL pattern -- one clear preferred direction.
%     A channel that responds to two opposite directions (bimodal), or
%     broadly to several adjacent directions but not in a way that sums
%     to a strong single vector, can have a perfectly real, structured,
%     significant direction effect and still fail this test, because
%     opposing vector contributions cancel in the vector sum.
%   - pKW : Kruskal-Wallis test across the 8 direction groups (omnibus:
%     does response differ across directions AT ALL, regardless of
%     shape). Catches bimodal/multimodal/complex patterns that pPerm
%     misses, at the cost of not saying anything about a preferred
%     direction.
%
% Combining them gives three categories instead of a binary tuned/not:
%   'directional' : pPerm significant -- has a clear single preferred
%                   direction (prefDir/resultantLen are meaningful).
%   'complex'     : pPerm NOT significant but pKW IS -- some real,
%                   direction-dependent structure exists, but not a
%                   single-peaked one (prefDir/resultantLen are NOT a
%                   meaningful summary for these; look at meanByDir
%                   directly instead).
%   'untuned'     : neither test significant.
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
category(kwSig) = {'complex'};
category(dirSig) = {'directional'}; % directional takes priority when both are significant

colDirectional = [0.85 0.10 0.10]; % red
colComplex     = [0.90 0.55 0.00]; % orange
colUntuned     = [0.70 0.70 0.70]; % gray

color = repmat(colUntuned, n, 1);
color(kwSig, :) = repmat(colComplex, sum(kwSig), 1);
color(dirSig, :) = repmat(colDirectional, sum(dirSig), 1);

if n == 1
    category = category{1};
end
end
