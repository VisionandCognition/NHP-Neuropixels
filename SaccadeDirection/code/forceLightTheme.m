function forceLightTheme(f)
% FORCELIGHTTHEME Force white figure/axes background and black text/lines
% on a figure, regardless of the MATLAB app's current light/dark theme
% setting. QC figures in this pipeline are saved as PNG via saveas, which
% follows the app theme at render time -- on a dark-themed installation
% this silently produced illegible (light-on-white) plots. Call this
% right before saveas().
%
% This must cover every graphics object type that carries its own color,
% not just axes/text: Colorbar and TiledChartLayout in particular have
% their OWN 'Color'/'Title' properties that are NOT reachable by walking
% axes children or findall(...,'Type','text') -- missing these was the
% direct cause of colorbar tick labels and tiledlayout titles staying
% pale/illegible even after axes/text were fixed.
set(f, 'Color', 'w');

axAll = [findall(f, 'Type', 'axes'); findall(f, 'Type', 'polaraxes')];
for i = 1:numel(axAll)
    if isprop(axAll(i), 'Color'), set(axAll(i), 'Color', 'w'); end
    if isprop(axAll(i), 'XColor'), set(axAll(i), 'XColor', 'k'); end
    if isprop(axAll(i), 'YColor'), set(axAll(i), 'YColor', 'k'); end
    if isprop(axAll(i), 'ZColor'), set(axAll(i), 'ZColor', 'k'); end
    if isprop(axAll(i), 'GridColor'), set(axAll(i), 'GridColor', [0.3 0.3 0.3]); end
    if isprop(axAll(i), 'MinorGridColor'), set(axAll(i), 'MinorGridColor', [0.3 0.3 0.3]); end
    if isprop(axAll(i), 'ThetaColor'), set(axAll(i), 'ThetaColor', 'k'); end
    if isprop(axAll(i), 'RColor'), set(axAll(i), 'RColor', 'k'); end
    if isprop(axAll(i), 'GridAlpha'), set(axAll(i), 'GridAlpha', 0.4); end
    ttl = get(axAll(i), 'Title');
    if ~isempty(ttl), set(ttl, 'Color', 'k'); end
    if isprop(axAll(i), 'LabelFontSizeMultiplier')
        xl = get(axAll(i), 'XLabel'); if ~isempty(xl), set(xl, 'Color', 'k'); end
        yl = get(axAll(i), 'YLabel'); if ~isempty(yl), set(yl, 'Color', 'k'); end
    end
end

% Colorbars carry their OWN Color property (controls tick labels + axis
% line) and their own Label (a Text child) -- neither is touched by the
% axes/text sweeps above.
cbAll = findall(f, 'Type', 'colorbar');
for i = 1:numel(cbAll)
    set(cbAll(i), 'Color', 'k');
    if isprop(cbAll(i), 'Label') && ~isempty(cbAll(i).Label)
        set(cbAll(i).Label, 'Color', 'k');
    end
end

% Every plain text object anywhere in the figure (free-floating text(),
% axes titles/labels already handled above but re-asserted here, and any
% annotation textboxes).
txtAll = [findall(f, 'Type', 'text'); findall(f, 'Type', 'textboxshape')];
for i = 1:numel(txtAll)
    set(txtAll(i), 'Color', 'k');
end

% Legends: text color + background.
legAll = findall(f, 'Type', 'legend');
for i = 1:numel(legAll)
    set(legAll(i), 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');
end

% TiledChartLayout title/subtitle (sgtitle-equivalent via title(tl,...))
% is stored directly on the layout object, not as an axes Title -- must
% be set explicitly, the generic 'text' findall above does not reliably
% reach inside a TiledChartLayout.
tlAll = findall(f, 'Type', 'tiledlayout');
for i = 1:numel(tlAll)
    if isprop(tlAll(i), 'Title') && ~isempty(tlAll(i).Title)
        set(tlAll(i).Title, 'Color', 'k');
    end
    if isprop(tlAll(i), 'Subtitle') && ~isempty(tlAll(i).Subtitle)
        set(tlAll(i).Subtitle, 'Color', 'k');
    end
end

% Any remaining line objects that were explicitly colored white (a common
% "visible on dark background" choice upstream) are invisible once the
% background is forced to white -- repaint pure-white lines/markers to
% black so nothing silently disappears.
lnAll = findall(f, 'Type', 'line');
for i = 1:numel(lnAll)
    if isequal(get(lnAll(i), 'Color'), [1 1 1])
        set(lnAll(i), 'Color', 'k');
    end
    mfc = get(lnAll(i), 'MarkerFaceColor');
    if ischar(mfc) == false && isequal(mfc, [1 1 1])
        set(lnAll(i), 'MarkerFaceColor', 'k');
    end
end
cstAll = findall(f, 'Type', 'constantline'); % xline/yline objects
for i = 1:numel(cstAll)
    if isequal(get(cstAll(i), 'Color'), [1 1 1])
        set(cstAll(i), 'Color', 'k');
    end
end
end
