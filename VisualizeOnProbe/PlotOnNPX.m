function PlotOnNPX(data, definitionPath, varargin)
% PLOTONNPX Plots scalar values on a Neuropixels probe layout.
%
%   PLOTONNPX(data, definitionPath) plots values on the spatial layout 
%   of a Neuropixels probe. Supports both standard (10mm) and NHP (45mm) 
%   probe configurations via .imro or .oebin files.
%
%   INPUTS:
%       data           - Scalar values to plot. Can be:
%                        1. A vector of values (indices 1:N map to channels 0:N-1).
%                        2. An Nx2 matrix [channelIndex, value] (0-based).
%       definitionPath - Path to an Open Ephys .oebin file or an .imro 
%                        configuration file.
%
%   OPTIONAL PARAMETERS:
%       'Colormap'      - Custom colormap matrix or name (e.g., 'hot', 'magma', 'jet').
%       'CLim'          - Color limits [min, max].
%       'Robust'        - Use robust scaling (median-based) if true (default: true).
%       'Symmetric'     - Force symmetric color limits around 0 (default: auto).
%       'clabel'        - Label for the colorbar (default: 'Value').
%       'Title'         - Title for the plot (default: 'Neuropixels Probe Layout').
%       'InactiveColor' - Color for non-supplied channels/electrodes (e.g., [0.8 0.8 0.8 0.1]).
%                         Providing this enables background shank visualization.
%       'ProbeRange'    - Range of electrodes to visualize.
%                         - Scalar N: Shows electrodes 0 to N.
%                         - Vector [min max]: Shows a specific segment.
%       'ProbeType'     - Preset for probe geometry (default: 'NHP-Long').
%                         Options: 'NHP-Long' (4416 electrodes), 'NHP-Short' (1920), 
%                         'NP1.0' (960).
%       'FullProbeExtent' - Manual override for total range of electrodes to draw.
%       'ShowShank'     - Boolean to show/hide the physical shank outline (default: true).
%       'TipCoordinate' - [x, y] coordinate for the probe tip (default: [0, 0]).
%       'XScale'        - Manual multiplier for the adaptive X-scaling (default: 1).
%       'YLim'          - Custom Y-axis limits in micrometers.
%       'CollapseColumns'- Boolean to collapse 4 columns into 2 interleaved columns (default: false).
%       'Compact'       - Boolean to show a dense 4-column grid with no staggered offsets (default: true).
%
%   Note: The X-axis automatically scales to maintain a "probe-like" 
%   aspect ratio regardless of the total height of the visualized shank.
%
%   EXAMPLES:
%       % Basic usage with an oebin file
%       PlotOnNPX(data, 'structure.oebin')
%
%       % Visualize the full 45mm NHP shank with grey placeholders
%       PlotOnNPX(data, 'structure.oebin', 'InactiveColor', [0.5 0.5 0.5])
%
%       % Zoom in on a specific bank and collapse to a 2-column interleaved view
%       PlotOnNPX(data, 'config.imro', 'ProbeRange', [500 1000], 'CollapseColumns', true)
%
%   COMPATIBILITY: 
%       Requires MATLAB R2016b or later (due to jsondecode and string support).

    p = inputParser;
    addRequired(p, 'data');
    addRequired(p, 'definitionPath');
    addParameter(p, 'Colormap', []);
    addParameter(p, 'TipCoordinate', [0, 0]);
    addParameter(p, 'Robust', true); 
    addParameter(p, 'CLim', []);
    addParameter(p, 'Symmetric', []);
    addParameter(p, 'clabel', 'Value');
    addParameter(p, 'Title', 'Neuropixels Probe Layout');
    addParameter(p, 'XScale', 1);
    addParameter(p, 'YLim', []); 
    addParameter(p, 'InactiveColor', []); 
    addParameter(p, 'ProbeRange', []); 
    addParameter(p, 'ShowShank', true); 
    addParameter(p, 'FullProbeExtent', []); 
    addParameter(p, 'ProbeType', 'NHP-Long'); 
    addParameter(p, 'CollapseColumns', false); 
    addParameter(p, 'Compact', true); 
    parse(p, data, definitionPath, varargin{:});
    
    cmap = p.Results.Colormap;
    if isempty(cmap)
        n = 256; n2 = n/2;
        blue = [0, 0, 1]; white = [1, 1, 1]; red = [1, 0, 0];
        cmap = [linspace(blue(1), white(1), n2)', linspace(blue(2), white(2), n2)', linspace(blue(3), white(3), n2)';
                linspace(white(1), red(1), n2)',  linspace(white(2), red(2), n2)',  linspace(white(3), red(3), n2)'];
    elseif ischar(cmap) || isstring(cmap)
        try
            cmap = feval(cmap, 256);
        catch
            warning('PlotOnNPX: Could not evaluate colormap name "%s". Falling back to default.', cmap);
            n = 256; n2 = n/2;
            blue = [0, 0, 1]; white = [1, 1, 1]; red = [1, 0, 0];
            cmap = [linspace(blue(1), white(1), n2)', linspace(blue(2), white(2), n2)', linspace(blue(3), white(3), n2)';
                    linspace(white(1), red(1), n2)',  linspace(white(2), red(2), n2)',  linspace(white(3), red(3), n2)'];
        end
    end
    
    tipPos = p.Results.TipCoordinate;
    robust = p.Results.Robust;
    clim_custom = p.Results.CLim;
    symmetric = p.Results.Symmetric;
    cbLabel = p.Results.clabel;
    userInputXScale = p.Results.XScale;
    ylim_custom = p.Results.YLim;
    inactiveColor = p.Results.InactiveColor;
    probRangeInput = p.Results.ProbeRange;
    showShank = p.Results.ShowShank;
    fullProbeExtent = p.Results.FullProbeExtent;
    probeType = p.Results.ProbeType;
    collapseColumns = p.Results.CollapseColumns;
    compactMode = p.Results.Compact;
    
    % Determine Probe Geometry
    if strcmpi(probeType, 'NHP-Long'); maxE_type = 4416;
    elseif strcmpi(probeType, 'NHP-Short'); maxE_type = 1920;
    elseif strcmpi(probeType, 'NP1.0'); maxE_type = 960;
    else; maxE_type = 960; end
    
    if ~isempty(fullProbeExtent); actualExtent = fullProbeExtent;
    else; actualExtent = maxE_type; end
    
    fprintf('PlotOnNPX: Using %s geometry (Full Extent: %d electrodes)\n', probeType, actualExtent);
    
    
    [~, ~, ext] = fileparts(definitionPath);
    chanToElectrode = [];
    
    if strcmpi(ext, '.oebin')
        try
            jsonStr = fileread(definitionPath);
            info = jsondecode(jsonStr);
            streamIdx = -1;
            
            % info.continuous could be a struct array or a cell array
            for i = 1:length(info.continuous)
                if iscell(info.continuous); stream = info.continuous{i};
                else; stream = info.continuous(i); end
                
                if isfield(stream, 'stream_name') && contains(upper(string(stream.stream_name)), 'AP')
                    streamIdx = i; break;
                end
            end
            
            if streamIdx == -1; error('No AP stream found.'); end
            
            if iscell(info.continuous); stream = info.continuous{streamIdx};
            else; stream = info.continuous(streamIdx); end
            
            channels = stream.channels;
            if iscell(channels); chanCell = channels; else; chanCell = num2cell(channels); end
            
            chanToElectrode = NaN(length(chanCell), 1);
            for i = 1:length(chanCell)
                chan = chanCell{i};
                if isfield(chan, 'channel_metadata')
                    meta = chan.channel_metadata;
                    if isstruct(meta); metaCell = num2cell(meta); else; metaCell = meta; end
                    for m = 1:length(metaCell)
                        item = metaCell{m};
                        if ~isstruct(item); continue; end
                        fnames = fieldnames(item);
                        for f = 1:length(fnames)
                            fieldVal = item.(fnames{f});
                            if contains(lower(string(fieldVal)), 'electrode_index')
                                if isfield(item, 'value')
                                    val = item.value;
                                    if iscell(val); val = val{1}; end
                                    chanToElectrode(i) = val(1);
                                    break;
                                end
                            end
                        end
                        if ~isnan(chanToElectrode(i)); break; end
                    end
                end
            end
            fprintf('PlotOnNPX: Parsed %d electrodes from .oebin\n', sum(~isnan(chanToElectrode)));
        catch ME
            fprintf('PlotOnNPX: Warning - Failed to parse .oebin (%s). Falling back to default mapping.\n', ME.message);
        end
    elseif strcmpi(ext, '.imro')
        try
            imroStr = fileread(definitionPath);
            tokens = regexp(imroStr, '\(([^)]+)\)', 'tokens');
            chanToElectrode = NaN(384, 1);
            for i = 1:length(tokens)
                tkStr = tokens{i}{1};
                if ~isempty(strfind(tkStr, ',')); continue; end
                vals = str2num(tkStr); 
                if length(vals) >= 2
                    chan = vals(1); bank = vals(2);
                    if chan < 384; chanToElectrode(chan+1) = bank * 384 + chan; end
                end
            end
        catch; end
    end
    
    if isempty(chanToElectrode) || all(isnan(chanToElectrode))
        chanToElectrode = (0:383)';
    end
    numTotalChannels = length(chanToElectrode);

    % Data processing
    fullValues = zeros(numTotalChannels, 1);
    isSupplied = false(numTotalChannels, 1);
    if size(data, 2) == 2
        inputChans = data(:, 1); inputValues = data(:, 2);
    else
        inputValues = data(:); inputChans = (0:length(inputValues)-1)';
    end

    for i = 1:length(inputChans)
        c = inputChans(i);
        if c >= 0 && c < numTotalChannels
            val = inputValues(i);
            fullValues(c+1) = val;
            if ~isnan(val); isSupplied(c+1) = true; end
        end
    end
    
    % Scaling Logic
    valsToScale = fullValues(isSupplied);
    if isempty(valsToScale); valsToScale = 0; end
    if ~isempty(clim_custom); minVal = clim_custom(1); maxVal = clim_custom(2);
    elseif robust
        m = median(valsToScale); s = median(abs(valsToScale - m)) * 1.4826; 
        if s == 0; minVal = min(valsToScale); maxVal = max(valsToScale);
        else; minVal = m - 3 * s; maxVal = m + 3 * s; end
        minVal = max(minVal, min(valsToScale)); maxVal = min(maxVal, max(valsToScale));
    else; minVal = min(valsToScale); maxVal = max(valsToScale); end
    if isempty(symmetric); symmetric = (minVal < 0 && maxVal > 0); end
    if symmetric; limit = max(abs(minVal), abs(maxVal)); minVal = -limit; maxVal = limit; end
    if maxVal == minVal; normValues = ones(numTotalChannels, 1) * 0.5;
    else; normValues = (fullValues - minVal) / (maxVal - minVal); end
    normValues = max(0, min(1, normValues));

    % --- Range Logic ---
    % 1. Determine VIEW range (what the axes will show initially)
    if isempty(probRangeInput)
        viewE = chanToElectrode(~isnan(chanToElectrode))';
    else
        if length(probRangeInput) == 1; viewE = 0:probRangeInput;
        else; viewE = probRangeInput(1):probRangeInput(2); end
    end
    if isempty(viewE); viewE = 0:383; end
    
    % 2. Determine DRAW range (all elements drawn in the background)
    if ~isempty(inactiveColor)
        drawE = 0:actualExtent;
    else
        drawE = viewE; 
    end
    if isempty(drawE); drawE = viewE; end

    % Adaptive X-Scaling based on VIEW span
    if compactMode
        % In compact mode, 4 electrodes are stacked into one 40um physical row.
        % So the visual span is half of the standard 2-column staggering.
        spanY_view = (max(viewE)/4 - min(viewE)/4) * 40;
        baseXScale = max(1, (0.15 * spanY_view) / 60); % Matrix view needs wider default
    else
        spanY_view = (max(viewE) - min(viewE)) * 10;
        baseXScale = max(1, (0.08 * spanY_view) / 70);
    end
    if spanY_view <= 0; spanY_view = 100; end
    
    actualXScale = baseXScale * userInputXScale;
    
    fprintf('PlotOnNPX: Y-Span (View) = %.0f um, Adaptive XScale = %.1f, CLim = [%.2f, %.2f]\n', ...
            spanY_view, actualXScale, minVal, maxVal);
    
    % Determine X-Grid offsets and vertical spacing
    if compactMode
        x_offsets = [0, 15, 30, 45] * actualXScale; 
        yStep = 40; % 4 electrodes per 40um to match physical depth exactly
        eW = 15 * actualXScale; eH = 40;
    elseif collapseColumns
        x_offsets = [11, 43, 11, 43] * actualXScale; 
        yStep = 20;
        eW = 12 * actualXScale; eH = 12;
    else
        x_offsets = [11, 43, 27, 59] * actualXScale; 
        yStep = 20;
        eW = 12 * actualXScale; eH = 12;
    end
    
    % Reverse map for fast lookup during plotting
    maxActiveE = max(chanToElectrode(~isnan(chanToElectrode)));
    if isempty(maxActiveE); maxActiveE = 0; end
    revMap = NaN(max(max(drawE), maxActiveE) + 1, 1);
    validChans = find(~isnan(chanToElectrode));
    for i = 1:length(validChans)
        cIdx = validChans(i);
        eIdx = chanToElectrode(cIdx);
        if eIdx >= 0; revMap(eIdx+1) = cIdx; end
    end

    hold on;
    
    % Draw shank background (over the DRAW range)
    if showShank && ~compactMode
        shankXL = tipPos(1)*actualXScale;
        if collapseColumns; shankXR = shankXL + 54 * actualXScale; 
        else; shankXR = shankXL + 70 * actualXScale; end
        
        yMinDraw = tipPos(2) + min(drawE)/2 * 20 - 10;
        yMaxDraw = tipPos(2) + max(drawE)/2 * 20 + 20;
        
        line([shankXL, shankXL], [yMinDraw, yMaxDraw], 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        line([shankXR, shankXR], [yMinDraw, yMaxDraw], 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        
        if min(drawE) <= 0
            line([shankXL, shankXR], [yMinDraw, yMinDraw], 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    minX = Inf; maxX = -Inf; % Used only for horizontal axis padding
    
    for S = drawE
        chanIdx = NaN;
        if (S+1) <= length(revMap); chanIdx = revMap(S+1); end
        
        if compactMode
            y_row = floor(S / 4);
            x_idx = mod(S, 4) + 1;
            posX = tipPos(1)*actualXScale + x_offsets(x_idx);
            posY = tipPos(2) + y_row * yStep;
        else
            x_idx = mod(S, 4) + 1;
            y_row = floor(S / 2);
            posX = tipPos(1)*actualXScale + x_offsets(x_idx);
            posY = tipPos(2) + y_row * yStep;
        end
        
        minX = min(minX, posX - eW); maxX = max(maxX, posX + eW);

        if ~isnan(chanIdx) && isSupplied(chanIdx)
            cIdx = max(1, min(size(cmap, 1), round(normValues(chanIdx) * (size(cmap, 1) - 1)) + 1));
            patch('XData', [posX-eW/2, posX+eW/2, posX+eW/2, posX-eW/2], ...
                  'YData', [posY, posY, posY+eH, posY+eH], ...
                  'FaceColor', cmap(cIdx, :), 'EdgeColor', 'none', 'FaceAlpha', 1.0);
        elseif ~isempty(inactiveColor)
            if length(inactiveColor) == 4; rgb = inactiveColor(1:3); alpha = inactiveColor(4);
            else; rgb = inactiveColor; alpha = 0.05; end
            patch('XData', [posX-eW/2, posX+eW/2, posX+eW/2, posX-eW/2], ...
                  'YData', [posY, posY, posY+eH, posY+eH], ...
                  'FaceColor', rgb, 'EdgeColor', [0.8 0.8 0.8], 'EdgeAlpha', 0.1, 'FaceAlpha', alpha);
        end
    end
    
    % Final axis control based on VIEW range
    if isempty(ylim_custom)
        if compactMode
            yMinView = tipPos(2) + floor(min(viewE)/4) * 40 - 100;
            yMaxView = tipPos(2) + floor(max(viewE)/4) * 40 + 100;
        else
            yMinView = tipPos(2) + min(viewE)/2 * 20 - 100;
            yMaxView = tipPos(2) + max(viewE)/2 * 20 + 100;
        end
    else
        yMinView = ylim_custom(1); yMaxView = ylim_custom(2);
    end
    
    ylabel('Y (\mum)'); title(p.Results.Title, 'Units', 'normalized', 'Position', [0.5, 1.02, 0]);
    set(gca, 'XTick', [], 'FontSize', 12); axis normal; grid on;
    if ~isinf(minX)
        % Use periodic padding based on the INTRINSIC width to make XScale responsive.
        % This ensures that as actualXScale increases, the probe takes up more 
        % physical space on the screen relative to the padding.
        intrinsicWidth = (maxX - minX);
        if userInputXScale > 0; intrinsicWidth = intrinsicWidth / userInputXScale; end
        
        axis([minX - 0.2*intrinsicWidth, maxX + 0.5*intrinsicWidth, yMinView, yMaxView]);
    end
    colormap(cmap); c = colorbar; c.Title.String = cbLabel; c.FontSize = 11; 
    % Make colorbar 1/3 height and top-aligned with safe tick landing
    c.Location = 'eastoutside';
    if maxVal ~= minVal; caxis([minVal maxVal]); end
    
    % Force explicit ticks and rounded labels
    tickSteps = [minVal, (minVal+maxVal)/2, maxVal];
    c.Ticks = tickSteps;
    % Round to 2 decimals for the display labels
    c.TickLabels = arrayfun(@(v) sprintf('%.2f', v), tickSteps, 'UniformOutput', false);
    
    pause(0.05); % Brief pause to let MATLAB calculate the auto-position
    cPos = c.Position;
    cPos(4) = cPos(4) * 0.3; % 30% of axis height
    cPos(2) = cPos(2) + cPos(4) * 2.2; % Shift to top with extra padding for title
    cPos(1) = cPos(1) - 0.08; % Pull in to the left more significantly
    c.Position = cPos;
    
    hold off;
end
