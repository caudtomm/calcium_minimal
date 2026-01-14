classdef PlotConfig
% Sets plotting parameters
    properties
        % Data handling

        % Visual appearance
        theme char = 'dark'          % plotting theme ('light', 'dark', etc.)
        colormapName char = 'jet'    % MATLAB colormap name
        showGrid logical = false     % display grid in internal plots
        axWidth double = 0.5         % default axis line width
        lineWidth double = 1         % default line width
        lineStyle char = '-'         % default line type
        c double = lines(100)        % default plot colors
        axcol double = [0 0 0]       % default axis color
        bgcol = [1 1 1]              % default figure background color
        textcol double = [0 0 0]     % default text color

        % Panel formatting for publication
        units char = 'centimeters'   % figure units
        panelWidth double = 8        % panel width in centimeters
        panelHeight double = 6       % panel height in centimeters
        tickLength double = 0.032    % tick length
        fontType char = 'Arial'      % font type
        fontSize double = 6          % base font size

        % Figure formatting for publication
        figSize char {mustBeMember(figSize, {'small', 'medium', 'large'})} = 'small' % panel size preset
        aspRatioType char {mustBeMember(aspRatioType, {'square', ...
                                                'tall', ...
                                                'wide'})} = 'square' % panel type determines aspect ratio
        figPos double = [10, 10, 8, 6] % figure position [left bottom width height]
        figDPI double = 600          % figure resolution in DPI
        figVectFormat char = 'svg'   % vector figure format for saving ('svg', 'pdf', etc.)
        figRasterFormat char = 'png' % raster figure format for saving ('png', 'jpg', etc.)

        % Technical settings
        renderingFactor = 72/96    % factor to convert line widths from screen (96 dpi) to print (72 dpi)

        % Optional: arbitrary metadata for caller use
        custom struct = struct()
    end

    methods
        function obj = PlotConfig(varargin)        
            if nargin == 1 && isstruct(varargin{1})
                % Initialize from struct
                s = varargin{1};
                fn = fieldnames(s);
                for i = 1:numel(fn)
                    if isprop(obj, fn{i})
                        obj.(fn{i}) = s.(fn{i});
                    end
                end
            elseif mod(nargin, 2) == 0
                % Name-value pair input
                for i = 1:2:nargin
                    name = varargin{i};
                    value = varargin{i+1};
                    if isprop(obj, name)
                        obj.(name) = value;
                    else
                        error('Invalid property name: %s', name);
                    end
                end
            elseif nargin > 0
                error('Unsupported PlotConfig constructor usage.');
            end
        end

        %% getters
        function c = get.axcol(obj)
            switch obj.theme
                case 'dark'
                    c = [1 1 1];  % white
                case 'light'
                    c = [0 0 0];  % black
                otherwise
                    c = [0.5 0.5 0.5];  % fallback gray
            end
        end
        
        function c = get.bgcol(obj)
            switch obj.theme
                case 'dark'
                    c = [0 0 0];  % dark gray
                case 'light'
                    c = [1 1 1];        % white
                otherwise
                    c = [0.95 0.95 0.95];  % light fallback
            end
        end
        
        function c = get.textcol(obj)
            switch obj.theme
                case 'dark'
                    c = [1 1 1];  % white
                case 'light'
                    c = [0 0 0];  % black
                otherwise
                    c = [0.2 0.2 0.2];  % fallback dark gray
            end
        end
        
        function colors = get.c(obj)
            % Returns 10 well-visible line/scatter colors adapted to theme
        
            switch obj.theme
                case 'dark'
                    base_colors = [
                        0.90 0.60 0.00; % amber
                        0.00 0.60 0.90; % cyan
                        0.90 0.00 0.60; % magenta
                        0.00 0.80 0.30; % green
                        0.80 0.40 0.00; % orange
                        0.50 0.50 1.00; % lavender blue
                        1.00 0.20 0.20; % red
                        0.60 0.20 1.00; % violet
                        0.20 1.00 1.00; % aqua
                        1.00 1.00 0.20  % yellow
                    ];
                    colors = repmat(base_colors, 10, 1);
                    jitter = 0.05 * randn(size(colors)); % small variation
                    colors = min(max(colors + jitter, 0), 1); % clamp between 0 and 1
                case 'light'
                    colors = lines(100); % MATLAB default
                otherwise
                    colors = colorcube(100); % fallback
            end
        end

        function cmap = getColormap(obj)
            % Return the specified colormap
            try
                cmap = feval(obj.colormapName);
            catch
                warning('Unknown colormap "%s". Falling back to "parula".', obj.colormapName);
                cmap = parula;
            end
        end

        function val = get.axWidth(obj)
            val = obj.axWidth * obj.renderingFactor;
        end

        function val = get.lineWidth(obj)
            val = obj.lineWidth * obj.renderingFactor;
        end

        function val = get.fontSize(obj)
            val = obj.fontSize * obj.renderingFactor;
        end
        
        function pos = get.figPos(obj)
            % Return the axis position within the figure
            pos = obj.figPos;

            switch obj.aspRatioType
                case 'square'
                    aspRatio = 1;
                case 'tall'
                    aspRatio = 2;
                case 'wide'
                    aspRatio = 1/2;
                otherwise
                    aspRatio = 1; % default to square
            end

            switch obj.figSize
                case 'small'
                    obj.panelWidth = 3.5;
                case 'medium'
                    obj.panelWidth = 3.5 * 2;
                case 'large'
                    obj.panelWidth = 3.5 * 3;
                otherwise
                    % default to preset
            end

            pos = [pos(1), pos(2), ...
                   obj.panelWidth * obj.renderingFactor * 0.8, ... % width
                   (obj.panelWidth * aspRatio) * obj.renderingFactor * 0.8]; % height
        end

        function saveFigure(obj, figHandle, filenameBase, format)
            % Save figure in specified formats
            arguments
                obj PlotConfig
                figHandle handle
                filenameBase char
                format char {mustBeMember(format, {'both', 'vector', 'raster'})} = 'both'
            end

            if strcmp(format, 'vector') || strcmp(format, 'both')
                % Save vector format
                print(figHandle, [filenameBase '.' obj.figVectFormat], ...
                ['-d' obj.figVectFormat], ['-r' num2str(obj.figDPI)], '-vector');
            end

            if strcmp(format, 'raster') || strcmp(format, 'both')
                % Save raster format
                print(figHandle, [filenameBase '.' obj.figRasterFormat], ...
                ['-d' obj.figRasterFormat], ['-r' num2str(obj.figDPI)], '-raster');
            end
        end
        
        function setFigure(obj)
            % Apply figure settings to current figure
            fig = gcf;
            set(fig, 'Color', obj.bgcol, ...
                     'Units', obj.units, ...
                     'Position', obj.figPos, ...
                     'PaperPositionMode', 'auto');

            % set axes for each axis (subplot) in the figure
            axesHandles = findall(fig, 'Type', 'axes');
            for i = 1:length(axesHandles)
                obj.setAxes(axesHandles(i));
            end
        end

        function setAxes(obj, axHandle)
            arguments
                obj PlotConfig
                axHandle = gca
            end

            %set(axHandle, 'Position', obj.axPos); % set position

            % find appropriate tick length
            AbsTickLength= obj.tickLength; % in cm
            pos=get(axHandle,'position');
            longaxis=max(pos(3:4));
            tickfactor=0.1/longaxis;
            Ticklengthvec=[AbsTickLength*tickfactor AbsTickLength*tickfactor];
            
            set(axHandle,'tickdir','out', ...
                    'fontsize',obj.fontSize, ...
                    'FontName', obj.fontType, ...   
                    'TitleFontWeight','normal', ...
                    'ticklength',Ticklengthvec, ...
                    'LineWidth',obj.lineWidth, ...
                    'box','off', ...
                    'XGrid', obj.showGrid, 'YGrid', obj.showGrid, 'ZGrid', obj.showGrid, ...
                    'XLimitMethod', 'tight', 'YLimitMethod', 'tight', 'ZLimitMethod', 'tight', ...
                    'LooseInset', get(gca, 'TightInset'), ...
                    'TitleFontSizeMultiplier',1,'TitleFontWeight','normal','LabelFontSizeMultiplier',1, ...
                    'XColor',obj.axcol,'YColor',obj.axcol,'ZColor',obj.axcol, 'Color',obj.bgcol);

            set(findobj(axHandle, 'Type', 'Line'), ...
                    'LineWidth', obj.lineWidth, ...
                    'LineStyle', obj.lineStyle); % set default line properties
        end
    end
end
