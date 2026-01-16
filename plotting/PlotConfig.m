classdef PlotConfig
% Sets plotting parameters
    properties
        % Data handling

        % Visual appearance
        theme char = 'light'         % plotting theme ('light', 'dark', etc.)
        colormapName char = 'lapaz'  % MATLAB colormap name
        showGrid logical = false     % display grid in internal plots
        axWidth double = 0.5         % default axis line width
        lineWidth double = 1         % default line width
        lineStyle char = '-'         % default line type
        c double = lines(100)        % default plot colors
        favouriteColors double = 1:10 % indices of favorite colors from colormap
        axcol double = [0 0 0]       % default axis color
        bgcol = [1 1 1]              % default figure background color
        textcol double = [0 0 0]     % default text color

        % Panel formatting for publication
        units char = 'centimeters'   % figure units
        axPos double = [.2, .2, .7, .7] % axis position within figure [left bottom width height]
        baseSize double = 3          % base panel size in centimeters
        panelWidth double = 8        % panel width in centimeters
        panelHeight double = 6       % panel height in centimeters
        tickLength double = 0.032    % tick length
        fontType char = 'Arial'      % font type
        fontSize double = 6          % base font size

        % Figure formatting for publication
        figSize char {mustBeMember(figSize, {'tiny', 'small', 'medium', 'large'})} = 'small' % panel size preset
        aspRatioType char {mustBeMember(aspRatioType, {'square', ...
                                                'tall', ...
                                                'wide'})} = 'square' % panel type determines aspect ratio
        figPos double = [10, 10, 8, 6] % figure position [left bottom width height]
        figDPI double = 600          % figure resolution in DPI
        figVectFormat char = 'svg'   % vector figure format for saving ('svg', 'pdf', etc.)
        figRasterFormat char = 'png' % raster figure format for saving ('png', 'jpg', etc.)
        savePath char = ''            % path to save figures

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

        function cmap = getColormap(obj, mode)
            arguments
                obj PlotConfig
                mode char {mustBeMember(mode, {'categorical', ...
                                                'discrete10', ...
                                                'discrete25', ...
                                                'discrete50', ...
                                                'discrete100', ...
                                                'continuous'})} = 'continuous'
            end

            % Return the colormap based on the theme and colormapName
            switch mode
                case 'categorical'
                    suffix = 'S';
                case 'discrete10'
                    suffix = '10';
                case 'discrete25'
                    suffix = '25';
                case 'discrete50'
                    suffix = '50';
                case 'discrete100'
                    suffix = '100';
                case 'continuous'
                    suffix = '';
                otherwise
                    suffix = '';
            end

            thisname = [obj.colormapName,suffix];

            cmap = feval(thisname);
        end
        
        function colors = get.c(obj)
            % Some well-visible line/scatter colors adapted to theme
        
            try
                colors = obj.getColormap('categorical');
                colors = colors(obj.favouriteColors, :);
                return
            catch
            end

            % the following functions as a catch
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
                    % Nature standard (colorblind friendly)
                    rgb_vals = [
                        0 0 0; % black
                        182 219 255; % light blue
                        123 176 223; % mid blue
                        25 100 176; % dark blue
                        0 201 146; % light teal
                        0 138 105; % teal
                        56 99 80; % dark teal
                        233 220 109; % yellow
                        244 166 55; % orange
                        219 88 41; % vermillion
                        137 75 69; % maroon
                        210 187 215; % light purple
                        174 117 162; % purple
                        136 45 113; % dark purple
                        222 222 222; % grey
                    ];
                    colors = rgb_vals / 255; % normalize to [0, 1]
                otherwise
                    colors = colorcube(100); % fallback
                    % or
                    % colors = lines(100); % MATLAB default
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
        
        function val = get.panelWidth(obj)
            switch obj.figSize
                case 'tiny'
                    val = obj.baseSize/2;
                case 'small'
                    val = obj.baseSize;
                case 'medium'
                    val = obj.baseSize * 1.5;
                case 'large'
                    val = obj.baseSize * 3;
                otherwise
                    % default to preset
            end
        end
        
        function val = get.panelHeight(obj)
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

            val = obj.panelWidth * aspRatio;
        end

        function pos = get.figPos(obj)
            % Return the axis position within the figure
            pos = obj.figPos;

            pos = [pos(1), pos(2), ...
                   obj.panelWidth * obj.renderingFactor, ... % width
                   obj.panelHeight * obj.renderingFactor]; % height
        end

        function saveFigure(obj, figHandle, filenameBase, format)
            % Save figure in specified formats
            arguments
                obj PlotConfig
                figHandle handle
                filenameBase char
                format char {mustBeMember(format, {'both', 'vector', 'raster'})} = 'both'
            end

            if ~isempty(obj.savePath) && ~isfolder(obj.savePath)
                mkdir(obj.savePath);
            end

            filenameBase = fullfiletol(obj.savePath, filenameBase);

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
                     'PaperUnits', obj.units, ...
                     'PaperPosition', obj.figPos, ...
                     'PaperPositionMode', 'manual');

            % set axes for each axis (subplot) in the figure
            axesHandles = findall(fig, 'Type', 'axes');
            for i = 1:length(axesHandles)
                obj.setAxes(axesHandles(i));
                obj.setAxes(axesHandles(i)); % for some reason needs to be called twice to apply properly
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
                    'LineWidth',obj.axWidth, ...
                    'box','off', ...
                    'XGrid', obj.showGrid, 'YGrid', obj.showGrid, 'ZGrid', obj.showGrid, ...
                    'XLimitMethod', 'tight', 'YLimitMethod', 'tight', 'ZLimitMethod', 'tight', ...
                    'LooseInset', get(gca, 'TightInset'), ...
                    'TitleFontSizeMultiplier',1,'TitleFontWeight','normal','LabelFontSizeMultiplier',1, ...
                    'XColor',obj.axcol,'YColor',obj.axcol,'ZColor',obj.axcol, 'Color',obj.bgcol);

            set(findobj(axHandle, 'Type', 'Line'), ...
                    'LineWidth', obj.lineWidth, ...
                    'LineStyle', obj.lineStyle); % set default line properties

            if ~isempty(findobj(axHandle, 'Type', 'Image'))
                colormap(axHandle, obj.getColormap()); % set colormap if images are present
            end

            textHandles = obj.getTextinAxis(axHandle);
            set(textHandles, 'FontSize', obj.fontSize, 'FontName', obj.fontType); % set text properties
        end

        function allHandles = getTextinAxis(~, axHandle)
            textHandles = findall(axHandle, 'Type', 'text', '-or', ...
                                 'Type', 'xlabel', '-or', ...
                                 'Type', 'ylabel', '-or', ...
                                 'Type', 'zlabel', '-or', ...
                                 'Type', 'title', '-or', ...
                                 'Type', 'legend', '-or', ...
                                 'Type', 'annotation');

            cbHandles = findall(ancestor(axHandle, 'figure'), 'Type', 'colorbar');

            cbLabelHandles = [];
            if ~isempty(cbHandles) && all(isgraphics(cbHandles))
                % Use [cbHandles.Label] to get the text objects for the labels
                cbLabelHandles = [cbHandles.Label]'; 
            end

            allHandles = [textHandles; cbLabelHandles; cbHandles];
        end
    end
end
