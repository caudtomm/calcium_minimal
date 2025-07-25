classdef PlotConfig
% Sets plotting parameters
    properties
        % Data handling

        % Visual appearance
        theme char = 'dark'          % plotting theme ('light', 'dark', etc.)
        colormapName char = 'jet'    % MATLAB colormap name
        showGrid logical = true      % display grid in internal plots
        lineWidth double = 1.5       % default line width
        c double = lines(100)        % default plot colors
        axcol double = [0 0 0]       % default axis color
        bgcol = [1 1 1]              % default figure background color
        textcol double = [0 0 0]     % default text color

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

    end
end
