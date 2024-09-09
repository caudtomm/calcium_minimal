classdef Movie
    % A Movie Object
    properties
        description string = ""
        detailed_description string = ""

        stack double
        path
        timeavg double
        badperiods double = []
        log cell = {}

        % initially copied from scanimage_meta for quicker access
        h double = 0    % frame height
        w double = 0    % frame width
        nfr double = 1  % number of frames
        scanimage_meta cell
        metadata struct
    end
    properties (SetAccess = immutable)
        fs double       % framerate
        frameperiod double = 0
        imagingFovUm double = [0 0; 0 0; 0 0; 0 0]  % frame size [um]
    end

    methods
        function obj = Movie(src) % Constructor
            % arg 'src' is filename or numerical matrix
            arguments
                src
            end

            if nargin>0 && ~isempty(src)
                [obj.stack, obj.scanimage_meta, obj.path] = readSource(src);
            end

            % copy properties from 'scanimage_meta'
            try
                obj.h = obj.scanimage_meta{1}.SI_hRoiManager_linesPerFrame;
                obj.w = obj.scanimage_meta{1}.SI_hRoiManager_pixelsPerLine;
                obj.nfr = numel(obj.scanimage_meta);
                obj.fs = obj.scanimage_meta{1}.SI_hRoiManager_scanFrameRate;
                obj.frameperiod = obj.scanimage_meta{1}.SI_hRoiManager_scanFramePeriod;
                obj.imagingFovUm = obj.scanimage_meta{1}.SI_hRoiManager_imagingFovUm;
            catch
                [obj.h,obj.w,obj.nfr] = size(obj.stack);
            end
        end
        
        function value = get.timeavg(obj)
            value = mean(obj.stack,3,'omitmissing');
        end

        function obj = checkAgain(obj)
            %initialize
            movie = obj.stack;
            
            % re-compute
            [height,width,framecount] = size(movie);
            zproj = mean(movie,3,'omitmissing');

            % store
            obj.h = height;
            obj.w = width;
            obj.nfr = framecount;
            obj.timeavg = zproj;
        end
        
        % Method to log operations and their hash address
        function obj = update_log(obj,charv)
            % append a new voice to log
            obj.log{end+1} = charv;
            % make log vertical for easier visualization (only needs to happen once)
            if numel(obj.log)==2; obj.log = transpose(obj.log); end
        end

        function save(obj, newpath, type)
            arguments
                obj Movie
                newpath char = ''
                type char = 'mat'
            end
            
            outpath = pwd;
            if ~isempty(newpath)
                outpath = newpath;
            elseif ~isempty(obj.path) && ~isempty(obj.path.orig_fpath)
                outpath = obj.path.orig_fpath;
            end

            outfname = 'movie';
            if ~isempty(obj.path) && ~isempty(obj.path.fname)
                outfname = obj.path.fname;
            end

            switch type
                case 'mat'
                    FileOut = fullfile(outpath,[outfname,'.mat']);
                    b = prompt_overwrite(FileOut);
        
                    obj.path = getFileNameSpecs(FileOut);
        
                    movie = obj;
                    if b
                        save(FileOut,'movie');
                    end
                    
                case 'tif'
                    FileOut = fullfile(outpath,[outfname,'.tif']);
                    b = prompt_overwrite(FileOut);
        
                    movie = obj.stack;
                    if b
                        saveastiff(uint16(movie), FileOut)
                    end

                otherwise
                    error(['File type ''',type,''' not supported for saving.'])
            end
        end
    end
end

function [stack,scanimage_meta,path] = readSource(src)
    % initialize
    scanimage_meta = {};
    [stack, path] = deal([]);

    % read
    switch class(src)
        case {'string','char'}  % src : filename
            src = char(src);
            disp(['Loading: ',src])
            if endsWith(src,'.mat')
                load(src,'movie');
                stack = movie.stack;
                scanimage_meta = movie.scanimage_meta;
                path = movie.path;
            else
                [stack, scanimage_meta] = loadTiffStack(src);
                stack = double(stack);
                scanimage_meta = convertScanimageMeta(scanimage_meta);
                path = getFileNameSpecs(src);
            end
            if isempty(path.orig_fpath)
                path.orig_fpath = pwd;
            end
        case {'double'}         % src : numeric matrix defining the movie
            stack = src;
    end
end

function result = convertScanimageMeta(char_meta)
    result = char_meta;
    for i = 1:numel(char_meta)
        result{i} = parseToStruct(char_meta{i});
    end
end

function result = parseToStruct(inputStr)
    % Initialize an empty structure
    result = struct();
    
    % Split the input string into lines
    lines = strsplit(inputStr, '\n');
    
    % Process each line
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if isempty(line)
            continue;
        end
        
        % Find the equals sign to split parameter names and values
        eqIndex = find(line == '=', 1);
        if isempty(eqIndex)
            continue;
        end
        
        % Extract the parameter name and value
        paramName = strtrim(line(1:eqIndex-1));
        paramValue = strtrim(line(eqIndex+1:end));
        
        % Check if the value is numeric, a string, or an array
        if ~isempty(paramValue) && all(ismember(paramValue(1), '0123456789+-.[],'))
            % Try to evaluate the expression to handle numbers and vectors
            try
                paramValue = eval(paramValue);
            catch
                % If eval fails, revert to treating it as a string
                paramValue = strrep(paramValue, '''', ''); % Remove any enclosing single quotes
            end
        else
            % It's a string, strip single quotes
            paramValue = strrep(paramValue, '''', '');
        end
        
        % Assign the value to the structure
        paramName = matlab.lang.makeValidName(paramName); % Ensure the field name is a valid MATLAB identifier
        result.(paramName) = paramValue;
    end
end
