classdef BatchProcess < Process
% while constructing, specify processor-specific initialization
% manually for your use case. Then, run the batch process.

    properties
        Processor Process = BasicMovieProcessor('subtract_baseline')  % Instance of a subclass of Process to apply
        includeFilter char = '.mat' % should have private SetAccess, and call setDataPath when updated              ## TODO ##
        excludeFilter char = '' % should have private SetAccess, and call setDataPath when updated                  ## TODO ##
        results cell
        OutFolder char = 'batch'
    end

    properties (SetAccess = private)
        DataPath = pwd;
        DataList  % List of data items (e.g., file paths, variables)
    end

    methods (Access = protected, Static)
        function [dataFiles,dataPath] = reformatDataList(item, includeFilter, excludeFilter)
            % Load data based on input type
            if isstruct(item) && isfield(item, 'name')  % Check for dir() output structure
                dataPath = item{1}.folder;
                dataFiles = {item.name};
            elseif ischar(item) || isstring(item)
                if exist(item, 'file') == 2  % Single file
                    dataFiles = {item};
                    dataPath = getFileNameSpecs(item).orig_fpath;
                elseif exist(item, 'dir') == 7  % Directory path
                    dataPath = item;
                    % Apply filters if specified
                    if exist('includeFilter', 'var') && ~isempty(includeFilter)
                        files = dir(fullfiletol(item, ['*',includeFilter]));  % List all files
                    else
                        files = dir(fullfiletol(item, '*'));
                    end
                    dataFiles = {files.name};

                    if exist('excludeFilter', 'var') && ~isempty(excludeFilter)
                        dataFiles = dataFiles(~endsWith(dataFiles, excludeFilter));
                    end
                else
                    error('Specified path does not exist: %s', item);
                end
            else
                error('Unsupported item type for loading data.');
            end
            
            % Filter for files
            dataFiles = dataFiles(isfile(fullfiletol(dataPath,dataFiles)));

            % check dataPath validity
            if isempty(dataPath) || ~isfolder(dataPath)
                dataPath = pwd;
            end
        end
        
        function processor = loadProcessor(processor,data)
            % assumes compatibility
            if isSubclassOf(processor, 'MovieProcessing')
                try
                    processor = processor.setDataRaw(data);
                catch
                    warning('this was not a Movie object! replaced with an empty Movie')
                    processor = processor.setDataRaw(Movie);            
                end
            else
                error('Behavior unspecified for this class!')
            end
        end
    end

    methods (Access = protected)
        function disp_runheader(obj)
            disp('')
            disp('###==============================###')
            disp('#------- Batch Run started --------#')
            disp('###==============================###')
            disp('')
            disp(['Hash : ', obj.hash])
            disp('')
            disp('')
        end    

        function str = genOutFolderName(obj)
            str = ['batch_', obj.hash];
        end
    end
    
    methods (Access = protected, Static)
        function data = loadData(item)            
            % Load data from one item of the file list
            data = robust_io('load',item);  % Assuming MATLAB .mat files or compatible types
        end
    end
    
    methods
        function obj = BatchProcess(processor, dataPath) % Constructor
            arguments
                processor Process
                dataPath = pwd;
            end
            obj = resetHash(obj);
            obj = setDataList(obj, dataPath);
            obj.results = {};
            obj.init.keep_heavy_res_in_memory = false;

            obj.Processor = processor;
            obj.Processor.init.superHash = obj.hash; % make sure hash 
            % attribution in the saved files is the batch hash;
        end

        function obj = resetHash(obj)
            obj = resetHash@Process(obj); % extend inherited method
            obj.OutFolder = obj.genOutFolderName();
        end

        function obj = run(obj)
            doapplyresults = obj.init.operation_precomputed;
            
            % Signal registration start
            obj.disp_runheader()
            if doapplyresults; disp('REAPPLYING PRECOMPUTED OPERATIONS TO A NEW DATASET'); end

            obj.Processor.init.batch_run = true;
            % initialize list of files, process,  ...
            nData = numel(obj.DataList);
            % Preallocate cell array to store results
            if obj.init.keep_heavy_res_in_memory && ~doapplyresults
                obj.results = cell(1, nData);
            end
            % create results directory
            outpath = fullfiletol(obj.DataPath,obj.OutFolder);
            mkdir(outpath);
            % create inits directory
            outpath_inits = fullfiletol(outpath,'inits');
            mkdir(outpath_inits);

            % call specific process FOR each file;
            for idx = 1:nData 
                data = obj.loadData(fullfiletol(obj.DataPath, obj.DataList{idx}));
                fields = fieldnames(data);
                data = data.(fields{1});
                try
                    disp([data.path.fname,'.mat'])
                catch
                    disp(['item #',num2str(idx)])
                end
                
                % load data into processor
                obj.Processor = obj.Processor.resetHash;
                obj.Processor = obj.loadProcessor(obj.Processor,data);
                if doapplyresults
                    obj.Processor.init.operation_precomputed = true;
                    obj.Processor.init.orig_computeHash = obj.results{idx}.init.hash;
                    obj.Processor.operation = obj.results{idx}.operation;
                end
                
                result = obj.Processor.run();  % Execute processing
                if obj.init.keep_heavy_res_in_memory
                    obj.results{idx} = result;
                else
                    obj.results{idx} = clearProperties(result,obj.Processor.heavyProperties);
                end
                
                % save result to disk
                result.data_processed.save(outpath,'mat') % always save results independently, 
                % including inits, in a new subfolder of the data folder;
                result.init.save(outpath_inits,'json')
                save(getLight(result),outpath_inits); % internal save fun
                
                obj.init.subHash{end+1} = result.hash;
                obj.init.subHash = obj.init.subHash(:); % make vertical
                
                % do this for each iteration to save memory!
                obj = getLight(obj,true);
            end
            
            obj.init.save(outpath,'json')
            try obj.Processor = getLight(obj.Processor); catch; end 
            save(getLight(obj,true),outpath); % internal save fun
        end        

        function obj = applyresults(obj,newdatapath)
            arguments
                obj
                newdatapath char
            end

            if isempty(obj.results)
                error('Precomputed results not found!')
            end
            obj.init.operation_precomputed = true;
            obj.init.orig_computeHash = obj.init.hash;
            obj = obj.resetHash;
            obj = obj.setDataPath(newdatapath);

            obj = obj.run();
        end
        
        function obj = setDataPath(obj, datapath)
            [obj.DataList, obj.DataPath] = obj.reformatDataList(datapath,obj.includeFilter,obj.excludeFilter);
        end

        function obj = setDataList(obj, dataList)
            obj = setDataPath(obj, dataList);
        end
    end
end
