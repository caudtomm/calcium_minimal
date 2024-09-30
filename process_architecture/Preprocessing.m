classdef Preprocessing
    properties
        sj Subject
        autosave logical = false
        pl logical = true % whether to plot stuff
    end

    methods
        function obj = Preprocessing()
            
        end

        function obj = createSubject(obj,subjectID,group,datafolder,odordelay)
            arguments
                obj 
                subjectID char
                group char
                datafolder char = 'scratch\gfriedri\caudtomm\testground'
                odordelay double = 0 % odor delay to nostril [sec]
            end

            % subjectID = 'TC_230910_TC0028_230906beh1b3_sxpDp_odorexp005_RPB3144501500AG';
            locations = Locations;
            locations = locations.setSubjectID(subjectID);
            locations = locations.setDataFolder(datafolder);
            
            % startup
            cd(locations.subject_datapath)
            
            % Create a new Subject object (this is a new fish)
            subject = Subject('fish1',locations, group);
            
            % set known attributes
            subject.odor_delay = odordelay;
            subject.notes = ""; % load notes file

            % store
            obj.sj = subject;

        end
    
        function obj = claheRawTrials(obj)
            datapath = obj.sj.locations.rawtrials;
            
            % initialize Batch Process
            b = BatchProcess(BasicMovieProcessor('clahe'));
            b.init.description = "CLAHE raw trials";
            b = b.setDataPath(datapath);            
            
            % run
            b = b.run;
            
            % move to layer 1 dir
            sourceFolder = fullfile(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.histeqtrials;
            movedirTC(sourceFolder,destinationFolder)
            
            % save a separate clahe reference
            cd(destinationFolder)
            ref = obj.sj.retrieveSnippetfromCoords(fullfile(obj.sj.locations.subject_datapath, ...
                obj.sj.locations.references.raw,'metadata.json'));
            ref.metadata.Name = "Raw Reference Image (CLAHE)";
            foutname = getFileNameSpecs(obj.sj.filelist(ref.metadata.Trial_relativenum).name).fname; % w/o extension
            fout = fullfile(obj.sj.locations.subject_datapath, ...
                obj.sj.locations.references.histeq,[foutname,'.mat']);
            movie = ref; clear ref
            save(fout,'movie','-mat'); clear movie
            cd(obj.sj.locations.subject_datapath)
        end
    
        function obj = rigidregHisteq2Raw(obj)
            % prep CLAHEd trials
            datapath = obj.sj.locations.histeqtrials;
            obj.sj.Tiff2Movie(datapath)
            
            % initialize Batch Process
            b = BatchProcess(RigidRegistration);
            b.init.description = "Rigid alignment of CLAHEd raw trials to the CLAHEd reference image";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
            b = b.setDataPath(datapath);
            imgref = obj.sj.loadReferenceImg(obj.sj.locations.references.histeq);
            b.Processor.reference_img=imgref;
            
            % run
            b = b.run;
            charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
            obj.sj.log{end+1} = charv;
            
            % move to layer 1 dir
            sourceFolder = fullfile(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.histeqtrials_rigidreg;
            movedirTC(sourceFolder,destinationFolder)
            
            
            % apply computed shifts to the raw data
            datapath = obj.sj.locations.rawtrials;
            b.init.detailed_description = "Result shifts applied to raw trials";
            b = b.applyresults(datapath);
            
            % move to layer 1 dir
            sourceFolder = fullfile(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.rawtrials_rigidreg;
            movedirTC(sourceFolder,destinationFolder)
            
            
            obj.sj = obj.sj.update_currentstate('Rigid alignment complete');
            cd(obj.sj.locations.subject_datapath)

        end
    
        function obj = rigidregRaw(obj)
            % prep CLAHEd trials
            datapath = obj.sj.locations.rawtrials;
            obj.sj.Tiff2Movie(datapath)
            
            % initialize Batch Process
            b = BatchProcess(RigidRegistration);
            b.init.description = "Rigid alignment of raw trials to the raw reference image";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
            b = b.setDataPath(datapath);
            imgref = obj.sj.reference_img;
            b.Processor.reference_img=imgref;
            
            % run
            b = b.run;
            charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
            obj.sj.log{end+1} = charv;
            
            % % move to layer 1 dir
            % sourceFolder = fullfile(datapath, b.OutFolder);
            % destinationFolder = obj.sj.locations.rawtrials_rigidreg;
            % movedirTC(sourceFolder,destinationFolder)
            
            obj.sj = obj.sj.update_currentstate('Rigid alignment complete');
            obj.sj.save2mat(obj.autosave)
            cd(obj.sj.locations.subject_datapath)
        end
    
        function obj = removeMovieBaseline(obj)
            datapath = obj.sj.locations.rawtrials_rigidreg;

            % initialize Batch Process
            b = BatchProcess(BasicMovieProcessor('subtract_baseline'));
            b.init.description = "Removing baseline from registered raw trials";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
            b = b.setDataPath(datapath);

            % run
            b = b.run;
            charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
            obj.sj.log{end+1} = charv;

            % move to layer 1 dir
            sourceFolder = fullfile(datapath, b.OutFolder);
            destinationFolder = [obj.sj.locations.rawtrials_rigidreg,'_rmbase'];
            movedirTC(sourceFolder,destinationFolder)
            obj.sj.locations.rawtrials_rigidreg = destinationFolder;


            obj.sj = obj.sj.update_currentstate('Movie baseline removed');
            cd(obj.sj.locations.subject_datapath)
        end
    
        function obj = updateSubject(obj, loc)
            arguments
                obj 
                loc char = obj.sj.locations.rawtrials_rigidreg
            end
            
            obj.sj.anatomy_imgs = obj.sj.retrieve_trial_anatomies(loc);
            obj.sj.reference_img = obj.sj.retrieve_ref_img(loc);
            obj.sj.localcorr_imgs = obj.sj.retrieve_localcorr_maps(loc);
            
            obj.sj = obj.sj.update_currentstate( ...
                'Trial anatomies and reference frame updated in Subject file');
            obj.sj.save2mat(obj.autosave)
            
            if obj.pl % plot anatomy
                SubjectViewer(obj.sj).visualize_anatomy_physFOV;
            end

            cd(obj.sj.locations.subject_datapath)
        end
    
        function obj = opticflowregRaw(obj)
            % prep raw trials
            datapath = obj.sj.locations.rawtrials;
            obj.sj.Tiff2Movie(datapath)
            
            % initialize Batch Process
            b = BatchProcess(OpticFlowRegistration);
            b.init.description = "OpticFlow alignment of raw trials to the rigidreg reference image";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
            b = b.setDataPath(datapath);
            imgref = obj.sj.reference_img;
            b.Processor.reference_img=imgref;
            b.Processor.init.original_path = fullfile(b.Processor.init.original_path, datapath);
            
            % run
            b = b.run;
            charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
            obj.sj.log{end+1} = charv;
            
            % move to layer 1 dir
            sourceFolder = fullfile(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.rawtrials_opticflowwarp;
            movedirTC(sourceFolder,destinationFolder)
            
            
            obj.sj = obj.sj.update_currentstate('Opticflow alignment complete');
            obj.sj.save2mat(obj.autosave)
            cd(obj.sj.locations.subject_datapath)
        end
    
        function obj = opticflowregHisteq(obj)
            % prep raw trials
            datapath = obj.sj.locations.histeqtrials;
            obj.sj.Tiff2Movie(datapath)
            
            % initialize Batch Process
            b = BatchProcess(OpticFlowRegistration);
            b.init.description = "OpticFlow alignment of raw trials to the rigidreg reference image";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
            b = b.setDataPath(datapath);
            imgref = obj.sj.loadReferenceImg(obj.sj.locations.references.histeq);
            b.Processor.reference_img=imgref;
            b.Processor.init.original_path = fullfile(b.Processor.init.original_path, datapath);
            
            % run
            b = b.run;
            charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
            obj.sj.log{end+1} = charv;
            
            % move to layer 1 dir
            sourceFolder = fullfile(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.histeqtrials_opticflowwarp;
            movedirTC(sourceFolder,destinationFolder)
            
            
            obj.sj = obj.sj.update_currentstate('Opticflow alignment complete');
            obj.sj.save2mat(obj.autosave)
            cd(obj.sj.locations.subject_datapath)
        end

        function obj = opticflowHisteq2Raw(obj)
            % prep CLAHEd trials
            datapath = obj.sj.locations.histeqtrials;
            obj.sj.Tiff2Movie(datapath)
            
            % initialize Batch Process
            b = BatchProcess(OpticFlowRegistration);
            b.init.description = "Opticflow alignment of CLAHEd raw trials to the CLAHEd reference image";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
            b = b.setDataPath(datapath);
            imgref = obj.sj.loadReferenceImg(obj.sj.locations.references.histeq);
            b.Processor.reference_img=imgref;
            
            % run
            b = b.run;
            charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
            obj.sj.log{end+1} = charv;
            
            % move to layer 1 dir
            sourceFolder = fullfile(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.histeqtrials_opticflowwarp;
            movedirTC(sourceFolder,destinationFolder)
            
            
            % apply computed shifts to the raw data
            datapath = obj.sj.locations.rawtrials;
            b.init.detailed_description = "Result warp applied to raw trials";
            b = b.applyresults(datapath);
            
            % move to layer 1 dir
            sourceFolder = fullfile(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.rawtrials_opticflowwarp_fromhisteq;
            movedirTC(sourceFolder,destinationFolder)
            
            
            obj.sj = obj.sj.update_currentstate('Opticflow alignment complete');
            cd(obj.sj.locations.subject_datapath)

        end

        function obj = selectROIs(obj)
            % select background ROIs
            obj.sj.backgroundROImap = imageSequenceGUI(...
                obj.sj.anatomy_imgs,...
                obj.sj.localcorr_imgs,...
                obj.sj.backgroundROImap,...
                'Select background ROIs (N->E->S->W)');
            
            % select cell ROIs
            if exist("plane",'var'); obj.sj.ROImap = plane{1}.ROI_map; end
            obj.sj.ROImap = imageSequenceGUI(...
                obj.sj.anatomy_imgs,...
                obj.sj.localcorr_imgs,...
                obj.sj.ROImap,...
                'Select neuronal somata');
            
            % save
            obj.sj = obj.sj.update_currentstate( ...
                'ROI selection complete');
            obj.sj.save2mat(obj.autosave)
        end

        function obj = extractCalciumTraces(obj)
            % address to find the source movies
            obj.sj.locations.traces_src = obj.sj.locations.rawtrials_rigidreg;
            
            loc = obj.sj.locations;
            
            PMToffmeta = fullfile(loc.subject_datapath, loc.rawtrials, 'PMToff_metadata.json');
            noLightmeta = fullfile(loc.subject_datapath, loc.rawtrials, 'noLight_metadata.json');
            obj.sj = obj.sj.defineTraces(PMToffmeta,noLightmeta);
            cd(obj.sj.locations.subject_datapath)
            
            refmovieforbackground = 'TC_240218_TC0028_240213beh1b2_sxpDp_odorexp004_RPB3144501500AG_00035_00001.tif';
            refmovieforbackground = fullfile(obj.sj.locations.general_datapath,refmovieforbackground);
            
            
            cd(obj.sj.locations.traces_src)
            
            interval = 960:1100;
            obj.sj.traces = obj.sj.traces.setPMToff(Snippet(refmovieforbackground,interval));
            interval = 1160:1300;
            obj.sj.traces = obj.sj.traces.setNoLight(Snippet(refmovieforbackground,interval));
            obj.sj.traces = obj.sj.traces.defineFundamentalProperties(obj.sj);
            obj.sj.traces = obj.sj.traces.setDerivativeProperties;
            
            cd(obj.sj.locations.subject_datapath)
            
            
            obj.sj = obj.sj.update_currentstate('Calcium traces extracted');
            obj.sj.save2mat(obj.autosave)
        end

    end
end