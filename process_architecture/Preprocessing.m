classdef Preprocessing
    properties
        sj Subject
        autosave logical = false
        pl logical = true % whether to plot stuff
    end

    methods (Static)
        function ROImap_out = pruneROIs(ROImap_in)
            arguments
                ROImap_in double
            end

            % initialize output
            ROImap_out = zeros(size(ROImap_in));

            labs = unique(ROImap_in); labs(labs==0) = [];
            nrois = numel(labs);

            for i_roi = 1:nrois
                thisroi = ROImap_in == labs(i_roi);
                
                ROImap_out(thisroi) = i_roi;
            end
        end
    end

    methods (Access = protected)
        function obj = tail_sequence(obj,destinationFolder)
            
            % save anatomy images for each trial
            obj.sj.saveTrialAvgs(destinationFolder);
            % save anatomical reference, based on the results of this
            % method
            obj.sj.retrieve_ref_img(destinationFolder,true);
            % save a side-by-side subsampled avi movie of the results for
            % all the trials, for easy supervicial visual inspection
            RegistrationViewer(obj.sj,destinationFolder). ...
                subsampledMovieCollage.save(destinationFolder,'avi','sidebyside');
        end
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
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            datapath = obj.sj.locations.rawtrials;
            
            % initialize Batch Process
            b = BatchProcess(BasicMovieProcessor('clahe'));
            b.init.description = "CLAHE raw trials";
            b = b.setDataPath(datapath);            
            
            % run
            b = b.run;
            
            % move to layer 1 dir
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.histeqtrials;
            movedirTC(sourceFolder,destinationFolder)
            
            % save a separate clahe reference
            cd(destinationFolder)
            ref = obj.sj.retrieveSnippetfromCoords(fullfiletol(obj.sj.locations.subject_datapath, ...
                obj.sj.locations.references.raw,'metadata.json'));
            ref.metadata.Name = "Raw Reference Image (CLAHE)";
            foutname = getFileNameSpecs(obj.sj.filelist(ref.metadata.Trial_relativenum).name).fname; % w/o extension
            fout = fullfiletol(obj.sj.locations.subject_datapath, ...
                obj.sj.locations.references.histeq,[foutname,'.mat']);
            movie = ref; clear ref
            save(fout,'movie','-mat'); clear movie
            cd(fullfiletol(obj.sj.locations.subject_datapath))

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
        end

        function [obj, b] = rigidregHisteq(obj)
            % prep CLAHEd trials
            cd(fullfiletol(obj.sj.locations.subject_datapath))
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
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.histeqtrials_rigidreg;
            movedirTC(sourceFolder,destinationFolder)

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
        end
    
        function [obj, b] = rigidregHisteq2Raw(obj)
            [obj, b] = obj.rigidregHisteq();
            
            % apply computed shifts to the raw data
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            datapath = obj.sj.locations.rawtrials;
            b.init.detailed_description = "Result shifts applied to raw trials";
            b = b.applyresults(datapath);
            
            % move to layer 1 dir
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.rawtrials_rigidreg_fromhisteq;
            movedirTC(sourceFolder,destinationFolder)

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
            
            
            obj.sj = obj.sj.update_currentstate('Rigid alignment complete');
            cd(fullfiletol(obj.sj.locations.subject_datapath))

        end
    
        function [obj, b] = rigidregRaw(obj)
            % prep CLAHEd trials
            cd(fullfiletol(obj.sj.locations.subject_datapath))
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
            
            % move to layer 1 dir
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.rawtrials_rigidreg;
            movedirTC(sourceFolder,destinationFolder)

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
            
            obj.sj = obj.sj.update_currentstate('Rigid alignment complete');
            obj.sj.save2mat(obj.autosave)
            cd(fullfiletol(obj.sj.locations.subject_datapath))
        end
    
        function [obj, b] = removeMovieBaseline(obj)
            cd(fullfiletol(obj.sj.locations.subject_datapath))
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
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = [obj.sj.locations.rawtrials_rigidreg,'_rmbase'];
            movedirTC(sourceFolder,destinationFolder)
            obj.sj.locations.rawtrials_rigidreg = destinationFolder;

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);


            obj.sj = obj.sj.update_currentstate('Movie baseline removed');
            cd(fullfiletol(obj.sj.locations.subject_datapath))
        end
    
        function obj = updateSubject(obj, loc)
            arguments
                obj 
                loc char = obj.sj.locations.rawtrials_rigidreg
            end

            cd(fullfiletol(obj.sj.locations.subject_datapath))
            
            obj.sj.anatomy_imgs = obj.sj.retrieve_trial_anatomies(loc);
            obj.sj.reference_img = obj.sj.retrieve_ref_img(loc);
            obj.sj.localcorr_imgs = obj.sj.retrieve_localcorr_maps(loc);

            % store as Subject traces source location
            obj.sj.locations.traces_src = loc;
            obj.sj = obj.sj.update_currentstate( ...
                'Trial anatomies and reference frame updated in Subject file');
            obj.sj.save2mat(obj.autosave)
            
            if obj.pl % plot anatomy
                SubjectViewer(obj.sj).visualize_anatomy_physFOV;
            end

            cd(fullfiletol(obj.sj.locations.subject_datapath))
        end
    
        function [obj, b] = opticflowregRaw(obj)
            % prep raw trials
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            datapath = obj.sj.locations.rawtrials;
            obj.sj.Tiff2Movie(datapath)
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            
            % initialize Batch Process
            b = BatchProcess(OpticFlowRegistration);
            b.init.description = "OpticFlow alignment of raw trials to the rigidreg reference image";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
            b = b.setDataPath(datapath);
            imgref = obj.sj.reference_img;
            b.Processor.reference_img=imgref;
            b.Processor.init.original_path = fullfiletol(b.Processor.init.original_path, datapath);
            
            % run
            b = b.run;
            charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
            obj.sj.log{end+1} = charv;
            
            % move to layer 1 dir
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.rawtrials_opticflowwarp;
            movedirTC(sourceFolder,destinationFolder)

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
            
            
            obj.sj = obj.sj.update_currentstate('Opticflow alignment complete');
            obj.sj.save2mat(obj.autosave)
            cd(fullfiletol(obj.sj.locations.subject_datapath))
        end
    
        function [obj, b] = opticflowregHisteq(obj)
            % prep raw trials
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            datapath = obj.sj.locations.histeqtrials;
            obj.sj.Tiff2Movie(datapath)
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            
            % initialize Batch Process
            b = BatchProcess(OpticFlowRegistration);
            b.init.description = "OpticFlow alignment of raw trials to the rigidreg reference image";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
            b = b.setDataPath(datapath);
            imgref = obj.sj.loadReferenceImg(obj.sj.locations.references.histeq);
            b.Processor.reference_img=imgref;
            b.Processor.init.original_path = fullfiletol(b.Processor.init.original_path, datapath);
            
            % run
            b = b.run;
            charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
            obj.sj.log{end+1} = charv;
            
            % move to layer 1 dir
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.histeqtrials_opticflowwarp;
            movedirTC(sourceFolder,destinationFolder)

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
            
            
            obj.sj = obj.sj.update_currentstate('Opticflow alignment complete');
            obj.sj.save2mat(obj.autosave)
            cd(fullfiletol(obj.sj.locations.subject_datapath))
        end

        function [obj, b] = opticflowHisteq2Raw(obj)
            % % prep CLAHEd trials
            % datapath = obj.sj.locations.histeqtrials;
            % obj.sj.Tiff2Movie(datapath)
            % 
            % % initialize Batch Process
            % b = BatchProcess(OpticFlowRegistration);
            % b.init.description = "Opticflow alignment of CLAHEd raw trials to the CLAHEd reference image";
            % b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
            % b = b.setDataPath(datapath);
            % imgref = obj.sj.loadReferenceImg(obj.sj.locations.references.histeq);
            % b.Processor.reference_img=imgref;
            % 
            % % run
            % b = b.run;
            % charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
            % obj.sj.log{end+1} = charv;
            % 
            % % move to layer 1 dir
            % sourceFolder = fullfiletol(datapath, b.OutFolder);
            % destinationFolder = obj.sj.locations.histeqtrials_opticflowwarp;
            % movedirTC(sourceFolder,destinationFolder)
            [obj,b] = obj.opticflowregHisteq();
            
            % apply computed shifts to the raw data
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            datapath = obj.sj.locations.rawtrials;
            b.init.detailed_description = "Result warp applied to raw trials";
            b = b.applyresults(datapath);
            
            % move to layer 1 dir
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.rawtrials_opticflowwarp_fromhisteq;
            movedirTC(sourceFolder,destinationFolder)

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
            
            
            obj.sj = obj.sj.update_currentstate('Opticflow alignment complete');
            cd(fullfiletol(obj.sj.locations.subject_datapath))

        end

        function obj = allRegistration(obj)
            % every registration protocol, done sequentially
            obj = obj.rigidregRaw;
            % then
            obj = obj.rigidregHisteq2Raw;
            % then
            obj = obj.opticflowregRaw;
            % then
            obj = obj.opticflowHisteq2Raw;
        end

        function obj = allRigidReg(obj)
            % every rigid registration protocol, done sequentially
            obj = obj.rigidregRaw;
            % then
            obj = obj.rigidregHisteq2Raw;
        end

        function obj = selectROIs(obj)
            % select background ROIs
            obj.sj.backgroundROImap = imageSequenceGUI(...
                obj.sj.anatomy_imgs,...
                obj.sj.localcorr_imgs,...
                obj.sj.backgroundROImap,...
                'Select background ROIs (N->E->S->W)');
            obj.sj.backgroundROImap = obj.pruneROIs(obj.sj.backgroundROImap);
            
            % select cell ROIs
            if exist("plane",'var'); obj.sj.ROImap = plane{1}.ROI_map; end
            obj.sj.ROImap = imageSequenceGUI(...
                obj.sj.anatomy_imgs,...
                obj.sj.localcorr_imgs,...
                obj.sj.ROImap,...
                'Select neuronal somata');
            obj.sj.ROImap = obj.pruneROIs(obj.sj.ROImap);
            
            % save
            obj.sj = obj.sj.update_currentstate( ...
                'ROI selection complete');
            obj.sj.save2mat(obj.autosave)
        end

        function obj = extractCalciumTraces(obj, do_save)    
            arguments
                obj 
                do_save logical = true
            end
            loc = obj.sj.locations;

            % temporary
            obj.sj.backgroundROImap = obj.pruneROIs(obj.sj.backgroundROImap);
            obj.sj.ROImap = obj.pruneROIs(obj.sj.ROImap);
            
            PMToffmeta = fullfiletol(loc.subject_datapath, loc.rawtrials, 'PMToff_metadata.json');
            noLightmeta = fullfiletol(loc.subject_datapath, loc.rawtrials, 'noLight_metadata.json');
            obj.sj = obj.sj.defineTraces(PMToffmeta,noLightmeta);
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            
            refmovieforbackground = 'TC_240218_TC0028_240213beh1b2_sxpDp_odorexp004_RPB3144501500AG_00035_00001.tif';
            refmovieforbackground = fullfiletol(obj.sj.locations.general_datapath,refmovieforbackground);
            
            
            cd(fullfiletol(obj.sj.locations.traces_src))
            
            interval = 960:1100;
            obj.sj.traces = obj.sj.traces.setPMToff(Snippet(refmovieforbackground,interval));
            interval = 1160:1300;
            obj.sj.traces = obj.sj.traces.setNoLight(Snippet(refmovieforbackground,interval));
            obj.sj.traces = obj.sj.traces.defineFundamentalProperties(obj.sj);
            obj.sj.traces = obj.sj.traces.setDerivativeProperties;
            
            cd(fullfiletol(obj.sj.locations.subject_datapath))

            if do_save
                obj.sj.traces.save('','full',obj.autosave);
                obj.sj.traces.save('','light',obj.autosave);
            end

            obj.sj.traces.Fpx = {};          
            
            obj.sj = obj.sj.update_currentstate('Calcium traces extracted');
            obj.sj.save2mat(obj.autosave)

        end

        function obj = TracesQC(obj, ignore_previous)
            arguments
                obj 
                ignore_previous logical = false 
            end            

            obj.sj.traces = obj.sj.traces.setManually(ignore_previous);

            % expand with quality metrics (make use of TraceViewer)
        end

    end
end