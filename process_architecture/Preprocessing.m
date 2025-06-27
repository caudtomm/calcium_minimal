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
            % all the trials, for easy superficial visual inspection
            RegistrationViewer(obj.sj,destinationFolder). ...
                subsampledMovieCollage.save('','avi','sidebyside');
        end
    end

    methods
        function obj = Preprocessing()
            
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
    
        function obj = createSubject(obj,subjectID,group,drive,datafolder,odordelay)
            arguments
                obj 
                subjectID char
                group char
                drive char = 'W:'
                datafolder char = 'scratch\gfriedri\caudtomm\testground'
                odordelay double = 0 % odor delay to nostril [sec]
            end

            locations = Locations;
            locations = locations.setSubjectID(subjectID);
            locations = locations.setDataFolder(datafolder);
            locations = locations.setDrive(drive);
            
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
    
        function [obj, b] = replaceBadPeriodsWithNans(obj)
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            datapath = obj.sj.locations.orig_trials;
            
            % initialize Batch Process
            b = BatchProcess(BasicMovieProcessor('replace_badperiods_with_nans'));
            b.init.description = "Replace bad periods with nans in raw trials";
            b = b.setDataPath(datapath);            
            
            % run
            b = b.run;
            
            % move to layer 1 dir
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.rawtrials;
            movedirTC(sourceFolder,destinationFolder)

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
        end

        function [obj, b] = claheRawTrials(obj)
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            datapath = obj.sj.locations.rawtrials;
            
            % initialize Batch Process
            b = BatchProcess(BasicMovieProcessor('clahe'));
            b.init.description = "CLAHE raw trials";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
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
            foutpath = fullfiletol(obj.sj.locations.subject_datapath, ...
                obj.sj.locations.references.histeq);
            if ~exist(foutpath, 'dir'); mkdir(foutpath); end
            fout = fullfiletol(foutpath,[foutname,'.mat']);
            movie = ref; clear ref
            s.movie = movie;
            robust_io('save',fout,s,'-mat'); clear movie
            cd(fullfiletol(obj.sj.locations.subject_datapath))

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
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
    
        function [obj, b] = correctBidiScanningHisteq(obj)
            % prep histeq trial movies
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            datapath = obj.sj.locations.histeqtrials;
            obj.sj.Tiff2Movie(datapath)
            
            % initialize Batch Process
            b = BatchProcess(CorrectBidiScanning);
            b.init.description = "Correction of the bidirectional alignment artefact";
            b.init.detailed_description = "Based on the consensus transform across lines of the histeq reference image";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
            b = b.setDataPath(datapath); % only needed to get b.DataList - b.applyresults will set the new path anyway
            imgref = obj.sj.loadReferenceImg(obj.sj.locations.references.histeq);
            b.Processor = b.Processor.setReference(imgref);

            % precompute shifts and apply to everything
            comp = CorrectBidiScanning(Movie(imgref)).run(true);
            b.results = repmat({comp}, size(b.DataList));
            
            % run
            b = b.applyresults(datapath);
            charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
            obj.sj.log{end+1} = charv;
            
            % move to layer 1 dir
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.histeqtrials;
            movedirTC(sourceFolder,destinationFolder)

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);

            cd(fullfiletol(obj.sj.locations.subject_datapath))
        end

        function [obj, b] = correctBidiScanningHisteq2Raw(obj)
            [obj, b] = correctBidiScanningHisteq(obj);
            
            % apply computed shifts to the raw data
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            datapath = obj.sj.locations.rawtrials;
            
            % run
            b = b.applyresults(datapath);
            charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
            obj.sj.log{end+1} = charv;
            
            % move to layer 1 dir
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.rawtrials;
            movedirTC(sourceFolder,destinationFolder)

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
            
            obj.sj = obj.sj.update_currentstate('Bidi correction complete');
            obj.sj.save2mat(obj.autosave)
            cd(fullfiletol(obj.sj.locations.subject_datapath))
        end

        % _________________________________________________
        % FRAME-BY-FRAME REGISTRATION

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
            obj = obj.allRigidReg;
            % then
            obj = obj.allWarpReg;
        end

        function obj = allRigidReg(obj)
            % every rigid registration protocol, done sequentially
            obj = obj.rigidregRaw;
            % then
            obj = obj.rigidregHisteq2Raw;
        end

        function obj = allWarpReg(obj)
            % every warp registration protocol, done sequentially
            obj = obj.opticflowregRaw;
            % then
            obj = obj.opticflowHisteq2Raw;
        end

        % _________________________________________________
        
        function [obj, b] = correctRegistrationResults(obj)
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            datapath = obj.sj.locations.rawtrials_opticflowwarp_fromhisteq; % < path to pre-egistered movies
            sourcepath = obj.sj.locations.rawtrials; % < path to raw movies

            % check for consistent filenum across reference and source
            % folders
            [fileList, numFiles] = obj.sj.getFileListInSubfolder(datapath);
            [~, numFiles2] = obj.sj.getFileListInSubfolder(sourcepath);
            if numFiles~=obj.sj.getNTrials; warning( ...
                    'Number of movies found  in the destination folder is different from the number of trials.'); end
            if numFiles2~=obj.sj.getNTrials; warning( ...
                    'Number of movies found  in the source folder is different from the number of trials.'); end

            % quality control of registration results in datapath
            th = -5;
            v = RegistrationViewer(obj.sj,datapath);
            [~,~,putative_fails,hf] = v.QC(th); % identifying putatively failed frames
            savefig(hf, fullfile(datapath,'pre-correction_qc.fig'));

            % load batchprocess json (to have a list of init subhashes for the transformation outputs of each file)
            bias = 0; if contains(datapath,'2'); bias = numFiles; end
            filename = dir(fullfile(datapath,'*_init.json'));
            batchinit = readJson(fullfile(filename(1).folder,filename(1).name));
            initlist = batchinit.subHash(bias+(1:numFiles));
            
            % build 'results' array for operation
            results = cell(numFiles,1);
            for i = 1:numFiles
                op.fail_idx = putative_fails(:,i);
                op.preregmovie = fileList{i};
                op.transforminit = fullfile(datapath,'inits',[initlist{i},'.mat']);
                op.zscorethresh = th;

                results{i}.operation = op;

                results{i}.init = readJson(fullfile(datapath,'inits',[initlist{i},'_init.json']));
            end

            % initialize Batch Process
            b = BatchProcess(RegistrationCorrection);
            b.init.description = "Corrected registration results";
            b.init.detailed_description = "re-registered failed frames + added uncorrectable frames to badperiods";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT           
            b.results = results;

            % run
            b = b.applyresults(sourcepath);
            
            % move to layer 1 dir
            sourceFolder = fullfiletol(sourcepath, b.OutFolder);
            destinationFolder = obj.sj.locations.rawtrials_opticflowwarp_fromhisteq_corrected;
            movedirTC(sourceFolder,destinationFolder)
            
            % new QC on the output
            v = RegistrationViewer(obj.sj,destinationFolder);
            [~,~,~,hf] = v.QC(th);
            savefig(hf, fullfile(destinationFolder,'post-correction_qc.fig'));

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
        end

        function obj = subtractIC(obj, idx)
            arguments
                obj 
                idx double = []
            end
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            datapath = obj.sj.locations.rawtrials_opticflowwarp_fromhisteq_corrected;

            %% check for pre-existing PCAICA output files
            FileIn{1} = 'pcaica_w.mat';
            FileIn{2} = 'pcaica.mat';
            if exist(FileIn{1},"file") && exist(FileIn{2},"file")
                pcaica_w = robust_io('load', FileIn{1}).pcaica_w;
                decomposition = robust_io('load', FileIn{2}).pcaica;
            else
                movie = PCAICA.getFull2PMovieSubsampled(obj.sj,datapath);
                decomposition = PCAICA(movie); 
                [decomposition, ~,~, pcaica_w, ~,~] = decomposition.run();

                % save figures (in a really stupid way, sorry)
                mkdir('figures')
                figs = findall(0, 'Type', 'figure');
                for i = 1:numel(figs)
                    fig = figs(i);
                    fname = sprintf('fig_%d.fig', i); % or use i instead of fig.Number
                    try; saveas(fig, fullfiletol('figures',fname), 'fig'); catch; end
                end
                close all
            end


            %%

            if isempty(idx)
                % select idx of ICs to reconstruct in isolation, so they can be
                % removed from the original data
                [~,idx] = sort(pcaica_w.spectral_bias.bias_scores,'descend');
                idx = idx(1); % isolate only the top low-freq-biased IC
            end
            disp(['Removing following ICs: ',mat2str(idx)]) % print IC #

            % resize weight maps to original size
            resize_factor = decomposition.operation.resize_factor;
            wmapsIC = pcaica_w.ICweightmaps.whole;
            wmapsIC = cellfun(@(x) imresize(x,'nearest','scale',resize_factor), ...
                wmapsIC,'UniformOutput',false);

            % re-linearize wmaps and eliminate pre-identified NaNs
            [wmapsIC, idx_nanfr, idx_nanpx] = cellfun(@(x) ...
                PCAICA.linearizeFrames(x, true,[],[],decomposition.init.method), ...
                wmapsIC, 'UniformOutput',false);
            idx_nanfr = idx_nanfr{1};
            idx_nanpx = idx_nanpx{1};

            % concatenate linearized maps to one 2d array
            wmapsIC = cell2mat(wmapsIC'); % [px x IC#]

            % useful later
            pxstd = pcaica_w.pxstdfullsize';
            pxstd(idx_nanpx) = [];
            pxmu = pcaica_w.pxmufullsize';
            pxmu(idx_nanpx) = [];

            % keep only idx
            to_elim = true(width(wmapsIC),1); to_elim(idx) = false;
            Wnew = wmapsIC; Wnew(:,to_elim) = zeros(height(wmapsIC),sum(to_elim));
            W = wmapsIC * pinv(Wnew);

            % -------

            % package process parameters
            op.pxmu = pxmu;
            op.pxstd = pxstd;
            op.transform = W;
            op.idx_nanpx = idx_nanpx;
            op.inputmethod = decomposition.init.method;
            op.mode = 'nooffset';
            op.subtract_from_inputdata = true;

            % single file run (test)
            %filename = 'OFreg_stacks_clahe2raw/TC_241213_TC0003_241209beh1B2_sxpDp_odorexp004_RPB3144501500AG_00002_00001.mat';
            %reconstructedIC = ComponentReconstruction(op,filename).run().data_processed;

            % initialize Batch Process
            b = BatchProcess(ComponentReconstruction(op));
            b.init.description = "Subtract top bias-scored IC from fully registered trials";
            b.includeFilter = [num2str(obj.sj.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
            b = b.setDataPath(datapath);            
            
            % run
            b = b.run;
            
            %% move to layer 1 dir
            sourceFolder = fullfiletol(datapath, b.OutFolder);
            destinationFolder = obj.sj.locations.rawtrials_opticflowwarp_fromhisteq_corrected_noIC;
            movedirTC(sourceFolder,destinationFolder)

            % save anatomy and visualization results to disk
            obj = obj.tail_sequence(destinationFolder);
        end
        
        function obj = selectROIs(obj, ROImap)
            arguments
                obj 
                ROImap double = [] 
            end
            % select background ROIs
            obj.sj.backgroundROImap = imageSequenceGUI(...
                obj.sj.anatomy_imgs,...
                obj.sj.localcorr_imgs,...
                obj.sj.backgroundROImap,...
                'Select background ROIs (N->E->S->W)');
            obj.sj.backgroundROImap = obj.pruneROIs(obj.sj.backgroundROImap);
            
            % select cell ROIs
            if isempty(ROImap); ROImap = obj.sj.ROImap; end
            obj.sj.ROImap = imageSequenceGUI(...
                obj.sj.anatomy_imgs,...
                obj.sj.localcorr_imgs,...
                ROImap,...
                'Select neuronal somata');
            obj.sj.ROImap = obj.pruneROIs(obj.sj.ROImap);
            
            % save
            obj.sj = obj.sj.update_currentstate( ...
                'ROI selection complete');
            obj.sj.save2mat(obj.autosave)
        end

        function obj = extractCalciumTraces(obj, outfolder, do_save)    
            arguments
                obj 
                outfolder char = pwd
                do_save logical = true
            end
            loc = obj.sj.locations;

            % temporary
            obj.sj.backgroundROImap = obj.pruneROIs(obj.sj.backgroundROImap);
            obj.sj.ROImap = obj.pruneROIs(obj.sj.ROImap);
            
            obj.sj = obj.sj.defineTraces; 
            cd(fullfiletol(obj.sj.locations.subject_datapath))
            
            defaultrefmovieforbackground = 'TC_240218_TC0028_240213beh1b2_sxpDp_odorexp004_RPB3144501500AG_00035_00001.tif';
            defaultrefmovieforbackground = fullfiletol(obj.sj.locations.general_datapath,defaultrefmovieforbackground);
            
            cd(fullfiletol(obj.sj.locations.traces_src))

            % set PMToff
            PMToffmeta = fullfiletol(loc.subject_datapath, loc.orig_trials, 'PMToff_metadata.json');
            PMToffmeta = readJson(PMToffmeta);
            interval = PMToffmeta.Frame_range(1):PMToffmeta.Frame_range(2);
            if PMToffmeta.Internal_reference
                refmovie = obj.sj.filelist(PMToffmeta.Trial_relativenum+obj.sj.min_trialnum-1).name;
                refmovie = fullfiletol(loc.subject_datapath,loc.orig_trials,refmovie);
                obj.sj.traces = obj.sj.traces.setPMToff(Snippet(refmovie,interval));
            else
                obj.sj.traces = obj.sj.traces.setPMToff(Snippet(defaultrefmovieforbackground,interval));
            end

            % set noLight
            noLightmeta = fullfiletol(loc.subject_datapath, loc.orig_trials, 'noLight_metadata.json');
            noLightmeta = readJson(noLightmeta);
            interval = noLightmeta.Frame_range(1):noLightmeta.Frame_range(2);
            if noLightmeta.Internal_reference
                refmovie = obj.sj.filelist(noLightmeta.Trial_relativenum+obj.sj.min_trialnum-1).name;
                refmovie = fullfiletol(loc.subject_datapath,loc.orig_trials,refmovie);
                obj.sj.traces = obj.sj.traces.setNoLight(Snippet(refmovie,interval));
            else
                obj.sj.traces = obj.sj.traces.setNoLight(Snippet(defaultrefmovieforbackground,interval));
            end

            % get the actual traces
            obj.sj.traces = obj.sj.traces.defineFundamentalProperties(obj.sj);
            obj.sj.traces = obj.sj.traces.setDerivativeProperties;
            
            cd(fullfiletol(obj.sj.locations.subject_datapath))

            if do_save
                mode = 'full';
                fname = ['traces',mode,'.mat'];
                FileOut = fullfiletol(outfolder,fname);
                obj.sj.traces.save(FileOut,mode,obj.autosave);

                mode = 'light';
                fname = ['traces',mode,'.mat'];
                FileOut = fullfiletol(outfolder,fname);
                obj.sj.traces.save(FileOut,mode,obj.autosave);

                obj.sj.traces.Fpx = {}; 
            end                                 
            
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