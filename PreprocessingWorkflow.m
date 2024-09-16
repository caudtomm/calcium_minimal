%%% Preprocessing workflow

%% Initialize variables

clear all

subjectID = 'TC_230910_TC0028_230906beh1b3_sxpDp_odorexp005_RPB3144501500AG';
locations = Locations;
locations = locations.setSubjectID(subjectID);
locations = locations.setDataFolder('scratch\gfriedri\caudtomm\testground');

% odor delay to nostril
odordelay = 0; % sec

%% make new fish

% startup
cd(locations.subject_datapath)

% Create a new Subject object (this is a new fish)
fish1 = Subject('fish1',locations, 'naive');

% set known attributes
fish1.odor_delay = odordelay;
fish1.notes = ""; % load notes file

cd(fish1.locations.subject_datapath)

%% clahe

datapath = fish1.locations.rawtrials;

% initialize Batch Process
b = BatchProcess(BasicMovieProcessor('clahe'));
b.init.description = "CLAHE raw trials";
b = b.setDataPath(datapath);


% run
b = b.run;

% move to layer 1 dir
sourceFolder = fullfile(datapath, b.OutFolder);
destinationFolder = fish1.locations.histeqtrials;
movedirTC(sourceFolder,destinationFolder)

% save a separate clahe reference
cd(destinationFolder)
ref = fish1.retrieveSnippetfromCoords(fullfile(fish1.locations.subject_datapath, ...
    fish1.locations.references.raw,'metadata.json'));
ref.metadata.Name = "Raw Reference Image (CLAHE)";
foutname = getFileNameSpecs(fish1.filelist(ref.metadata.Trial_relativenum).name).fname; % w/o extension
fout = fullfile(fish1.locations.subject_datapath, ...
    fish1.locations.references.histeq,[foutname,'.mat']);
movie = ref; clear ref
save(fout,'movie','-mat'); clear movie
cd(fish1.locations.subject_datapath)


%% rigid frame-by-frame registration

% prep CLAHEd trials
datapath = fish1.locations.histeqtrials;
fish1.Tiff2Movie(datapath)

% initialize Batch Process
b = BatchProcess(RigidRegistration);
b.init.description = "Rigid alignment of CLAHEd raw trials to the CLAHEd reference image";
b.includeFilter = [num2str(fish1.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
b = b.setDataPath(datapath);
imgref = fish1.loadReferenceImg(fish1.locations.references.histeq);
b.Processor.reference_img=imgref;

% run
b = b.run;
charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
fish1.log{end+1} = charv;

% move to layer 1 dir
sourceFolder = fullfile(datapath, b.OutFolder);
destinationFolder = fish1.locations.histeqtrials_rigidreg;
movedirTC(sourceFolder,destinationFolder)


% apply computed shifts to the raw data
datapath = fish1.locations.rawtrials;
b.init.detailed_description = "Result shifts applied to raw trials";
b = b.applyresults(datapath);

% move to layer 1 dir
sourceFolder = fullfile(datapath, b.OutFolder);
destinationFolder = fish1.locations.rawtrials_rigidreg;
movedirTC(sourceFolder,destinationFolder)


fish1 = fish1.update_currentstate('Rigid alignment complete');
cd(fish1.locations.subject_datapath)


%% rigid frame-by-frame registration (directly on raw trials)

% prep CLAHEd trials
datapath = fish1.locations.rawtrials;
fish1.Tiff2Movie(datapath)

% initialize Batch Process
b = BatchProcess(RigidRegistration);
b.init.description = "Rigid alignment of raw trials to the raw reference image";
b.includeFilter = [num2str(fish1.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
b = b.setDataPath(datapath);
imgref = fish1.reference_img;
b.Processor.reference_img=imgref;

% run
b = b.run;
charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
fish1.log{end+1} = charv;

% % move to layer 1 dir
% sourceFolder = fullfile(datapath, b.OutFolder);
% destinationFolder = fish1.locations.rawtrials_rigidreg;
% movedirTC(sourceFolder,destinationFolder)

fish1 = fish1.update_currentstate('Rigid alignment complete');
fish1.save2mat()
cd(fish1.locations.subject_datapath)


%% remove baseline

% datapath = fish1.locations.rawtrials_rigidreg;
% 
% % initialize Batch Process
% b = BatchProcess(BasicMovieProcessor('subtract_baseline'));
% b.init.description = "Removing baseline from registered raw trials";
% b.includeFilter = [num2str(fish1.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
% b = b.setDataPath(datapath);
% 
% % run
% b = b.run;
% charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
% fish1.log{end+1} = charv;
% 
% % move to layer 1 dir
% sourceFolder = fullfile(datapath, b.OutFolder);
% destinationFolder = [fish1.locations.rawtrials_rigidreg,'_rmbase'];
% movedirTC(sourceFolder,destinationFolder)
% fish1.locations.rawtrials_rigidreg = destinationFolder;
% 
% 
% fish1 = fish1.update_currentstate('Movie baseline removed');
% cd(fish1.locations.subject_datapath)


%% update Subject
fish1.anatomy_imgs = fish1.retrieve_trial_anatomies(fish1.locations.rawtrials_rigidreg);
fish1.reference_img = fish1.retrieve_ref_img(fish1.locations.rawtrials_rigidreg);
fish1.localcorr_imgs = fish1.retrieve_localcorr_maps(fish1.locations.rawtrials_rigidreg);

fish1 = fish1.update_currentstate( ...
    'Trial anatomies and reference frame updated in Subject file');
fish1.save2mat()

fish1.visualize_anatomy_physFOV
cd(fish1.locations.subject_datapath)

%% opticflow registration

% prep raw trials
datapath = fish1.locations.rawtrials;
fish1.Tiff2Movie(datapath)

% initialize Batch Process
b = BatchProcess(OpticFlowRegistration);
b.init.description = "OpticFlow alignment of raw trials to the rigidreg reference image";
b.includeFilter = [num2str(fish1.min_trialnum, '%05.f'),'.mat']; % include only actual trials, not the last batch run result MAT
b = b.setDataPath(datapath);
imgref = fish1.reference_img;
b.Processor.reference_img=imgref;

% run
b = b.run;
charv = [char(datetime('today')), ' - hash:', b.hash, ' - ',b.init.description];
fish1.log{end+1} = charv;

% move to layer 1 dir
sourceFolder = fullfile(datapath, b.OutFolder);
destinationFolder = fish1.locations.rawtrials_opticflowwarp;
movedirTC(sourceFolder,destinationFolder)


fish1 = fish1.update_currentstate('Opticflow alignment complete');
fish1.save2mat()
cd(fish1.locations.subject_datapath)


%% ROI Selection

% select background ROIs
fish1.backgroundROImap = imageSequenceGUI(...
    fish1.anatomy_imgs,...
    fish1.localcorr_imgs,...
    fish1.backgroundROImap,...
    'Select background ROIs (N->E->S->W)');

% select cell ROIs
if exist("plane",'var'); fish1.ROImap = plane{1}.ROI_map; end
fish1.ROImap = imageSequenceGUI(...
    fish1.anatomy_imgs,...
    fish1.localcorr_imgs,...
    fish1.ROImap,...
    'Select neuronal somata');

fish1.save2mat()

%%
fish1 = fish1.update_currentstate( ...
    'ROI selection complete');
fish1.save2mat()

%% Extraction of calcium traces

% address to find the source movies
fish1.locations.traces_src = fish1.locations.rawtrials_rigidreg;

loc = fish1.locations;

PMToffmeta = fullfile(loc.subject_datapath, loc.rawtrials, 'PMToff_metadata.json');
noLightmeta = fullfile(loc.subject_datapath, loc.rawtrials, 'noLight_metadata.json');
fish1 = fish1.defineTraces(PMToffmeta,noLightmeta);
cd(fish1.locations.subject_datapath)

refmovieforbackground = 'TC_240218_TC0028_240213beh1b2_sxpDp_odorexp004_RPB3144501500AG_00035_00001.tif';
refmovieforbackground = fullfile(fish1.locations.general_datapath,refmovieforbackground);


cd(fish1.locations.traces_src)

interval = 960:1100;
fish1.traces = fish1.traces.setPMToff(Snippet(refmovieforbackground,interval));
interval = 1160:1300;
fish1.traces = fish1.traces.setNoLight(Snippet(refmovieforbackground,interval));
fish1.traces = fish1.traces.defineFundamentalProperties(fish1);
fish1.traces = fish1.traces.setDerivativeProperties;

cd(fish1.locations.subject_datapath)


fish1 = fish1.update_currentstate('Calcium traces extracted');
fish1.save2mat()

%% unwarping 1
% fish1 = fish1.registration('opticflow');
% or run StackDewobbler.m

% 
% %%
% 
% fish1.extract_data
% fish1 = fish1.reload
% fish1 = fish1.update_currentstate('Generated DATA Structure');
% fish1.save2mat(true)
% fish1.loadDATAfile
% 
% 






