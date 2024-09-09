%%% Preprocessing workflow

%% Initialize variables

clear all

subjectID = 'TC_230809_TC0028_naive01_sxpDp_odorexp004_RPB3144501500AG';
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
% 
% % move to layer 1 dir
% sourceFolder = fullfile(datapath, b.OutFolder);
% destinationFolder = [fish1.locations.rawtrials_rigidreg,'_rmbase'];
% movedirTC(sourceFolder,destinationFolder)
% fish1.locations.rawtrials_rigidreg = destinationFolder;
% 
% 
% fish1 = fish1.update_currentstate('Movie baseline removed');


%% update Subject
fish1.anatomy_imgs = fish1.retrieve_trial_anatomies(fish1.locations.rawtrials_rigidreg);
fish1.reference_img = fish1.retrieve_ref_img(fish1.locations.rawtrials_rigidreg);
fish1.localcorr_imgs = fish1.retrieve_localcorr_maps(fish1.locations.rawtrials_rigidreg);

fish1.save2mat()
fish1 = fish1.update_currentstate( ...
    'Trial anatomies and reference frame updated in Subject file');

fish1.visualize_anatomy_physFOV

%% ROI Selection

% select background ROIs
fish1.backgroundROImap = imageSequenceGUI(...
    fish1.anatomy_imgs,...
    fish1.localcorr_imgs,...
    fish1.backgroundROImap,...
    'Select background ROIs (N->E->S->W)');

% select cell ROIs
fish1.backgroundROImap = imageSequenceGUI(...
    fish1.anatomy_imgs,...
    fish1.localcorr_imgs,...
    fish1.ROImap,...
    'Select neuronal somata');

fish1.save2mat()

%%
fish1 = fish1.update_currentstate( ...
    'ROI selection complete');

%% Extraction of calcium traces

% address to find the source movies
fish1.locations.traces_src = fish1.locations.rawtrials_rigidreg;

loc = fish1.locations;

PMToffmeta = fullfile(loc.subject_datapath, loc.rawtrials, 'PMToff_metadata.json');
noLightmeta = fullfile(loc.subject_datapath, loc.rawtrials, 'noLight_metadata.json');
fish1 = fish1.defineTraces(PMToffmeta,noLightmeta);

fish1 = fish1.update_currentstate('Calcium traces extracted');

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






