%%% Preprocessing workflow

%% Initialize variables

clear all

% locations.drive = '\\tungsten-nas.fmi.ch\tungsten';
locations.drive = 'W:';

locations.subject_ID = 'TC_230519_TC0003_230516beh2b2_sxpDp_odorexp004_RPB3144501500AG';
locations.general_datapath = fullfile(locations.drive,"scratch\gfriedri\caudtomm\ev_data");
locations.subject_datapath = fullfile(locations.general_datapath,locations.subject_ID);
locations.datafile_ext = 'tif';

locations.rawtrials = 'trials';
locations.histeqtrials = 'trials_clahe';
locations.rawtrials_rigidreg = 'reg_stacks3_raw';
locations.histeqtrials_rigidreg = 'reg_stacks3_clahe';

locations.warp.rawtrials_opticflow = 'trials_warp_of';
locations.warp.rawtrials_opticflow_unchecked = 'flow_registration_results';
locations.warp.rawtrials_bunwarpj = 'trials_warp_bunwarpj';

locations.references.raw = fullfile('references','rawtrials');
locations.references.histeq = fullfile('references','histeqtrials');

locations.defRois = 'defROIs';



% odor delay to nostril
odordelay = 0; % sec

%% make new fish

% startup
cd(locations.subject_datapath)

% Create a new Subject object (this is a new fish)
fish1 = Subject('fish1',locations, 'trained2');

% set known attributes
fish1.odor_delay = odordelay;

%% rigid frame-by-frame registration
[fish1,~,~,stats] = fish1.registration('rigid');
anatomy = fish1.retrieve_trial_anatomies(fish1.locations.rawtrials_rigidreg);
fish1.anatomy_imgs = anatomy;
fish1.save2mat()

%% unwarping 1
fish1 = fish1.registration('opticflow');
% or run StackDewobbler.m


%% unwarping 2

ReplaceBadWarpPeriods; % NEEDS TESTING

%% ROI SELECTION
input_folder = fish1.locations.rawtrials_rigidreg;
fish1 = fish1.define_ROIs(input_folder,'manual')
fish1 = fish1.reload
%%
fish1.define_ROIs(input_folder,'automatic')
fish1 = fish1.update_currentstate('ROI selection complete');
fish1.save2mat(true)

%%
fish1.calculate_dFoverF
fish1 = fish1.update_currentstate('dF/F calculation complete');
fish1.save2mat(true)

%%

fish1.extract_data
fish1 = fish1.reload
fish1 = fish1.update_currentstate('Generated DATA Structure');
fish1.save2mat(true)
fish1.loadDATAfile








