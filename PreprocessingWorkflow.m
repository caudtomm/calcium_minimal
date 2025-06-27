%%% Preprocessing workflow

%% Initialize variables

clear all

%% startup preprocessing
p = Preprocessing;

%optional
p.autosave = true;

%% create or load fish (autosaving never applies to subject creation)
drive = '\\tachyon.fmi.ch\tachyon\';
datafolder = 'groups\scratch\gfriedri\processed-data\250129_TC_invivoCaIMG_odorexp004_005';

% --------
subjectID = 'TC_241213_TC0003_241209beh1B2_sxpDp_odorexp004_RPB3144501500AG';
group = 'trained1';
odordelay = 0;
p = p.createSubject(subjectID,group,drive,datafolder,odordelay);

% or

load('fish1.mat')
fish1.locations = fish1.locations.setDrive(drive);
fish1.locations = fish1.locations.setDataFolder(datafolder);
p.sj = fish1; clear fish1
% --------

%% early preprocessing

% remove bad periods and replace them with nan frames
p = p.replaceBadPeriodsWithNans;

% histogram equalization
p = p.claheRawTrials;

%% correct bidirectional scan artefact
p = p.correctBidiScanningHisteq2Raw;

%% registration

% FOR TESTING
% % can run anywhere
% p = p.allRigidReg;
% 
% % needs to run on windows (best use FAIM Workstations)
% p = p.allWarpReg;

% FOR DEPLOYMENT (run on Windows with sufficient RAM (128 GB is optimal))
p = p.opticflowHisteq2Raw;

% registration correction
p = p.correctRegistrationResults;

%% PCAICA
% to remove the most highly low-freq-biased IC (usually also the most variant in time)

p = p.subtractIC;

%% select best source of registration, after visual inspection of the results
src = p.sj.locations.rawtrials_opticflowwarp_fromhisteq_corrected_noIC;
p = p.updateSubject(src);

%% traces

p = p.selectROIs; % manual input required

p.sj = p.sj.setManually; % manual input required
p.sj.save2mat(p.autosave);

tracesfolder = 'traces\manIC_conservative';
p = p.extractCalciumTraces(tracesfolder,true);

p = p.TracesQC; % manual input required
p.sj.save2mat(p.autosave);

% save traces after QC
mode = 'light';
fname = ['traces',mode,'.mat'];
FileOut = fullfiletol(tracesfolder,fname);
p.sj.traces.save(FileOut,mode,p.autosave);

%% construct Experiment
clear all

thisdrive = '\\tachyon.fmi.ch\tachyon\';
datafolder = 'groups\scratch\gfriedri\processed-data\250129_TC_invivoCaIMG_odorexp004_005';

tracesfolder = 'traces\top1biasedIC';

a = Experiment(fullfiletol(thisdrive,'\groups\scratch\gfriedri\caudtomm\data_record2.xlsx'),tracesfolder,'odorexp004_analysis');
a.locations = a.locations.setDrive(thisdrive);
a.locations = a.locations.setDataFolder(datafolder);


% load 'light' traces (without single px values)
a = a.loadSubjectTraces;

a.name = "odorexp004_IC1_130625";
save(a.name,'a','-v7.3')

%% output backward compatible experiment structure and save to file

exp_name = char(strcat('exp_',a.name,'_iSR'));
experiment = a.convert2BackwardCompatibleStruct(exp_name);

% save to pwd
s.experiment = experiment;
robust_io('save',[exp_name,'.mat'],s,'-mat','-v7.3')

%%



