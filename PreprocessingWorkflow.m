%%% Preprocessing workflow

%% Initialize variables

clear all

%% startup preprocessing
p = Preprocessing;

%optional
p.autosave = true;

%% create or load fish (autosaving never applies to subject creation)
subjectID = 'TC_240217_TC0028_240213beh1b4_sxpDp_odorexp004_RPB3144501500AG';
group = 'uncoupled';
drive = '/tungstenfs/';
datafolder = 'scratch\gfriedri\caudtomm\testground';
odordelay = 0;
p = p.createSubject(subjectID,group,drive,datafolder,odordelay);
% or
p.sj = fish1; clear fish1

%% early preprocessing

% remove bad periods and replace them with nan frames
p = p.replaceBadPeriodsWithNans;

% histogram equalization
p = p.claheRawTrials;

%% correct bidirectional scan artefact
p = p.correctBidiScanningHisteq2Raw;



%% registration

% can run anywhere
p = p.allRigidReg;

% needs to run on windows (best use FAIM Workstations)
p = p.allWarpReg;

%% select best source of registration, after visual inspection of the results
src = p.sj.locations.rawtrials_opticflowwarp_fromhisteq
p = p.updateSubject(src);

%% traces

p = p.selectROIs; % manual input required

p.sj = p.sj.setManually; % manual input required
p.sj.save2mat(p.autosave);

p = p.extractCalciumTraces(false);

p = p.TracesQC; % manual input required

p.sj.traces.save('','full',p.autosave);
p.sj.traces.save('','light',p.autosave);

%% construct Experiment

thisdrive = '\\tungsten-nas.fmi.ch\tungsten';
a = Experiment(fullfiletol(thisdrive,'\scratch\gfriedri\caudtomm\data_record2.xlsx'),'odorexp004_analysis');
a.locations = a.locations.setDrive(thisdrive);
a.locations = a.locations.setDataFolder('scratch\gfriedri\caudtomm\testground'); % or whatever


% load 'light' traces (without single px values)
a = a.loadSubjectTraces;

%% output backward compatible experiment structure and save to file

exp_name = 'new_analysis';
experiment = a.convert2BackwardCompatibleStruct(exp_name);

% save to pwd
save([expname,'.mat'],'experiment','-mat','-v7.3')

%%



