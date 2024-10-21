%%% Preprocessing workflow

%% Initialize variables

clear all

%% startup preprocessing
p = Preprocessing;

%optional
p.autosave = true;

%% create or load fish (autosaving never applies to subject creation)
subjectID = 'TC_230910_TC0028_230906beh1b3_sxpDp_odorexp005_RPB3144501500AG';
group = 'naive';
datafolder = 'scratch\gfriedri\caudtomm\testground';
odordelay = 0;
p = p.createSubject(subjectID,group,datafolder,odordelay);
% or
p.sj = fish1; clear fish1

%% histogram equalization
p = p.claheRawTrials;


%% registration

p = p.allRegistration;

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



