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

p = p.rigidregHisteq2Raw;
% or
p = p.rigidregRaw;
% or
p = p.opticflowregRaw;
% or (invalid - doesn't touch raw trials)
p = p.opticflowregHisteq;

p = p.updateSubject(p.sj.locations.rawtrials_opticflowwarp_fromhisteq);

%% traces

p = p.selectROIs;

p.sj = p.sj.setManually;

p = p.extractCalciumTraces;


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






