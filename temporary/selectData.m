function [a, X, stims,traces,reltrialn, trialidx] = selectData(experiment, todo_fish, stims2use, trialn2usein, observations)

%% initialize global vars
global ntrials L ncells labs s C pathout tframe fs trialn2use stims

%%

% pathout = 'figures';
% if ~isfolder(pathout); mkdir(pathout); end
% cd(pathout)

% all stimuli
% stims = {'spont.','ACSF','Trp','Ser','Ala','Food'};
% stims = {'spont.','Arg','Ala','His','Trp','Ser','ACSF','Leu'};
stims = {'Arg','Ala','His','Trp','Ser','Leu'};
% stimuli we want to include
% stims2use = {'Trp'};
s2u_idx = find(ismember(stims,stims2use));

% trial nums to include
trialn2use = trialn2usein;

% select observations : "trialavg" | "time"
if ~exist('observations','var'); observations = "time"; end

% % time period to include (relative to [stim_on, stim_on) [s]
% tframe = [.5, 20];

% save?
s = true;

% tag for saving (creates a new folder in cd\selectData)
tag = strcat(cell2mat(stims2use),'_trial',mat2str(trialn2use),'_fish',mat2str(todo_fish));


%% initialize params


%
data = experiment.series{todo_fish}.data;
fs = data.meta.framerate;
stim_on = data.stim_on_sec;
stim_off = data.stim_off_sec;
labs = stims;

traces = data.tracesdesdn;
% FindTopUnits
% traces = selectCells(traces,data.L,alltopunits');
traces = traceFormat(traces,data.L);
interval = floor((stim_on + tframe(1))*fs):floor((stim_off + tframe(2))*fs);
traces = traces(interval,:,:);
T = diff([stim_on, stim_off]) + diff(tframe);
% traces(traces<-.2) = -.2;

% stimulus attribution
odortrial_idx = [];
for i_stim = 1:numel(stims)
    odortrial_idx(i_stim,:) = ismember(data.stim_type,stims{i_stim});
end
reltrialn = max(cumsum(odortrial_idx')'.*odortrial_idx);
[~,X] = max(odortrial_idx);
C = [1 1 1;
     1 1 0;
     0 0 1;
     1 0 0;
     0 1 1;
     0 1 0;
     1 0 1;
     .5 .8 .2];

% save folder
pathout = fullfile(mfilename,tag)
if ~exist(pathout,'dir'); mkdir(pathout); end

%% select data

% select for odors
idx1 = ~ismember(data.stim_type,stims2use);

% select for relative trial number
idx2 = ~ismember(reltrialn,trialn2use);

idx = idx1 | idx2;
traces(:,:,idx) = [];
X(idx) = [];
reltrialn(idx) = [];


ntrials = size(traces,3);
L = size(traces,1);
ncells = size(traces,2);

trialidx = find(~idx);

%% define dataset

a = [];
switch observations
    case "trialavg"
        a = squeeze(nanmean(traces,1))'; whos a
        Xrep = X;
    case "time"
        a = nan(size(traces,1)*size(traces,3),size(traces,2));
        for i_trial = 1:size(traces,3)
            idx = (i_trial-1) * size(traces,1) + [1:size(traces,1)];
            a(idx,:) = traces(:,:,i_trial);
        end; whos a
        Xrep = repelem(X,size(traces,1));
    otherwise
        error('Observation format not recognised.')
end

