% all stimuli
% stims = {'spont.','ACSF','Trp','Ser','Ala','Food'};
% stims2useIDX = {[3:6]; ...
%                 [3]; ...
%                 [4]; ...
%                 [5]; ...
%                 [3,5]};
% 
% trialn2useM = {[1:12]; ...
%                [2:12]};

% stims = {'spont.','Arg','Ala','His','Trp','Ser','ACSF','Leu','Food'};
stims = {'Arg','Ala','His','Trp','Ser','Leu'};
% stims2useIDX = {[1:9]};
stims2useIDX = {[1:6]};

trialn2useM = {[1:5]; ...
               [2:5]};

%%

for i_fish = 1:numel(experiment.series)
    for i_stim = 1:numel(stims2useIDX)
        stims2use = stims(stims2useIDX{i_stim});
        for i_trials = 1:numel(trialn2useM)
            trialn2use = trialn2useM{i_trials};
            UMAPandPCA(experiment,i_fish, stims2use, trialn2use,'time') % or 'trialavg'
            close all
        end
    end
end