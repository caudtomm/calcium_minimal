%% how many clusters/response? order by trial n

a = data.singletrial;
ntrials = numel(a);

nclust = [];

for i_trial = 1:ntrials
    thisnclust = [a{i_trial}.intime.resp_stage_clust{1}.optNClust, ...
        a{i_trial}.intime.resp_stage_clust{1}.optNClust2];
    
    nclust = [nclust; thisnclust];
end

sortorder = 1:ntrials;
hf = figure(896);
h = barh(nclust,'stacked'); axis square; box off
xlabel('number of response stages')
lab = {'optimal', 'second optimal'};
legend(h, lab)
yticks(1:size(nclust,1)); yticklabels(data.stim_type(sortorder))

% save fig
fnameout = 'n_reponse_stages';
FileOut = [fnameout, '.jpeg']
saveas(hf,FileOut,'jpeg')                       % JPEG
FileOut = [fnameout, '.fig']
savefig(hf,FileOut)                             % FIG




nclust = [];

for i_trial = 1:ntrials
    thisnclust = [a{i_trial}.intime.resp_stage_clust{1}.nclustcorrect, ...
        a{i_trial}.intime.resp_stage_clust{1}.nclustcorrect2];
    
    nclust = [nclust; thisnclust];
end

sortorder = 1:ntrials;
hf = figure(897);
h = barh(nclust,'stacked'); axis square; box off
xlabel('number of response stages')
lab = {'optimal', 'second optimal'};
legend(h, lab)
yticks(1:size(nclust,1)); yticklabels(data.stim_type(sortorder))

% save fig
fnameout = 'n_reponse_stages_corrected';
FileOut = [fnameout, '.jpeg']
saveas(hf,FileOut,'jpeg')                       % JPEG
FileOut = [fnameout, '.fig']
savefig(hf,FileOut)                             % FIG