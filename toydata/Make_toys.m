p = ToyParams;

p.D = 30;
p.A = 1;
p.lambda_A = 6;
p.tau_A = 1.5;
p.B = 9;
p.tau = 1;
p.gamma = 0;
p.C = 2.5;
p.theta = .06;
p.alpha = 0;
p.rho = .1;
p.rho_e = .2;
p.eta = 0;
p.rho_d = 0;

p.rotation_mix = .7;
p.fano_factor = 1.2;
p.fr_mean = .07;
p.fr_cv = .9;
p.noise_corr = .1;

p.stim_idx = [1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 6 1 2 3 4 5 6 6 6 6];
p.stimulus_names = {'Arg','Ala','His','Trp','Ser','Leu'};

%% plots patterns onto provided axes

tg = ToyDataGenerator(p);
[labs, idx, mus, fr_avg] = getPatterns(tg);

plotPatterns(labs, idx, mus, fr_avg)

% helper functions

function [labs, idx, mus, fr_avg] = getPatterns(tg)
arguments
    tg ToyDataGenerator
end

p = tg.params;

% stimulus window
t = p.odor_frames(1):p.odor_frames(2);
L = length(t);

% labels and indices
ntrials = p.K*p.R;
labs = p.stimulus_names(p.stim_idx);
[~,idx] = sort(labs);

% central patterns
mus = reshape(tg.cp.mu,p.N,ntrials);

% average traces
fr = reshape(tg.nm.firing_rates(:,t,:,:),p.N,L,ntrials);
fr_avg = squeeze(mean(fr,2));

end

function plotPatterns(labs, idx, mus, fr_avg)

% central patterns
subplot(141); imagesc(mus)

% correlations among central patterns
corrs = 1-squareform(pdist(mus',"correlation"));
subplot(342); imagesc(corrs); axis square; % chronological
subplot(346); imagesc(corrs(idx,idx)); axis square; % by stim

% euclidean distance among central patterns
corrs = squareform(pdist(mus',"euclidean"));
subplot(3,4,10); imagesc(corrs(idx,idx)); axis square; % by stim


% firing rates
subplot(143); imagesc(fr_avg)

% correlations among firing rates
corrs = 1-squareform(pdist(fr_avg',"correlation"));
subplot(344); imagesc(corrs); axis square; % chronological
subplot(348); imagesc(corrs(idx,idx)); axis square; % by stim

% euclidean distance among firing rates
corrs = squareform(pdist(fr_avg',"euclidean"));
subplot(3,4,12); imagesc(corrs(idx,idx)); axis square; % by stim

set(gcf,'Color','w')
end

%% batch generate toys

outdir = 'toys\hyp';
rmdir(outdir,'s')

ToyDataGenerator.generate_batch(p,outdir,5)
% ToyDataGenerator.generate_batch(p.makeNull,'toys\null',5)

experiment = ToyDataLoader.fromDir(outdir, 'trained1');

% load experiment into viewer
v = ExperimentViewer(experiment);
v.plotConfig = cfg;
v.dataFilter.traceType = 'pSpike';
v.dataFilter.interval = [1,20];
v.dataFilter.trial_sorting = 'stim_id';
dft = v.dataFilter;

% plot corrs
% figure; v.plotDistancesHead;

% plot PCA
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [.5, 20];
v.dataFilter.repetitions = 1:5;
[~,events,labs] = ModeSelector(v).extract;
proj = computeLDE(events,labs,'pooldata',true,'nans2zeros',true, 'method','pca');
figure; out = plotLDE(proj.embedding{1}.reduction(:,1:2,:),'lines',proj.labs{1},cfg, 'ldetype',proj.name); % plot
xlabel('PC 1'); ylabel('PC 2');

% plot UMAP
v.dataFilter.stims_allowed = 'all stimuli';
v.dataFilter.interval = [1, 20];
v.dataFilter.repetitions = 1:5;
[~,events,labs] = ModeSelector(v).extract;
proj = computeLDE(events,labs,'pooldata',true,'nans2zeros',true, ...
    'method','umap','n_components',2, 'metric','euclidean');
figure; out = plotLDE(proj.embedding{1}.reduction,...
    'lines',proj.labs{1},cfg, 'ldetype',proj.name); % plot
xlabel('UMAP 1'); ylabel('UMAP 2');