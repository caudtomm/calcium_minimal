function R = driftMetricsFigure(v, varargin)
%DRIFTMETRICSFIGURE Compute and visualize repetition drift metrics.
%
% This wraps your analysis into a single call that:
%   1) pulls tuning-curve data via v.plotUnitActivityMetricHead
%   2) computes several drift-related metrics
%   3) optionally renders a set of figures
%
% INPUT
%   v                struct/object with a .dataFilter field and a method:
%                    out = v.plotUnitActivityMetricHead('method','tuning curves')
%
% Name-Value Pairs (all optional)
%   'plot'           logical, make figures (default true)
%
% OUTPUT (struct R)
%   .out             cell, per subject data; each cell: [units x stims x reps]
%   .scalingf        [units_all x 1], mean over stims & reps per unit (omit NaNs)
%   .driftv          cell, per subject drift tensors: diff(out,[],3)
%   .corrMat         [stim_pairs_all x (reps-1)], 1 - correlation distance
%   .driftStrength   [units_all*stims x (reps-1)], normalized drift magnitudes
%   .interStimdist   [stim_pairs_all x reps], Euclidean distances per repetition
%   .interRepAngle   [triplets_all x 1], angle between successive drift steps
%   .anglecorr       [pairs_across_stims x 1], cosine similarity of angles
%   .edges           vector, histogram bin edges used for CDF plots
%   .signrank_p      scalar, p-value comparing corrMat(:,1) vs mean(:,2:end)
%   .figures         struct of figure handles (if plot==true)
%
% Example
%   R = driftMetricsFigure(v,'repetitions',1:5,'stimname','all novel','cfg',cfg);

% -------- parse args
p = inputParser;
addParameter(p,'plot',true);
parse(p,varargin{:});
o = p.Results;
cfg = v.plotConfig;
dft = v.dataFilter;

% -------- load data: out{sj} -> [units x stims x reps]
out = v.plotUnitActivityMetricHead('method','tuning curves');

% -------- per-unit scaling (mean across stims & reps, NaN-safe)
sf = cellfun(@(x) mean(x,[2,3],"omitmissing"), out, 'UniformOutput', false);
sf = cell2mat(sf);

% -------- per-subject drift tensors (difference across adjacent reps)
dv = cellfun(@(x) diff(x,[],3), out, 'UniformOutput', false);
[~,ns,nr_d] = size(dv{1});  % ns: #stims, nr_d: #rep-differences
nr  = nr_d + 1;              % #repetitions
nsj = numel(out);            % #subjects

% -------- correlation between drift vectors across stimuli
cm = zeros(nsj*(ns*(ns-1)/2), nr_d);
for r = 1:nr_d
    dat = cellfun(@(x) x(:,:,r), dv, 'UniformOutput', false);            % [units x stims]
    thismat = cell2mat(cellfun(@(x) 1-pdist(x',"correlation"), dat, ...
                               'UniformOutput', false));                 % pairwise across stims
    cm(:,r) = thismat(:);
end

% -------- normalized drift magnitudes per unit & stim
ds3 = zeros(sum(cellfun(@(x) size(x,1),out)), ns, nr_d);
i0 = 0;
for sj = 1:nsj
    n = size(out{sj},1);
    for r = 1:nr_d
        A = dv{sj}(:,:,r);                 % [units x stims]
        ds3(i0+(1:n),:,r) = A ./ sf(i0+(1:n)); % normalize by per-unit scale
    end
    i0 = i0 + n;
end

% -------- inter-stimulus distances per repetition (Euclidean)
isd = zeros(nsj*(ns*(ns-1)/2), nr);
for r = 1:nr
    dat = cellfun(@(x) x(:,:,r), out, 'UniformOutput', false);           % [units x stims]
    thismat = cell2mat(cellfun(@(x) pdist(x',"euclidean"), dat, ...
                                'UniformOutput', false));                % pairwise across stims
    isd(:,r) = thismat(:);
end

% -------- inter-repetition angles for each stim (per subject)
ira = [];
for sj = 1:nsj
    for s = 1:ns
        X = squeeze(out{sj}(:,s,:));       % [units x reps] trajectory over reps
        ira = [ira; triplet_angles(X)];    %#ok<AGROW>
    end
end

% -------- similarity of angle sequences across stimuli
ac = [];
for sj = 1:nsj
    T = [];
    for s = 1:ns
        X = squeeze(out{sj}(:,s,:));
        T = [T; triplet_angles(X)];        %#ok<AGROW>
    end
    ac = [ac; 1 - pdist(T,'cosine')];      %#ok<AGROW>
end

% -------- summary stats / bins for CDF plot
edges = -10:.05:10;
pval  = NaN;
if nr_d >= 2
    % Wilcoxon signed-rank: step 1–2 vs mean of later steps
    pval = signrank(cm(:,1), mean(cm(:,2:end),2));
end

% -------- pack results
R.out            = out;
R.scalingf       = sf;
R.driftv         = dv;
R.corrMat        = cm;
R.driftStrength  = reshape(ds3,[],nr_d);   % [units*stims x rep-steps]
R.interStimdist  = isd;
R.interRepAngle  = ira;
R.anglecorr      = ac;
R.edges          = edges;
R.signrank_p     = pval;

% -------- plots (plot functions accept precomputed inputs)
if o.plot
    F.angle_box      = plot_interRepAngle_box(ira,cfg);
    F.anglecorr_hist = plot_anglecorr_hist(ac,dft.stims_allowed,cfg);
    F.interstim_box  = plot_interStimdist_box(isd,dft.stims_allowed,cfg);
    F.drift_hists    = plot_driftStrength_hists(R.driftStrength,cfg);
    F.drift_cdf      = plot_driftStrength_cdf(R.driftStrength,edges,nr_d,dft.stims_allowed,cfg);
    F.corr_box       = plot_corrMat_box(cm,nr,dft.stims_allowed,cfg);
    R.figures        = F;
end
end

% =======================================================================
% Helpers
% =======================================================================

function a = triplet_angles(X)
%TRIPLET_ANGLES Angle between successive repetition steps per unit.
% INPUT:  X [units x reps]
% OUTPUT: a [1 x (reps-2)], angles (radians), NaN where step has zero norm
m = size(X,2);
v1 = X(:,2:m-1) - X(:,1:m-2);
v2 = X(:,3:m)   - X(:,2:m-1);
d  = sum(v1.*v2,1);
n1 = sqrt(sum(v1.^2,1));
n2 = sqrt(sum(v2.^2,1));
c  = d./(n1.*n2);
c  = max(-1,min(1,c));       % numeric guard
a  = acos(c);
a(n1==0 | n2==0) = NaN;
end

function f = plot_interRepAngle_box(ira,cfg)
% ira is [N x K] where K = nr-2
K = size(ira,2);
lbl = arrayfun(@(a,b) sprintf('%d-%d',a,b), 1:K, (1:K)+2, 'uni', false);

% Boxplot of inter-repetition angles (in units of pi)
f = figure; boxplot(ira./pi, 'Labels',lbl);
ylabel('inter-repetition angle (units of pi)');
axis square; box off; tcol(f,cfg);
end

function f = plot_anglecorr_hist(ac,ttl,cfg)
% Histogram of cosine similarity between angle sequences across stimuli
f = figure; histogram(ac,50,'FaceColor','k');
xlabel('cosine similarity btw drift trajectories across stimuli');
ylabel(''); axis square; box off; title(ttl); tcol(f,cfg);
end

function f = plot_interStimdist_box(isd,ttl,cfg)
% Boxplot of inter-stimulus Euclidean distances across repetitions
f = figure; boxplot(isd);
xlabel('repetitions'); ylabel('interstimulus eucl. distance');
box off; title(ttl); tcol(f,cfg);
end

function f = plot_driftStrength_hists(ds,cfg)
% Overlaid histograms of normalized drift magnitudes per rep-step
f = figure; hold on;
for r = 1:size(ds,2)
    histogram(ds(:,r));
end
hold off; tcol(f,cfg);
end

function f = plot_driftStrength_cdf(ds,edges,nr_d,ttl,cfg)
% Empirical CDFs of normalized drift magnitudes per rep-step
k = min(nr_d, size(ds,2));
f = figure; hold on; H = gobjects(k,1);
for r = 1:k
    c = histcounts(ds(:,r),edges);
    y = cumsum(c);
    if ~isempty(y) && y(end) > 0, y = y./y(end); end
    H(r) = plot(edges(1:end-1), y, 'LineWidth',1.5);
end
lbl = arrayfun(@(i) sprintf('%d-%d',i,i+1), 1:k, 'UniformOutput', false);
legend(H, lbl, 'Location','best');
hold off
axis square; xlabel('cell-wise norm. drift modulus'); ylabel('cumsum');
box off; title(ttl); tcol(f,cfg);
end

function f = plot_corrMat_box(cm,nr,ttl,cfg)
k   = min(size(cm,2), nr-1);
lbl = arrayfun(@(i) sprintf('%d-%d',i,i+1), 1:k, 'UniformOutput', false);

f = figure;
boxplot(cm(:,1:k), 'Labels', lbl);
xlabel('repetitions'); ylabel('drift vector correlation (btw. stimuli)');
box off; title(ttl); tcol(f,cfg);
end

function tcol(f,cfg)
% Apply dark/light theme
set(f,'color',cfg.bgcol);
ax = get(f,'CurrentAxes');
if ~isempty(ax)
    set(ax,'color',cfg.bgcol,'XColor',cfg.axcol,'YColor',cfg.axcol,'ZColor',cfg.axcol);
end
end
