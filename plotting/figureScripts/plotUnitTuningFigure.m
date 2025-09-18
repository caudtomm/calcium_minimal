function [hf, out, groups, odor_sets] = plotUnitTuningFigure(v, method, do_animate)  
% plotUnitTuningFigure - Plot single unit tuning distributions and low-dimensional embeddings.
%
% Usage:
%   [hf, out, groups, odor_sets] = plotUnitTuningFigure(v, method, do_animate)
%
% Inputs:
%   v          - ExperimentViewer object containing data and plotting methods.
%   method     - String specifying embedding method ('pca' or 'isomap', default: 'pca').
%   do_animate - Logical flag to animate embedding transitions (default: false).
%
% Outputs:
%   hf         - Array of figure handles.
%   out        - Cell array of output data from plots.
%   groups     - Cell array of subject group names.
%   odor_sets  - Cell array of odor set names.
    arguments
        v ExperimentViewer
        method = 'pca'
        do_animate logical = false
    end
    
    % Knobs
    groups = {...'naïve', ...
                'trained'}%, ...
                % 'uncoupled'};
    odor_sets = {'all stimuli', ...
                'all familiar', ...
                'all novel', ...
                };

    % useful metrics
    ngroups = numel(groups);
    nodor_sets = numel(odor_sets);
    nplots = ngroups*nodor_sets;

    % Initialize output
    hf = gobjects(10,1);
    i_hf = 0;
    out = cell(nplots,1);

    %% Figure 1: comparison of single unit tuning distributions across repetitions (histogram outlines)
    % i_hf = i_hf+1;
    % plotType = 'boxplot_repetitions';
    % [data,~,n_subplots] = plotGrid;
    % set(gcf,'Position',[1 1 300 1000])

    %% Figure 2: low dimensional embedding of single unit tuning across repetitions,
    % based on common axes ( one figure for each subplot from Figure 1)
    
    n = 0;
    for i_g = 1:ngroups
        thisgroup = groups{i_g};
        v.dataFilter.subjectGroup = thisgroup;
        
        for i_os = 1:nodor_sets
            n = n+1;
            thisodorset = odor_sets{i_os};
            
            i_hf = i_hf+1;
            hf(i_hf) = figure;
            cfg = v.plotConfig; % # TODO: use cfg to set up the figure

            [thisout.embedding,thisout.cols,thisout.gammas,thisout.d] = plotLDEmbedding();
            out{n} = thisout;
            if isempty(thisout.embedding); continue; end

            title([thisgroup,', ',thisodorset])
        end
    end

    function [Y,allc,allgamma,D] = plotLDEmbedding()
        % # TODO : add no_plot flag to only retrieve data (plotType = 'none')
        v.dataFilter.stims_allowed = thisodorset;
        s = v.plotUnitActivityMetricHead('plotType', 'none', ...
                                'method', 'selectivity of tuning repetitions', ...
                                'n_equals', 'cells', ...
                                'do_normalize', false); % cell, each [units x reps]
        tc = v.plotUnitActivityMetricHead('plotType', 'none', ...
                                'method', 'tuning curves', ...
                                'n_equals', 'cells', ...
                                'do_normalize', false); % cell, each [units x stims x reps]
    
        % metrics
        nSubjects = numel(tc);
    
        % pool subjects
        thistc = []; thiss = [];
        for i = 1:nSubjects
            thistc = [thistc ; tc{i}]; % [nUnits,nStims,nReps]
            thiss = [thiss; s{i}]; % [nUnits,nReps]
        end
        [nUnits,nStims,nReps] = size(thistc);

        if isempty(thistc); return; end
        
        D = [];
        Y = [];
        allgamma = [];
        for i_ref = 1:nReps 
            rep_oi = i_ref; % rep of interest
            refdata = thistc(:,:,rep_oi); % [units x stims]
    
            % get low-dimensional embedding axes based on reference repetition
            % then transform all repetitions to this reference
            y = nan(nUnits,2,nReps);
            switch method
                case 'pca'
                    refdata = nanzscore(refdata,[],2); % zscore
                    thistc = nanzscore(thistc,[],2); % zscore
                    refdata(isnan(refdata))=0;
                    thistc(isnan(thistc))=0;

                    mdl = pca(refdata); % mdl is PCA coefficients [nStims x nStims]
                    mdl = mdl(:,1:2); % cut only first 2 PCs
    
                    % transform all repetitions
                    for i_rep = 1:nReps
                        y(:,:,i_rep) = thistc(:,:,i_rep)*mdl; % project new samples (nUnits × nPCs)
                    end

                case 'isomap'
                    disp('isomap fit') 

                    refdata = nanzscore(refdata,[],2); % zscore
                    thistc = nanzscore(thistc,[],2); % zscore
                    refdata(isnan(refdata))=0;
                    thistc(isnan(thistc))=0;
                    
                    mdl = isomap_fit(refdata,'k',10,2);  % (k=10 neighbors, 2D embedding)
                    
                    % transform all repetitions
                    for i_rep = 1:nReps
                        % verbose
                        disp(['isomap transform #',num2str(i_rep)])
                        y(:,:,i_rep) = isomap_transform(thistc(:,:,i_rep),mdl); % project new samples (nUnits × 2)
                    end
                    
                otherwise
                    error('Unknown method for low-dimensional embedding');
            end

            Y(:,:,:,i_ref) = y;
            
    
            % plotting
    
            % call figure handle
            figure(hf(i_hf));
            
            % define graphics:
            % cell color: preferred stimulus
            % gamma: selectivity index
            % based on one given rep
            c = parula(nStims);
            [~,pref_stim] = max(thistc,[],2);
            pref_stim = squeeze(pref_stim); % [N x reps]
            thisc = c(pref_stim(:,rep_oi),:);
            thisgamma = thiss(:,rep_oi)/2 + .5; % lim 0.5->1

            allc(:,:,rep_oi) = thisc;
            allgamma(:,rep_oi) = thisgamma;
    
            % plot a row of sublots
            for i_rep = 1:nReps
                subplot(nReps,nReps,nReps*(i_ref-1)+i_rep)
        
                b = scatter(y(:,1,i_rep), ...
                    y(:,2,i_rep), ...
                    20, thisc, 'filled');
                b.AlphaData = thisgamma;
                b.MarkerFaceAlpha = 'flat';
        
                axis square
                set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol);
                set(gcf, 'color', cfg.bgcol); 
                set(gcf,'Position',[1 1 1000 2000])
                hold off
            end

            % quantify distance traveled in embedding space between
            % repetitions for each unit
            % y : [units x 2 x reps]
            y_diff = diff(y,1,3); % dx and dy
            d = squeeze(sqrt(y_diff(:,1,:).^2 + y_diff(:,2,:).^2)); % distances [units x reps-1]
            D(:,:,i_ref) = d;
    
            % (optional) plot animated scatter
            if ~do_animate; continue; end
    
            figure
            FileOut = ['refrep_',num2str(i_ref),'.gif'];
            animate_reps(permute(y,[2,1,3]),thisc,thisgamma,60,2,gca,'g',.05,false,[1 2])
    
    
        end

        % accessory figure
        figure
        for i = 1:nReps
            subplot(nReps,1,i); violin(D(:,:,i));
            legend off
            box off
            set(gca, 'color', cfg.bgcol, 'XColor',cfg.axcol, 'YColor',cfg.axcol, 'ZColor',cfg.axcol);
            set(gcf, 'color', cfg.bgcol); 
        end
    end


end


function h=animate_reps2(Y,c,a,fps,dur,ax,track,pad,eq,dims)
% Y : [2 x units x reps]
[~,nUnits,nReps] = size(Y);
if nargin<2||isempty(c),   c= zeros(nUnits,3); end % all black by default
if nargin<3||isempty(a),   a=1; end
if nargin<4||isempty(fps), fps=60; end
if nargin<5||isempty(dur), dur=100; end
if nargin<6||isempty(ax),  ax=gca; end
if nargin<7||isempty(track), track='g'; end
if nargin<8||isempty(pad), pad=.05; end
if nargin<9||isempty(eq),  eq=false; end
if nargin<10||isempty(dims), dims=[1 2]; end

v=size(Y,1); n=size(Y,2); r=size(Y,3); assert(all(dims<=v)&&numel(dims)==2,'dims bad')
X=cell(r,1); Yy=cell(r,1);
for k=1:r
    X{k}=squeeze(Y(dims(1),:,k));
    Yy{k}=squeeze(Y(dims(2),:,k));
end
mn=@(u)min(u,[],'omitnan'); mx=@(u)max(u,[],'omitnan');
gx=[mn(cellfun(@(u) mn(u),X)) mx(cellfun(@(u) mx(u),X))];
gy=[mn(cellfun(@(u) mn(u),Yy)) mx(cellfun(@(u) mx(u),Yy))];
rx=gx(2)-gx(1); ry=gy(2)-gy(1); if rx==0, rx=1; end; if ry==0, ry=1; end
gx=gx+[-1 1]*rx*pad; gy=gy+[-1 1]*ry*pad;

axes(ax)
h=scatter(X{1},Yy{1},[],c,'filled');
if numel(a)>1, h.AlphaData=a; h.MarkerFaceAlpha='flat'; else, h.MarkerFaceAlpha=a; end
if eq, axis(ax,'equal'); end
xlim(gx); ylim(gy); box(ax,'on');

F=max(1,round(fps*dur));
for k=1:r-1
    x0=X{k}; y0=Yy{k}; x1=X{k+1}; y1=Yy{k+1};
    if track=='p'
        xl=[mn([x0(:);x1(:)]) mx([x0(:);x1(:)])];
        yl=[mn([y0(:);y1(:)]) mx([y0(:);y1(:)])];
        rx=xl(2)-xl(1); ry=yl(2)-yl(1); if rx==0, rx=1; end; if ry==0, ry=1; end
        xlim(xl+[-1 1]*rx*pad); ylim(yl+[-1 1]*ry*pad);
    else
        xlim(gx); ylim(gy);
    end
    for f=1:F
        t=f/F; s=.5-.5*cos(pi*t);
        set(h,'XData',x0+(x1-x0)*s,'YData',y0+(y1-y0)*s);
        drawnow limitrate
    end
end
end

function h=animate_reps(Y,c,a,fps,dur,ax,track,pad,eq,dims,gif)
% Y : [2 x units x reps]
[~,nUnits,nReps] = size(Y);
if nargin<2||isempty(c),   c= zeros(nUnits,3); end
if nargin<3||isempty(a),   a=1; end
if nargin<4||isempty(fps), fps=60; end
if nargin<5||isempty(dur), dur=100; end
if nargin<6||isempty(ax),  ax=gca; end
if nargin<7||isempty(track), track='g'; end
if nargin<8||isempty(pad), pad=.05; end
if nargin<9||isempty(eq),  eq=false; end
if nargin<10||isempty(dims), dims=[1 2]; end
if nargin<11, gif=''; end

v=size(Y,1); r=size(Y,3); 
assert(all(dims<=v)&&numel(dims)==2,'dims bad')

X=cell(r,1); Yy=cell(r,1);
for k=1:r
    X{k}=squeeze(Y(dims(1),:,k));
    Yy{k}=squeeze(Y(dims(2),:,k));
end

mn=@(u)min(u,[],'omitnan'); mx=@(u)max(u,[],'omitnan');
gx=[mn(cellfun(@(u) mn(u),X)) mx(cellfun(@(u) mx(u),X))];
gy=[mn(cellfun(@(u) mn(u),Yy)) mx(cellfun(@(u) mx(u),Yy))];
rx=gx(2)-gx(1); ry=gy(2)-gy(1); if rx==0, rx=1; end; if ry==0, ry=1; end
gx=gx+[-1 1]*rx*pad; gy=gy+[-1 1]*ry*pad;

axes(ax)
h=scatter(X{1},Yy{1},[],c,'filled');
if numel(a)>1, h.AlphaData=a; h.MarkerFaceAlpha='flat'; else, h.MarkerFaceAlpha=a; end
if eq, axis(ax,'equal'); end
xlim(gx); ylim(gy); box(ax,'on');

F=max(1,round(fps*dur));

% initialize gif if requested
if ~isempty(gif)
    fr=getframe(ax); im=frame2im(fr); [A,map]=rgb2ind(im,256);
    imwrite(A,map,gif,'gif','LoopCount',inf,'DelayTime',1/fps);
end

for k=1:r-1
    x0=X{k}; y0=Yy{k}; x1=X{k+1}; y1=Yy{k+1};
    if track=='p'
        xl=[mn([x0(:);x1(:)]) mx([x0(:);x1(:)])];
        yl=[mn([y0(:);y1(:)]) mx([y0(:);y1(:)])];
        rx=xl(2)-xl(1); ry=yl(2)-yl(1); if rx==0, rx=1; end; if ry==0, ry=1; end
        xlim(xl+[-1 1]*rx*pad); ylim(yl+[-1 1]*ry*pad);
    else
        xlim(gx); ylim(gy);
    end
    for f=1:F
        t=f/F; s=.5-.5*cos(pi*t);
        set(h,'XData',x0+(x1-x0)*s,'YData',y0+(y1-y0)*s);
        drawnow limitrate
        if ~isempty(gif)
            fr=getframe(ax); im=frame2im(fr); [A,map]=rgb2ind(im,256);
            imwrite(A,map,gif,'gif','WriteMode','append','DelayTime',1/fps);
        end
    end
end
end




function m = isomap_fit(X, nfun, nsz, d)
% Fit Isomap on a reference dataset and cache everything needed
% to project new points later on the *same* axes (out-of-sample).
%
% X    : N×P matrix (rows = samples, cols = variables)
% nfun : 'k' or 'epsilon' (neighborhood rule)
% nsz  : k (if 'k') or epsilon radius (if 'epsilon')
% d    : target embedding dimension (default 2)
%
% Returns a struct m with fields:
%   m.mu, m.sg   : z-score params from X
%   m.X          : reference samples used (possibly reduced to largest component), z-scored
%   m.A          : symmetric weighted adjacency of largest connected component
%   m.V, m.L     : top eigenvectors and eigenvalues of centered G^2 (for MDS)
%   m.r, m.g     : row means and grand mean of G^2 (for out-of-sample centering)
%   m.d          : actual dimensionality kept (>0 eigenvalues)
%   m.idx        : indices of X kept after selecting the largest component
%   m.Y          : reference embedding coordinates (N_kept × m.d)
%   m.nfun, m.nsz: neighborhood rule & size (for consistent projection)

if nargin<4, d=2; end

% 1) Standardize variables (z-score) so distances are comparable.
mu = mean(X,1);
sg = std(X,0,1); sg(sg==0)=1;
Z = (X - mu) ./ sg;
Z(isnan(Z)) = 0;                 % mean-impute missing features in z-space

% 2) Pairwise ambient distances among reference points (Euclidean).
D = pdist2(Z,Z);
N = size(D,1);

% 3) Build a local neighborhood graph (weighted by ambient distances).
if strcmp(nfun,'k')
    % Directed kNN (exclude self), then symmetrize later.
    [~,ix] = sort(D,2);
    k = nsz;
    I = repmat((1:N)',1,k);      % row indices i
    J = ix(:,2:k+1);             % k nearest neighbors per row
    V = D(sub2ind([N N],I(:),J(:)));
    A = sparse(I(:),J(:),V,N,N); % weighted adjacency (directed for now)
elseif strcmp(nfun,'epsilon')
    % Epsilon graph: connect pairs within radius nsz (exclude diagonal).
    M = (D<=nsz) & ~eye(N);
    [I,J] = find(M);
    V = D(M);
    A = sparse(I,J,V,N,N);
else
    error('nfun must be ''k'' or ''epsilon''');
end

% 4) Make the graph undirected by taking the max weight in each direction.
A = max(A,A');

% 5) Keep only the largest connected component (Isomap requires connectivity).
assert(~any(isnan(nonzeros(A))), 'A contains NaN weights');
assert(~any(isinf(nonzeros(A))), 'A contains Inf weights');
G = graph(A);
b = conncomp(G);
c = mode(b);             % label of the largest component
idx = find(b==c);        % indices that belong to it
A = A(idx,idx);
Z = Z(idx,:);

% 6) Geodesic distances on the graph (shortest paths with edge weights = ambient distances).
A(~isfinite(A)) = 0;             % no NaN/Inf weights
A = max(A,A'); % ensure symmetry, no self-loops, and use 'upper'
A = A - diag(diag(A));
G = graph(triu(A,1),'upper');
DG = distances(G);       % N_kept × N_kept geodesic distance matrix

% 7) Classical MDS on geodesic distances: double-center D^2 to get inner-product matrix B,
%    then take its top eigenpairs (this is the Isomap embedding).
G2 = DG.^2;
r = mean(G2,2);          % row means
g = mean(G2(:));         % grand mean
B = -0.5*(G2 - r*ones(1,size(G2,1)) - ones(size(G2,1),1)*r' + g);

% Symmetrize numerically and take top eigenpairs
[V,L] = eigs((B+B')/2, d, 'largestreal');
lv = diag(L);

% 8) Keep only strictly positive eigenvalues (valid Euclidean axes).
p  = lv>0;
V  = V(:,p);
lv = lv(p);

% 9) Reference embedding coordinates (rows = samples, cols = dims).
Y = V*diag(sqrt(lv));

% 10) Pack the "model" used for out-of-sample projections.
m.X   = Z;
m.A   = A;
m.V   = V;
m.L   = lv;
m.r   = r;
m.g   = g;
m.mu  = mu;
m.sg  = sg;
m.nfun = nfun;
m.nsz  = nsz;
m.d    = numel(lv);
m.idx  = idx;
m.Y    = Y;
end

function Y = isomap_transform(Xnew, m)
% Project new samples onto the fixed Isomap axes learned by isomap_fit.
%
% Xnew : M×P new data (same variables & scaling as the fit set)
% m    : struct returned by isomap_fit
%
% Returns:
%   Y : M×m.d out-of-sample coordinates on the same axes as m.Y
%       (rows that could not be connected will be NaN)

% 1) Apply the *same* z-score used for the reference set.
Z = (Xnew - m.mu)./m.sg;
Z(isnan(Z)) = 0;                  % mean-impute new samples

N = size(m.X,1);   % number of reference nodes (largest component)
M = size(Z,1);     % number of new samples
Y = nan(M, m.d);   % preallocate (NaN for failures to connect)

for t=1:M
    % 2) Distances from the new point to all reference nodes (ambient space).
    d0 = pdist2(Z(t,:), m.X);

    % 3) Attach the new point to the reference graph using the same rule.
    if strcmp(m.nfun,'k')
        [~,ix] = sort(d0);
        nb = ix(1:m.nsz);           % k nearest reference neighbors
    else
        nb = find(d0<=m.nsz);       % all neighbors within epsilon
    end
    if isempty(nb)
        % No neighbors found → cannot connect → leave as NaN.
        continue;
    end

    % 4) Form a temporary augmented graph with the new node as (N+1).
    S = sparse(N+1,N+1);
    S(1:N,1:N) = m.A;               % existing reference graph
    S(N+1,nb)  = d0(nb);            % connect new node to neighbors (weights = ambient distances)
    S(nb,N+1)  = d0(nb);            % undirected

    % 5) Single-source shortest paths from the new node to all references (geodesics).
    S = max(S,S');% ensure symmetry + no self-loops, then treat as undirected
    S = S - diag(diag(S));
    assert(~any(isnan(nonzeros(S))), 'S contains NaN weights');
    assert(~any(isinf(nonzeros(S))), 'S contains Inf weights');
    G = graph(triu(S,1),'upper');
    del = distances(G, N+1, 1:N);   % 1×N vector
    d2 = del'.^2;                   % N×1 squared geodesic distances

    % 6) Nyström-style out-of-sample formula for classical MDS:
    %    c = centered cross-term between the new point and reference points.
    c = -0.5*(d2 - m.r - mean(d2) + m.g);   % N×1

    % 7) Project onto the learned axes: y = Λ^{-1/2} V^T c
    y = (m.V' * c) ./ sqrt(m.L);    % m.d × 1

    % 8) Store row-wise.
    Y(t,1:numel(y)) = y';
end
end
