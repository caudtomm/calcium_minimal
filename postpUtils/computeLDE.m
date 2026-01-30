function out = computeLDE(data, labs, varargin)
arguments
    data cell % cell array of [variables x samples] or [cells x trials]; or [time x variables x samples]
    labs cell % cell array of cell arrays {samples x 1} or {trials x 1} of strings
end
arguments (Repeating)
    varargin
end

% Set default values
params.method = 'pca'; % 'pca' or 'umap'
params.poolData = false; % pool data across subjects
params.zscoreData = true; % z-score data before LDE computation
params.zscoreDim = 2; % dimension along which to z-score
params.nans2zeros = false; % replace NaNs with zeros before z-scoring
% UMAP-specific parameters
params.metric = 'cosine'; % metric for UMAP
params.min_dist = 0.8; % min_dist for UMAP
params.n_components = 2; % number of components for UMAP
params.n_neighbors = 199; % n_neighbors for UMAP
params.init = 'spectral'; % initialization for UMAP

% Parse name-value pairs
if ~isempty(varargin)
    for k = 1:2:length(varargin)
        switch lower(varargin{k})
            case 'method'
                params.method = varargin{k+1};
            case 'pooldata'
                params.poolData = varargin{k+1};
            case 'zscoredata'
                params.zscoreData = varargin{k+1};
            case 'zscoredim'
                params.zscoreDim = varargin{k+1};
            case 'nans2zeros'
                params.nans2zeros = varargin{k+1};
            case 'metric'
                params.metric = varargin{k+1};
            case 'min_dist'
                params.min_dist = varargin{k+1};
            case 'n_components'
                params.n_components = varargin{k+1};
            case 'n_neighbors'
                params.n_neighbors = varargin{k+1};
            case 'init'
                params.init = varargin{k+1};
            otherwise
                error('Unknown parameter name: %s', varargin{k});
        end
    end
end


% Preprocess data
data = cellfun(@(x) preprocessData(x, params), data, 'UniformOutput', false);

if params.poolData
    % Pool data across subjects    
    mytraces = cellfun(@(x) x.traces,data,'UniformOutput',false);
    mytraces = mytraces(:)'; % ensure horizontal cell array
    mytraces = cell2mat(mytraces);

    data = data(1); % init to retain ntrials param
    data{1}.traces = mytraces;

    labs = labs(1); % use labs from first subject only
end

% Compute LDE
switch params.method
    case 'pca'
        embedding = cellfun(@(x) computePCA(x.traces), data, 'UniformOutput', false);
    case 'umap'
        embedding = cellfun(@(x) computeUMAP(x.traces, params), data, 'UniformOutput', false);
    case 'returnonly'
        embedding = data; % skip computation, return preprocessed data only
        for i = 1:numel(embedding)
            embedding{i}.reduction = embedding{i}.traces;
        end
    otherwise
        error('Requested LDE method is unknown.')
end

% Reshape embedding to [length x dimensions x trials]
for i = 1:numel(embedding)
    ntrials = data{i}.ntrials;
    embedding{i}.reduction = ActivityTraces.format(embedding{i}.reduction, ntrials);
end

% return
out.embedding = embedding;
out.inputData = data;
out.labs = labs;
out.params = params;
out.name = sprintf('LDE_%s', params.method);

end

function out = preprocessData(data, params)
    % Replace NaNs with zeros if specified
    if params.nans2zeros
        data(isnan(data)) = 0;
    end

    % Concatenate trials and output ntrials parameter
    ntrials = size(data, 3);
    if ntrials>1
        data = ActivityTraces.format(data);
    end

    % Z-score data if specified
    if params.zscoreData
        data = nanzscore(data, 0, params.zscoreDim);
    end

    out.traces = data;
    out.ntrials = ntrials;
end

function embedding = computePCA(data)
    % Compute PCA on the preprocessed data
    [coeff, score, latent, tsquared, explained, mu] = pca(data);

    embedding.coeff = coeff;
    embedding.reduction = score;
    embedding.latent = latent;
    embedding.tsquared = tsquared;
    embedding.explained = explained;
    embedding.mu = mu;
end

function embedding = computeUMAP(data, params)
    % remove lines with NaNs
    idx_nan = any(isnan(data),2);
    data(idx_nan,:) = [];

    % Compute UMAP on the preprocessed data
    % Requires UMAP package to be installed
    [reduction, umap, clusterIdentifiers] = run_umap(data, ...
        'metric',params.metric, ...
        'min_dist',params.min_dist, ... % .8 for cosine, .25 for euclidean
        'n_components',params.n_components, ...
        'n_neighbors',params.n_neighbors, ...
        'init',params.init);
    close % close UMAP figure

    reduction = insert_rows(reduction,find(idx_nan));

    embedding.reduction = reduction;
    embedding.umap = umap;
    embedding.clusterIdentifiers = clusterIdentifiers;
end

function mat_out = insert_rows(mat_in,original_idx,original_rows)
    % original_idx of type double (not boolean indices)

    ndim = size(mat_in,2);

    if ~exist('original_idx',"var"); mat_out = mat_in; return;  end
    if ~exist('original_rows',"var"); original_rows = nan(length(original_idx),ndim); end
    if size(original_rows,2)~=ndim; error('number of columns is inconsistent'); end

    mat_out = nan( size(mat_in,1)+size(original_rows,1) , ndim );

    for i = 1:size(mat_out,1)
        thisrow=[];
        position_in_idxmat = find(original_idx==i);
        if isempty(position_in_idxmat)
            thisrow = mat_in(1,:);
            mat_in(1,:) = [];
        else
            thisrow = original_rows(position_in_idxmat,:);
        end
        mat_out(i,:) = thisrow;
    end
end