function out = doModeDecomposition(data, varargin)
% doModeDecomposition performs dimensionality reduction on input data.
%
% Syntax:
%   out = doModeDecomposition(data, 'Name', Value, ...)
%
% Inputs:
%   data    - A double matrix of size [samples x variables] or [time x cells].
%
% Name-Value Pair Arguments:
%   'method'   - A string specifying the decomposition method. Options are:
%                'nmf' (default), 'pca', 'ica', 'rastermap', or 'dpca'.
%   'nfactors' - An integer specifying the number of factors to extract (default: 15).
%   'labs'     - Optional labels for the data (default: []).
%
% Outputs:
%   out - A structure with the following fields:
%         out.vals   - The reduced data representation.
%         out.coeffs - The coefficients or basis vectors.
%
% Example:
%   data = rand(100, 20); % Example data
%   out = doModeDecomposition(data, 'method', 'pca', 'nfactors', 5);
%
% Notes:
%   - The 'rastermap' and 'dpca' methods are not implemented and will throw an error.
%   - For 'ica', the FastICA algorithm is used.
arguments
    data double % [samples x variables] or [time x cells]
end
arguments (Repeating)
    varargin
end

% Set default values
method = 'nmf'; % 'nmf' or 'pca' or 'ica' or 'rastermap' or 'dpca'
nfactors = 15; % Number of factors to extract
knownLatents = [];

% Parse name-value pairs
if ~isempty(varargin)
    for k = 1:2:length(varargin)
        switch lower(varargin{k})
            case 'method'
                method = varargin{k+1};
            case 'nfactors'
                nfactors = varargin{k+1};
            case 'knownlatents'
                knownLatents = varargin{k+1};
        end
    end
end

% initialize output
out.vals = [];
out.coeffs = [];

switch lower(method)
    case 'nmf'
        [out.vals, out.coeffs] = nnmf(data, nfactors);
    case 'pca'
        [coeff, score, ~] = pca(data, 'NumComponents', nfactors);
        out.vals = score;
        out.coeffs = coeff;
    case 'ica'
        [out.vals, out.coeffs] = fastica(data', 'numOfIC', nfactors);
        out.vals = out.vals';
    case 'rastermap' % # TODO: Implement rastermap
        % Placeholder for rastermap implementation
        error('Rastermap method not implemented yet.');
    case 'dpca' % # TODO: Implement dPCA
        % Placeholder for dpca implementation (dPCA requires latents)
        error('dPCA method not implemented yet.');
    otherwise
        error('Unknown method: %s', method);

end