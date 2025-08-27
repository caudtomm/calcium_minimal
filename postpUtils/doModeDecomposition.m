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
%   - The 'ica', 'rastermap' and 'dpca' methods are not implemented and will throw an error.
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
SuS = []; % Suppression Score
SoT = []; % Selectivity of Tuning
SoTzero = []; % Selectivity of Tuning 'zero'

% Parse name-value pairs
if ~isempty(varargin)
    for k = 1:2:length(varargin)
        switch lower(varargin{k})
            case 'method'
                method = varargin{k+1};
            case 'nfactors'
                nfactors = varargin{k+1};
            
            % additional variables needed for some methods
            case 'knownlatents'
                knownLatents = varargin{k+1};
            case 'suppression score'
                SuS = varargin{k+1};
            case 'selectivity of tuning'
                SoT = varargin{k+1};
            case 'sot_zero'
                SoTzero = varargin{k+1};
        end
    end
end

% initialize output
out.vals = [];
out.coeffs = [];

switch lower(method)
    case 'nmf'
        % NMF employs the concept of low-rank approximation with respect to the feature space, 
        % and it would require that the matrix be void of empty elements/values.
        % Therefore, we interpolate missing values using a spline alg.
        data = fillmissing(data,'constant',0);

        options = statset; options.Display = 'final';
        [out.vals, out.coeffs, out.sqresiduals] = ...
            nnmf(nanzscore(data), nfactors,options=options);
        out.coeffs = out.coeffs';
    case 'pca'
        data = fillmissing(data,'constant',0);
        [coeff, score, ~] = pca(nanzscore(data), 'NumComponents', nfactors);
        out.vals = score;
        out.coeffs = coeff;
    case 'ica' % # TODO: Implement ICA
        % Placeholder for ICA implementation
        error('ICA method not implemented yet.');
    case 'rastermap' % # TODO: Implement rastermap
        % Placeholder for rastermap implementation
        error('Rastermap method not implemented yet.');
    case 'dpca' % # TODO: Implement dPCA
        % Placeholder for dpca implementation (dPCA requires latents)
        error('dPCA method not implemented yet.');
    case 'general suppression score vs tuning selectivity'
        % This method requires specific data structure and is not a standard decomposition method.
        if isempty(SuS) || isempty(SoT)
            error('For "general suppression score vs tuning selectivity", both SuS and SoT must be provided.');
        end
        out.suppressionScores = SuS;
        out.selectivityOfTuning = SoT;
        
        % four non-overlapping cell subset ('modes') are identified, based on the four quadrants of the scatter plot
        % of the suppression score vs tuning selectivity, centered around the origin.
        % Additionally, a score is associated with each cell, indicating how 'extreme' the cell is in its quadrant.

        SuSzero = 0; % definition of separating value for Suppression Score
        if isempty(SoTzero) % default value
            SoTzero = median(SoT,'omitmissing'); % definition of separating value for Selectivity of Tuning
        end
        out.SuSzero = SuSzero;
        out.SoTzero = SoTzero;

        out.upperLeft = SuS > SuSzero & SoT < SoTzero;
        out.upperLeftScore = SuS ./ SoT;

        out.upperRight = SuS > SuSzero & SoT > SoTzero;
        out.upperRightScore = SuS .* SoT;

        out.lowerLeft = SuS < SuSzero & SoT < SoTzero;
        out.lowerLeftScore = 1 ./ (SuS .* SoT);

        out.lowerRight = SuS < SuSzero & SoT > SoTzero;
        out.lowerRightScore = SoT ./ SuS;

        out.unassigned = ~(out.upperLeft | out.upperRight | out.lowerLeft | out.lowerRight);

        % coeff treats the subsets as modes
        out.coeffs = [out.upperLeft, out.upperRight, out.lowerLeft, out.lowerRight];
        
        out.coeffs2 = [out.upperLeftScore, out.upperRightScore, out.lowerLeftScore, out.lowerRightScore];
        out.coeffs2 = out.coeffs .* out.coeffs2;

        % vals is the average data value for each mode
        out.vals = fillmissing(data, 'constant',0) * out.coeffs;
        
    otherwise
        error('Unknown method: %s', method);

end