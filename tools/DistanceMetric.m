classdef DistanceMetric
    % Encapsulates a distance metric with its symmetry metadata.
    % Use static factory methods to construct instances.
    %
    % Examples:
    %   metric = DistanceMetric.fromPdist('correlation');  % symmetric
    %   metric = DistanceMetric.mahalanobis(labs);         % asymmetric
    %   metric = DistanceMetric.custom('mymetric', @myFn, false);

    properties
        name        (1,1) string
        isSymmetric (1,1) logical = true
        computeFn         function_handle
    end

    methods
        function D = compute(obj, X)
            % compute  Apply the metric to trial matrix X.
            % X : [n_trials x n_cells]
            % D : [n_trials x n_trials]
            D = obj.computeFn(X);
        end
    end

    methods (Static)
        function dm = fromPdist(method)
            % Wraps pdist + squareform. Symmetric by construction.
            dm = DistanceMetric;
            dm.name = string(method);
            dm.isSymmetric = true;
            dm.computeFn = @(X) squareform(pdist(X, char(method)));
        end

        function dm = mahalanobis(labs)
            % Class-conditional Mahalanobis. Asymmetric: D(i,j) is the
            % distance from trial i to the distribution of trial j's class,
            % using the covariance of that class.
            % labs : per-trial stimulus labels (cell array or string array)
            dm = DistanceMetric;
            dm.name = 'mahalanobis';
            dm.isSymmetric = false;
            dm.computeFn = @(X) classConditionalMahalanobis(X, labs);
        end

        function dm = custom(name, fn, isSymmetric)
            % Escape hatch for arbitrary metrics.
            % fn : function_handle (X [n_trials x n_cells]) -> [n_trials x n_trials]
            dm = DistanceMetric;
            dm.name = string(name);
            dm.isSymmetric = isSymmetric;
            dm.computeFn = fn;
        end
    end
end


function D = classConditionalMahalanobis(X, labs)
% D(i,j) = Mahalanobis distance from trial i to the class distribution
% of the stimulus that trial j belongs to (i.e., uses that class's cov).
% X    : [n_trials x n_cells]
% labs : {n_trials} cell of per-trial stimulus label strings

labs = cellstr(labs);
n    = size(X, 1);
D    = zeros(n);

unique_stims = unique(labs, 'stable');
mu       = struct;
SigmaInv = struct;

for s = unique_stims(:)'
    key = matlab.lang.makeValidName(s{1});
    idx = strcmp(labs, s{1});
    X_s = X(idx, :);
    mu.(key) = mean(X_s, 1);
    C = cov(X_s) + 1e-6 * eye(size(X_s, 2));  % Tikhonov regularization
    SigmaInv.(key) = inv(C);
end

for i = 1:n
    for j = 1:n
        if i == j; continue; end
        key    = matlab.lang.makeValidName(labs{j});
        diff   = X(i,:) - mu.(key);
        D(i,j) = sqrt(max(0, diff * SigmaInv.(key) * diff'));
    end
end
end
