function [samples, A, tau, C] = fitExponentialDecay(dataMatrix,samplerange)
    % Ensure that the input is a 2D matrix with [observations x time]
    [observations, timePts] = size(dataMatrix);

    % Define the exponential decay function
    expDecay = @(params, t) params(1) * exp(-t/params(2)) + params(3);

    % Initial guess (very naive)
    A0 = max(dataMatrix(:));
    tau0 = timePts / 2; 
    C0 = min(dataMatrix(:));
    initialParams = [A0, tau0, C0];

    % Time vector
    t = 1:timePts;

    % Store the estimated parameters for each observation
    A = zeros(observations, 1);
    tau = zeros(observations, 1);
    C = zeros(observations, 1);

    % Set options for lsqcurvefit
    options = optimoptions('lsqcurvefit', 'Display', 'off');

    % Loop through each observation and fit the exponential decay
    for obs = 1:observations
        data = dataMatrix(obs, :);

        % Interpolate missing values (NaNs) using linear interpolation
        validIndices = ~isnan(data);
        if sum(validIndices) > 1 % Ensure there are at least two non-NaN points for interpolation
            data = interp1(t(validIndices), data(validIndices), t, 'linear', 'extrap');
        else
            warning('Insufficient data for observation %d. Skipping.', obs);
            continue;
        end

        % Estimate parameters using lsqcurvefit
        fitparams = lsqcurvefit(expDecay, initialParams, t, data, [], [], options);

        % Store the estimated parameters
        A(obs) = fitparams(1);
        tau(obs) = fitparams(2);
        C(obs) = fitparams(3);
    end
    
    if ~exist('samplerange','var'); samplerange = t; end
    samples = generateExponentialDecaySample(samplerange, A, tau, C);
end

function y = generateExponentialDecaySample(t, A, tau, C)
    % Exponential decay sample
    y = A .* exp(-t./tau) + C;
end


