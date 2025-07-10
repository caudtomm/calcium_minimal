function [params, y_fit, hf] = fitExpSaturation(t, y, pl)
% FITEXPSATURATION Fit y = a*(1 - exp(-b*t)) + c using fminsearch (no toolboxes needed)
% Inputs:
%   t - time vector
%   y - data vector
%   pl - boolean, optional plotting
% Outputs:
%   params - fitted parameters [a, b, c]
%   y_fit  - model values at input t

t = t(:);
y = y(:);

% Model function
model = @(p, t) p(1) * (1 - exp(-p(2) * t)) + p(3);

% Loss: sum of squared errors
loss = @(p) sum((model(p, t) - y).^2);

% Initial guess: [a, b, c]
a0 = max(y) - min(y);
b0 = 1 / (max(t)/2);
c0 = min(y);
p0 = [a0, b0, c0];

% Fit using fminsearch
params = fminsearch(loss, p0);

% Compute fitted curve
y_fit = model(params, t);

% Optional plot
hf = gobjects(1,1);
if ~pl; return; end
hf = figure;
plot(t, y, 'ko', 'DisplayName', 'Data'); hold on;
plot(t, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fit');
xlabel('Time');
ylabel('Value');
legend;
title('Exponential Saturation Fit (fminsearch)');
end
