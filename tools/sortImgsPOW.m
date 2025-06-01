function [results, hf] = sortImgsPOW(images, pl)
% output: low->high freq power biasobj
arguments
    images cell % containing double images
    pl logical = false % optional plotting
end

% knobs
th = .1; % normalized power thresh -> for bias score

%init vars
N = numel(images);

%% resize images to be all the same (upscale to largest)


%% 
bias_scores = zeros(N,1);
[Ps, rpps] = deal(cell(N,1));

for i = 1:N
    img = images{i};

    % Handle NaNs
    if any(isnan(img(:)))
        m = mean(img,'all','omitmissing');
        img(isnan(img)) = m;
    end
    
    % FFT
    F = fft2(img);
    F = fftshift(F);
    P = abs(F).^2;

    % Radial frequency map
    [h, w] = size(P);
    [X, Y] = meshgrid(1:w, 1:h);
    cx = (w+1)/2; cy = (h+1)/2;
    R = sqrt((X - cx).^2 + (Y - cy).^2);

    % Radial power profile
    Rint = round(R);
    rpp = accumarray(Rint(:)+1, P(:)) ./ accumarray(Rint(:)+1, 1);
    rpp = rpp(2:end); % first value is usually NaN or otherwise unreliable

    % Normalize radius to [0,1]
    R = R ./ max(R(:));

    % Define low/high threshold (e.g. 0.5)
    low = P(R <= th);
    high = P(R > th);
    
    % Compute bias (low power ratio)
    bias_scores(i) = sum(low(:)) / (sum(low(:)) + sum(high(:)));

    % export spectra (for plotting)
    Ps{i} = P;
    rpps{i} = rpp;
end

% check for Ps and rpps all of the same size across images
maxlen = max(cellfun(@numel, rpps));
minlen = min(cellfun(@numel, rpps));
if diff([minlen maxlen])~=0
    Ps
    rpps
    error(['something went wrong: power spectra aren''t ' ...
        'of equal size across images'])
end

% store to double arrays
[h,w] = size(Ps{1});
Ps = cell2mat(Ps);
Ps = permute(reshape(Ps',w,h,N),[2,1,3]); % [h,w,N]

[h,w] = size(rpps{1});
rpps = cell2mat(rpps);
rpps = squeeze(permute(reshape(rpps',w,h,N),[2,1,3])); % [radii,N]

%% extract additional metrics

rpps_noneg = rpps; rpps_noneg(rpps_noneg<0)=0;
ref = mean(rpps_noneg,2,'omitmissing'); ref = ref / sum(ref,"omitmissing");

% rpp = rpps(:,i);
n = height(rpps);

% calc centroids
centroids = zeros(N,1);
for i = 1:N
    thiscentroid = regionprops(true(height(rpps),1), rpps(:,i), ...
        'WeightedCentroid').WeightedCentroid;    
    centroids(i) = thiscentroid(2);
end

% calc spectral slope (in loglog space)
slopes = zeros(N,1);
r = (1:n)';
for i = 1:N
    rpp = rpps(:,i);
    valid = rpp>0;
    p = polyfit(log(r(valid)), log(rpp(valid)),1);
    slopes(i) = p(1);
end
slopes = slopes.^2;

% calc KL divergence from reference spectrum
% normalize
rppnorm = rpps_noneg ./ sum(rpps_noneg,'omitmissing');
KL_divergence = sum(rppnorm .* log(rppnorm ./ ref), 'omitmissing')';

%% return results

spectral_bias = table(bias_scores, ...
                      centroids, ...
                      slopes, ...
                      KL_divergence);

results.spectral_bias = spectral_bias;
results.pow_spectra = Ps;
results.mesh_normalized = R;
results.rad_pow_profiles = rpps;
results.norm_profiles = rppnorm;
results.KLref = ref;

%% PLOTTING
hf = gobjects(1,4);
if ~pl; return; end

hf(1) = figure; % power spectra
ncols = floor(N/2)+1;
nrows = floor(N/(ncols-1));

t = tiledlayout(nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');

[~,idx] = sort(bias_scores, 'descend');
for i = 1:N
    nexttile
    imagesc(log(1+Ps(:,:,idx(i)))); axis image off
    colormap('jet'); %colorbar
    title(['Log PowSp (#',num2str(idx(i)),')'])
end

hf(2) = figure; % distribution of scores
x = 1+ .5*(rand(numel(bias_scores),1)-.5);
scatter(x,bias_scores,"filled");
ylim([0 1]); xlim([.5 1.5])
xticks(1); xticklabels({'Spectral Bias Scores'});
for i = 1:numel(bias_scores) % Add index labels
    text(x(i) + 0.03, bias_scores(i), num2str(i), 'FontSize', 8,'Color','b');
end
set(hf(2),"Position",[1 1 200 300])

hf(3) = figure; % radial power profiles (normalized) and spectral centroid
subplot(141)
n = height(rpps);
x = 0:1/n:1-1/n;
semilogy(x,rpps)
xlabel('Radius')
ylabel('Average Power')
title('Radial Power Profiles')

subplot(142)
semilogy(x,rppnorm); hold on
b = semilogy(x,ref,'r--','LineWidth',2.5);
legend(b,'reference');
xlabel('Radius')
ylabel('Average Power')
title('Norm. Profiles')

subplot(143) % plot again with centroids
imagesc(log(rpps(:,idx)')); colorbar % spectral profiles
hold on
y = 1:N;
b = scatter(centroids(idx), y, 50, 'r' , 'filled'); % centroids
ax = gca; try; ax.XScale = 'log'; catch; end
axis tight
legend(b,'centroids')
xticks([])
xlabel('Radius')
yticks(1:numel(idx)); yticklabels(idx)
ylabel('Image #')
title('Log Radial Power Profile')

subplot(144) % spectral slope distribution
x = 1+ .5*(rand(numel(slopes),1)-.5);
scatter(x,slopes,"filled");
xlim([.5 1.5])
xticks([]); title('Squared Spectral Slopes')
for i = 1:N % Add index labels
    text(x(i) + 0.03, slopes(i), num2str(i), 'FontSize', 8,'Color','b');
end
set(hf(3),"Position",[1 1 800 500])

hf(4) = figure; % comparison between metrics
subplot(131)% bias scores vs spectral slope
x = bias_scores;
y = slopes;
scatter(x,y, 'filled'); axis square
xlabel('Bias Score')
ylabel('Sq. Spectral Slope')
hold on
rsquared = fitlm(x,y).Rsquared.Ordinary;
text(min(x),mean(y),['R^2 = ',num2str(rsquared)], ...
    'FontSize',8,'Color','r');
for i = 1:N % Add index labels
    text(x(i) + 0.03*diff(xlim), y(i), num2str(i), 'FontSize', 8,'Color','b');
end

subplot(132) % bias scores vs centroids
x = bias_scores;
y = centroids;
scatter(x,y, 'filled'); axis square
xlabel('Bias Score')
ylabel('Profile centroid')
hold on
rsquared = fitlm(x,y).Rsquared.Ordinary;
text(min(x),mean(y),['R^2 = ',num2str(rsquared)], ...
    'FontSize',8,'Color','r');
for i = 1:N % Add index labels
    text(x(i) + 0.03*diff(xlim), y(i), num2str(i), 'FontSize', 8,'Color','b');
end

subplot(133) % bias scores vs KL divergence
x = bias_scores;
y = KL_divergence;
scatter(x,y, 'filled'); axis square
xlabel('Bias Score')
ylabel('KL divergence')
hold on
rsquared = fitlm(x,y).Rsquared.Ordinary;
text(min(x),mean(y),['R^2 = ',num2str(rsquared)], ...
    'FontSize',8,'Color','r');
for i = 1:N % Add index labels
    text(x(i) + 0.03*diff(xlim), y(i), num2str(i), 'FontSize', 8,'Color','b');
end
ax = gca; try; ax.YScale='log'; catch; end

set(hf(4),'Position',[1 1 800 300])
end
