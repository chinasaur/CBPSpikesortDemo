function [centroids, assignments, X, XProj, PCs, peak_idx] = ...
    EstimateInitialWaveforms(datastruct, params)

% Use PCA + K-means clustering to get an intial estimate of the waveform
% shapes. Also used to benchmark performance of clustering for spike
% sorting.
%

% TODO: Describe params

data = datastruct.data;
pars = params.clustering;
gen_pars = params.general;
pars.align_mode = datastruct.polarity;

if isempty(pars.window_len)
    pars.window_len = gen_pars.waveform_len;
end

if isempty(pars.peak_len)
    pars.peak_len = floor(pars.window_len / 2);
end


% Root-mean-squared (across channels) of data.
data_rms = sqrt(sum(data .^ 2, 1));

% Clustering parameters
%threshold = 4 * median(data_rms) ./ 0.6745; % robust
threshold = 4 * std(data_rms);

% Identify time indices of candidate peaks.
peak_idx = FindPeaks(data_rms(:), threshold, pars);
fprintf('Found %d candidates.\n', length(peak_idx));

% Construct a matrix with these windows, upsampled and aligned.
X = ConstructSnippetMatrix(data, peak_idx, pars);

[PCs, XProj] = TruncatePCs(X, pars.percent_variance);

% Do K-means clustering
assignments = DoKMeans(XProj, pars.num_waveforms);
centroids = GetCentroids(X, assignments);

% Put them in a canonical order (according to increasing 2norm);
[centroids, clperm] = OrderWaveformsByNorm(centroids);
assignments = PermuteAssignments(assignments, clperm);


%% Subroutines

% Find peaks (i.e. values greater than any other value within
% pars.peak_len samples).
function peak_idx = FindPeaks(data_rms, threshold, pars)
if (size(data_rms, 2) > 1)
    error('FindPeaks: can only find peaks in a vectorized signal!');
end
peak_idx = data_rms > threshold;

% Don't include borders
peak_idx(1 : pars.window_len) = false;
peak_idx(end - pars.window_len : end) = false;
for i = -pars.peak_len : pars.peak_len
    peak_idx = peak_idx & data_rms >= circshift(data_rms, i);
end
peak_idx = find(peak_idx);


% Get principal components accounting for desired percent variance
% Return the PCs as well as the projections of the data on to these PCs
% (PX)
function [PCs, XProj] = TruncatePCs(X, percent_variance)
fprintf('Doing PCA...');
% Get PCs
[PCs, Xproj, latent] = princomp(X');
[latent sorted_idx] = sort(latent, 'descend');
PCs = PCs(:,sorted_idx);
Xproj = Xproj(:, sorted_idx);

% Figure out how many PCs we need to account for
% the desired percent of total variance
cutoff = find(cumsum(latent) ./ sum(latent) * 100 > ...
    percent_variance, 1);
npcs = max(2, cutoff);
fprintf('%d PCs account for %.2f percent variance\n', ...
    npcs, percent_variance);
% Project onto leading npcs PCs
PC = PCs(:, 1 : npcs);
% Project on to PCs
XProj = Xproj(:, 1 : npcs);



% Do K-means clustering on the PC representation of the snippets
function assignments = DoKMeans(Xproj, nwaveforms)
fprintf('Clustering using k-means...\n');

% Number of times to try with random initialization
num_reps = 25;

% Use default K-means settings for now.
distance_mode = 'sqEuclidean';
start_mode = 'sample'; % centroid initialization = random sample
empty_mode = 'error'; % throw error when clusters are empty
opts = statset('MaxIter', 1e3);
% Run K-means
assignments = kmeans(Xproj, nwaveforms,...
    'Replicates', num_reps, ...
    'Distance', distance_mode, ...
    'Start', start_mode, ...
    'EmptyAction', empty_mode, ...
    'Options', opts);
fprintf('Done.\n');


% Get centroids from the assignments
function centroids = GetCentroids(X, assignments)
N = max(assignments);
centroids = zeros(size(X, 1), N);
for i = 1 : N
    idx = (assignments == i);
    centroids(:, i) = mean(X(:, idx), 2);
end

