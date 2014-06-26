%% Test script for running CBP spike sorting algorithm
%
% Chaitanya Ekanadham 12/26/2012
% Lab for Computational Vision, NYU & HHMI
%
% Updated by Peter H. Li, fall/winter 2013
%
% This script demonstrates our spike sorting method which relies on a
% sparse inference method known as Continuous Basis Pursuit (CBP).
%
% For more information see our publication:
%
% C. Ekanadham, D. Tranchina and E. P. Simoncelli. A blind sparse
%   deconvolution method for neural spike identification. Published in
%   Adv. Neural Information Processing Systems (NIPS*11), vol.24 
%   pp. 1440--1448, Dec 2011. 
%
% Available at: http://www.cns.nyu.edu/pub/lcv/Ekanadham11c-preprint.pdf


%% Setup
% Navigate Matlab to the demo directory, e.g.:
%   cd spikesort_demo
%
% Then run the setup function below from the demo directory.
%
% The setup function will print warnings or errors if there are issues 
% detected with setup, for example if there are mex files that need to be
% compiled for your system.
%
CBProotPath = '/home/peterli/MATLAB/matlab-standard/private/phli/CBP';
% CBProotPath = '/Users/eero/matlab/CBP';
cd(CBProotPath)
cd CBPSpikesortDemoPackage/spikesort_demo/
spikesort_demo_setup(pwd());

%% Load data, default parameters

% The input file path should be in "data_filename"
% This .mat file should contain:
%   data: channel x time matrix of trace values
%   dt : timestep (in seconds)
%
% There are a few parameters that should be set based on the data set:
%
%   filt_type : type of filter for preprocessing. Currently supports
%               "fir1" and "butter"
%
%   filt_freq_range : boundaries of freq band (in Hz) for designing filter
%                     Set as [] to avoid any pre-filtering.
%
%   filt_order : order of the filter
%   
%   noise_threshold : threshold applied to root-mean-squared signal value
%                     (taken across electrodes) to approximately discern
%                     noise from signal. This will be used to estimate
%                     noise statistics and design a whitening filter
%                     The threshold will be visualized over the signal
%                     during preprocessing. It is recommended you try a few
%                     values to determine how to reasonably estimate the
%                     noise.
%
%   whiten_trace : whether or not to whiten w.r.t. estimated noise
%                  statistics. This uses noise_threshold.
%
%   align_mode : During the clustering stage used to initialize spike
%                waveforms, segments of the trace are identified using a
%                threshold. These segments are then aligned based on this
%                parameter. Supported values are:
%                   maxrms : align wrt max L2 norm (RMS) value in segment
%                   maxabs : "" "" max L1 norm value in segment
%                   max : max signed sum across electrodes
%                   min : min "" "" "" 
%
%                For example, if the spike waveforms always have large
%                positive (negative) peaks, use max (min).
%
%   NUM_WAVEFORMS : estimated number of cells.
%
%   WINDOW_LEN : estimated length (in samples) of spike waveforms.

% Here are two basic example datasets.  There are additional options in the
% example_data directory.

switch 0 
    case 0
        % Quiroga simulated data examples
        data_filename = 'C_Easy1_noise015.mat';
        
    case 1
        % Harris data example
        %     data_filename = 'example_data/harris_d533101_v2.mat';
        data_filename = 'harris_d533101.mat';
end

% Load data
data = load(data_filename, 'data', 'dt');
data.filename = data_filename;
clear data_filename;

% Load default parameters
params = load_default_parameters();
[data, params] = load_data_defaults(data, params);


%% Preprocess the trace

% Filtering:
%   - pad with constant values (to avoid border effects when filtering)
%   - filter with band/highpass FIR or Butterworth filter
%   - remove mean from each channel
%   - scale globally (across all channels) to be in [-1, 1]


% Noise-whitening (assumes space/time separable noise):
%   - estimate "noise zones" (regions of length min_zone_len or more 
%     samples whose Lp-norm are under noise_threshold).
%   - compute noise ACF from these zones
%   - whiten each channel in time by taking inverse matrix sqrt of ACF
%   - whiten across channels by left-multiplying each time slice 
%     by inverse matrix sqrt of the covariance matrix).

filtdata = FilterData(data, params.filtering);
data_pp = PreprocessTrace(filtdata, params, 'plot', false);
clear filtdata; % Probably don't need to keep this around

% ADD : instead of L_inf norm, use an estimate of the waveform length
%       take a sliding L2 norm and cut the histogram based on this.


% ADD : histogram of samples and plot noise threshold (color-coded)
%       make noise_pars.num_acf_lags in units of seconds
%       General p-norm (p >= 1)

% 1. Noise threshold should selected to isolate the background activity
% from foreground spikes. (In Figure 2, red zones should look like noise)
% If not, adjust noise_pars.threshold accordingly.
%
% 2. If original ACF is not going to zero, increase noise_pars.num_acf_lags.
% If original ACF is noisy, there may not be enough data samples to
% robustly estimate it. You can increase the amount of data in several ways:
%   a. increase noise_pars.threshold ( to allow more samples)
%   b. decrease noise_pars.num_acf_lags
%   c. decrease noise_pars.min_zone_len (to allow shorter noise zones)

% 3. Spatial covariance matrix of whitened signal will be the identity matrix
% and ideally whitened ACF should be a delta function. In practice, this may
% not be the case because of (1) non-separability of noise in space/time,
% and/or (2) finite temporal window for estimating correlations. 
% If whitened temporal/spatial correlations are not white, try the following:
%   a. increase noise_pars.num_acf_lags
%   b. increase noise_pars.min_zone_len
% Note that this trades off with the quality the estimates (see #2 above).

% 4. In Figure 142, distribution of noise samples (blue) should be reasonably fit
%    by Generalized Gaussian distribution (red). We assume as default an
%    exponent of 2 (Gaussian). If this is not a good match, consider
%    changing the exponent accordingly (noise_pars.p_norm).


%% Plot preprocessed trace
% TODO: Replace with scrolling GUI


%% Windowed 2-Norm histograms
% TODO: Replace with higher level wrapping call taking just data and
% params.general
wndata    = windowed_norm(data.data,    params.general.waveform_len);
wndata_pp = windowed_norm(data_pp.data, params.general.waveform_len);
subplot(1,2,1);
hist(wndata', 100);
title('Windowed 2-Norm of original data');
subplot(1,2,2);
hist(wndata_pp', 100);
title('Windowed 2-Norm of filtered and whitened data');

clear wndata wndata_pp;


%% Use PCA + K-means to estimate initial waveforms
%  (and benchmark spike sorting performance using clustering)

[centroids, assignments, X, XProj, PCs, snippet_centers_cl] = ...
    EstimateInitialWaveforms(data_pp, params);
spike_times_cl = GetSpikeTimesFromAssignments(snippet_centers_cl, assignments);


%% Visualization of clustering

% Plot the putative spikes w.r.t. the 2 leading principal components.
% PC's are computed across all (aligned) windows which pass the threshold 
% test. K-means clustering is performed using the PC's accounting for
% cluster_pars.percent_variance portion of the total variance.

VisualizeClustering(XProj, assignments, X, data_pp.nchan);


%% For efficiency, chop signal into chunks for parallelization for CBP

% Data is partitioned by identifying "noise zones" in which the RMS of the
% signal does not surpass "threshold" for "min_separation_len" consecutive
% samples. Snippets are also constrained to have length between
% min/max_snippet_len.
%
% NOTE: This is for computational speedup only, 
%       so choose a very conservative (low) threshold to avoid missing
%       spikes.

% ADD : choose partion_pars.threshold based on sliding Lp-norm.
% The threshold should be the minimum window Lp-norm containing any overlap
% of any pair of waveforms. Display a histogram of sliding Lp-norms with
% this automated choice for the user.

data_rms = sqrt(sum(data_pp.data .^ 2, 1));% root-mean-squared across electrodes
%cluster_threshold = 4 * median(data_rms) ./ 0.6745; % robust
threshold = 4 * std(data_rms);

partition_pars = struct('threshold', threshold, ...
                        'smooth_len', 1, ...
                        'min_separation_len', ...
                        floor(params.general.waveform_len / 2), ...
                        'min_snippet_len', ...
                        params.general.waveform_len, ...
                        'max_snippet_len', ...
                        1001, ... % not enforced, only warnings
                        'min_pad_size', 5);
                    
[snippets, breaks, snippet_lens, snippet_centers, snippet_idx] = ...
    PartitionSignal(data_pp.data, partition_pars);
fprintf('Chopped up signal into %d chunks\n', length(snippets));

clear data_rms threshold

% MAKE TIME PARS IN UNITS OF SECONDS?


%% Plot windowed norm of snippets versus "silences"
wnsnip = cellfun(@(s) windowed_norm(s', params.general.waveform_len), snippets, 'UniformOutput', false);
wnbreak = windowed_norm(cell2mat(breaks')', params.general.waveform_len);
figure();
subplot(1,2,1);
hist(cell2mat(wnsnip')', 100);
title('Windowed 2-Norm of snippets');
subplot(1,2,2);
hist(wnbreak', 100);
title('Windowed 2-Norm of putative silences')


%% Visualization of partition of the voltage trace

% Just use scrolling GUI.  Maybe keep the histogram of snippet lengths?


%% CBP setup
% Should be able to leave most of these defaults
num_waveforms = size(centroids,2);

% The polar_1D version seemed to be cutting things down too much...
adjust_wfsize_fn = @(w) w; %polar_1D_adjust_wfsize(w, 0.1, 0.025, 301), ...

% cbp_outer_pars are parameters for learning the waveforms.
cbp_outer_pars = struct( ...
    'num_iterations', 2e2, ... % number of learning iterations
	'batch_size', 125, ... % batch size for learning
	'step_size', 5e-2, ... % step size for updating waveform shapes
	'step_size_decay_factor', 1, ... % annealing
	'plotevery',1, ... % plot interval
	'stop_on_increase', false, ... % stop when objective function increases
	'check_coeff_mtx', true, ... % sanity check (true to be safe)
	'adjust_wfsize', ... % called each iteration to adjust waveform size
        adjust_wfsize_fn, ...
	'rescale_flag', false, ... % always FALSE 
	'renormalize_features', false, ... % always FALSE
	'reestimate_priors', false, ... % always FALSE
    'CoeffMtx_fn', @polar_1D_sp_cnv_mtx, ... % convolves spikes w/waveforms
    'plot_every', 1 ... % plotting frequency    
);

% Set initial estimates of spike waveforms to the clustering centroids
cbp_outer_pars.init_features = cell(size(centroids, 2), 1);
cbp_outer_pars.num_chan = size(data_pp.nchan, 1);
for i = 1 : num_waveforms
    cbp_outer_pars.init_features{i} = ...
        reshape(centroids(:, i), [], cbp_outer_pars.num_chan);
end 

% cbp_pars are parameters for doing sparse inference.
% ADD : incorporate estimates of prior firing rates from clustering result
% ADD : pnorm of noise (determined above)
cbp_pars = struct ( ...
    'noise_sigma',  data_pp.noise_sigma, ... % Optimization parameters
    'firing_rates', 1e-3 .* ones(num_waveforms, 1), ... % prior firing rate
    'cbp_core_fn', @polar_1D_cbp_core, ... % CBP core interpolation
    'solve_fn', @cbp_ecos_2norm, ... % Optimization solver function
    'debug_mode', false, ... % debug mode
    'num_reweights', 1e3, ... % MAX number of IRL1 iterations
    'magnitude_threshold', 1e-2, ... % amplitude threshold for deleting spikes
    'parfor_chunk_size', Inf, ... % parallelization chunk size
    'num_features', num_waveforms ... % number of "cells"
);


%% CBP parameters that user should pick

% For picking template delta
cbp_pars.accuracy = 0.1;

% Try single-spike soln first
cbp_pars.compare_greedy = false; 
cbp_pars.greedy_p_value = 0;
% cbp_pars.greedy_p_value = 1 - 1e-5; % tol. to accept initial greedy soln

% Corr. threshold below which atoms will not be used during CBP.
% For speedup only; set to 0 to disable
cbp_pars.prefilter_threshold = 0; %0.01, ... 


%% Pick solver and reweighting parameters

% Setup CVX if needed
% addpath(fullfile(sstpath, '../cvx/'));
% cvx_setup;
% cvx_solver sedumi; % Can also use ECOS if using CVX2 and ECOS shim installed.

cbp_pars.solve_fn = @cbp_ecos_2norm;
reweight_exp = 1.5 * [1 1 1];
% cbp_pars.solve_fn = @cbp_cvx;
% cbp_pars.solve_fn = @cbp_qcml_sumsquare;
% reweight_exp = 25 * [1 1 1];

cbp_pars.lambda = reweight_exp(:); % multiplier for sparsity weight

% FIXME: Move function definition inside CBP, but leave explanation of
% reweight_exp for users.
%
% Set the reweighting function for the Iteration Reweighted L1 optimization
% for inferring the spikes. Theoretically, the new weight for a coefficient
% x should be set to -d/dz(log(P(z)))|z=x where P(z) is the prior density
% on the spike amplitudes. Here we employ a power-law distribution for 
% 0 <= x <= M with exponent=reweight_exp and offset=eps
cbp_pars.reweight_fn = cell(num_waveforms, 1);
for i = 1 : num_waveforms
    cbp_pars.reweight_fn{i} = @(x) reweight_exp(i) ./ (eps + abs(x));    
end



%% CBP run
% matlabpool open
starttime = tic;
cbp_pars.progress = true; % Set false if having Java errors from progress bar
[spike_times, spike_amps, recon_snippets] = ...
    SpikesortCBP(snippets, ...
                 snippet_centers, ...
                 cbp_outer_pars, ...
                 cbp_pars);
toc(starttime);

% Postprocess params
% TODO: Just make this a default value, either 30 or 40 inside relevant
% functions.  But first have to decide whether to try to reshift spike
% times to remove consistent bias; tricky as depends on how ground truth is
% defined, e.g. simulation versus intracellular electrode.
spike_location_slack = 30; % For live updating ground truth feedback


%% Pick amplitude thresholds and visualize effect on ACorr/XCorr
% NB: Much faster if mex trialevents.c is compiled

% Clustering XCorr and ACorr plots
% Optional
% spike_times_cl = GetSpikeTimesFromAssignments(snippet_centers_cl, assignments);
% f = figure();
% % ACorrs
% for i = 1:NUM_WAVEFORMS
%     subplot(NUM_WAVEFORMS, NUM_WAVEFORMS, sub2ind([NUM_WAVEFORMS NUM_WAVEFORMS], i, NUM_WAVEFORMS));
%     psthacorr(spike_times_cl{i}.*dt)
% end
% 
% % XCorrs
% for i = 1:NUM_WAVEFORMS
%     for j = i+1 : NUM_WAVEFORMS
%         subplot(NUM_WAVEFORMS, NUM_WAVEFORMS, sub2ind([NUM_WAVEFORMS NUM_WAVEFORMS], j, i));
%         psthxcorr(spike_times_cl{i}.*dt, spike_times_cl{j}.*dt)
%     end
% end


% True spike XCorr and ACorr plots (TST defined after load_data above)
% Optional
% ntruth = length(true_sp);
% figure();
% % ACorrs
% for i = 1:ntruth
%     sp = true_sp{i};
%     if isempty(sp), continue; end
%     
%     subplot(ntruth, ntruth, sub2ind([ntruth ntruth], i, ntruth));
%     psthacorr(sp.*dt)
% end
% % XCorrs
% for i = 1:ntruth
%     spi = true_sp{i};
%     if isempty(spi), continue; end
% 
%     for j = i+1 : ntruth
%         spj = true_sp{j};
%         if isempty(spj), continue; end
% 
%         subplot(ntruth, ntruth, sub2ind([ntruth ntruth], j, i));
%         psthxcorr(spi.*dt, spj.*dt)
%     end
% end

[atgf amp_threshold] = AmplitudeThresholdGUI(spike_amps, spike_times, 'dt', data_pp.dt, 'location_slack', spike_location_slack);


%% Histogram of windowed norm for data, whitened data, residual

data_recon = cell(size(snippets));
for i = 1:numel(snippets)
    data_recon{i} = snippets{i} - recon_snippets{i};
end
wnresid = cellfun(@(s) windowed_norm(s', params.general.waveform_len), data_recon, 'UniformOutput', false);

figure();
sanesubplot(1, 3, {1 1});
hist(cell2mat(wnsnip')', 100);
title('Windowed 2-Norm of snippets')
sanesubplot(1, 3, {1 2});
hist(wnbreak', 100);
title('Windowed 2-Norm of putative silences')
sanesubplot(1, 3, {1 3});
hist(cell2mat(wnresid')', 100);
title('Windowed 2-Norm of snippet residuals after CBP')


%% Calculate new waveform estimates with interpolation (defaults to cubic spline)
nlrpoints = (params.general.waveform_len-1)/2;
waveforms = cell(size(spike_times));
for i = 1:numel(spike_times)
    % Have to shift spiketimes by 1 because existing code treats data as 1 indexed.
    sts = spike_times{i}(spike_times{i} > amp_threshold(i)) - 1;
    waveforms{i} = CalcSTA(data_pp.data', sts, [-nlrpoints nlrpoints]);
end

% Show diff with initial estimates
figure();
nc = ceil(sqrt(num_waveforms));
nr = ceil(num_waveforms / nc);
for i = 1:numel(waveforms)
    subplot(nr, nc, i);
    inith = plot(cbp_outer_pars.init_features{i}, 'b');
    hold on
    finalh = plot(waveforms{i}, 'r');

    err = norm(cbp_outer_pars.init_features{i} - waveforms{i});
    title(sprintf('Norm of diff over norm of final: %.2f', ...
        err / norm(waveforms{i})));
	legend([inith(1) finalh(1)], {'Initial', 'Final'});
end


%% Set init_features to new waveforms and go back up to rerun CBP with new waveforms estimates.
% Optional
% cbp_outer_pars.init_features = waveforms;


%% Info on true spikes available?
ground_truth = load(data_pp.filename, 'true_spike_times', 'true_spike_class', 'dt');
ground_truth.filename = data_pp.filename;

% Acceptable slack for considering two spikes a match.  In units of samples.
% Currently two-sided, but this should probably be changed.
ground_truth.spike_location_slack = 30;

if isfield(ground_truth, 'true_spike_times') && isfield(ground_truth, 'true_spike_class')
    % Reformat as 1-cellarr per cell of spike times.
    ground_truth.true_sp = GetSpikeTimesFromAssignments(ground_truth.true_spike_times, ground_truth.true_spike_class);
    
    % Reorder to match cell numbering from clustering.
    [ground_truth.true_sp, ...
     ground_truth.true_spike_times, ... 
     ground_truth.true_spike_class] = ReorderCells( ...
        ground_truth.true_sp, spike_times_cl, ground_truth.spike_location_slack);
end


%% Plot various snippet subpopulations
[est_matches true_matches] = GreedyMatchTimes(spike_times, ground_truth.true_sp, ground_truth.spike_location_slack);

% Complete misses (missed even at 0 threshold)
completemisses = cell(size(ground_truth.true_sp));
for i = 1:length(ground_truth.true_sp)
    completemisses{i} = ground_truth.true_sp{i}(true_matches{i} == 0);
end

% All FPs (with 0 threshold)
allfps = cell(size(spike_times));
for i = 1:length(spike_times)
    allfps{i} = spike_times{i}(est_matches{i} == 0);
end

% For example; more complicated selections of spikes like misses or FPs
% with given threshold can also be calculated
% TODO PHLI: automate this.
% desiredspiketimes = cell2mat(allfps);         % Pick one or the other...
desiredspiketimes = cell2mat(completemisses); % Pick one or the other...
% desiredspiketimes = cell2mat(spike_times);

% Now find and plot the relevant snippets
snipindices = FindSnippets(desiredspiketimes, snippet_centers, snippets); % snippets just given as shorthand for calculating widths
ScrollSnippets(snippets, snippet_centers, ...
    'snipindices',  unique(snipindices(snipindices > 0)),  ...
    'cbp',          spike_times,        ...
    'cbpamp',       spike_amps,         ...
... %     'cbpampthresh', amp_threshold,      ... % Could use amp_threshold if we used that to pick snippets...
    'clust',        spike_times_cl,     ...
... %     'recons',       recon_snippets,     ...
    'true',         ground_truth.true_sp);


%% Visualize true spike assignments in PC-space

if isfield(ground_truth, 'true_spike_class') && isfield(ground_truth, 'true_spike_times')
    cluster_pars = params.clustering;
    if isempty(cluster_pars.window_len), cluster_pars.window_len = params.general.waveform_len; end
    cluster_pars.align_mode = data_pp.polarity;

    Xstar = ConstructSnippetMatrix(data_pp.data, ground_truth.true_spike_times, cluster_pars);
    % Remove mean component and project onto leading PC's
    XProjstar = (Xstar - repmat(mean(X, 2), 1, size(Xstar, 2)))' * PCs;
    VisualizeClustering(XProjstar, ground_truth.true_spike_class, Xstar, size(data_pp, 1), ...
                        figure(8), figure(9), '.');    
    
    % Plot the putative spikes w.r.t. the 2 leading principal components.
    % PC's are computed across all (aligned) windows which pass the threshold
    % test. K-means clustering is performed using the PC's accounting for
    % cluster_pars.percent_variance portion of the total variance.
    figure(8); title('True clustering');
end

% ADD NEW% Plot the putative spikes w.r.t. the 2 leading principal components.
% PC's are computed across all (aligned) windows which pass the threshold 
% test. K-means clustering is performed using the PC's accounting for
% cluster_pars.percent_variance portion of the total variance.
% FIGURE OF HISTOGRAM OF WHITENEDE RMS SAMPLES cluster sthreshold
% ADD multiple PC plots (optional) 

% What to check:
%
% 1. Adjust cluster_threshold to properly separate background activity from
%   spike data (shold cleanly separate histogram in Fig 200).
%
% 2. Adjust NUM_WAVEFORMS so that separate clusters in Fig 6 are identified
%    with separate colors (i.e. waveforms in Fig 7 should all have distinct
%    shapes).


%% Get greedy spike matches and plot RoC-style
% NB: Much faster if mex greedymatchtimes.c is compiled
PlotCBPROC(spike_times, spike_amps, ground_truth.true_sp, ground_truth.spike_location_slack);


%% Evaluation (ONLY if truth is available)

if isfield(ground_truth, 'true_spike_times') && isfield(ground_truth, 'true_spike_class')
    % Since we already permuted ground truth to match clustering, this is true
    % by definition
    best_ordering_cl = 1:length(spike_times_cl);
    best_ordering = 1:length(spike_times);
    
    % Evaluate clustering sorting
    [total_misses_cl, total_false_positives_cl, misses_cl, ...
        false_positives_cl] = ...
        evaluate_sorting(spike_times_cl, ground_truth.true_sp, ground_truth.spike_location_slack);
    fprintf('Clust %s', SortingEvaluationStr(ground_truth.true_sp, spike_times_cl, total_misses_cl, total_false_positives_cl));
    
    % Evaluate CBP sorting
    [total_misses, total_false_positives, prune_est_times, misses, ...
        false_positives] = ...
        EvaluateSorting(spike_times, spike_amps, ground_truth.true_sp, 'threshold', amp_threshold, 'location_slack', spike_location_slack);
    fprintf('  CBP %s', SortingEvaluationStr(ground_truth.true_sp, prune_est_times, total_misses, total_false_positives));
    
end


%% Things to add

% ADD : visualization of CBP spikes in PC-space alongside clustering
