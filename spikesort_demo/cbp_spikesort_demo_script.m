%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Demonstration script for CBP spike sorting.  
%%% Please see the file README.md for a description of this code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the current working directory to the directory containing this file:
% cd spikesort_demo

% Run the setup function, which sets paths and prints warnings or errors if
% there are issues detected (for example, mex/C files that need to be compiled
% for your system).
spikesort_demo_setup(pwd());

%% -----------------------------------------------------------------
% Step 0: Load raw electrode data

% Load an example data set, including raw data, the timestep, and (optionally) ground
% truth spike times.

% Simulated data: Single electrode, from: Quiroga et. al., Neural Computation,
% 16:1661-1687, 2004:
datasetName = 'Quiroga1';

% Real data: Tetrode + one ground-truth intracellular electrode, rat hippocampus,
% from: Harris et. al., J. Neurophysiology, 84:401-414, 2000:  
% datasetName = 'Harris1';

params = load_default_parameters();
[data, params] = load_raw_data(datasetName, params);

%% -----------------------------------------------------------------
% Preprocessing Step 1: Temporal filtering

% Remove low and high frequencies - purpose is to eliminate non-signal parts of the
% frequency spectrum, and enable crude identification and removal of segments
% containing spikes via local amplitude thresholding, after which the background
% noise covariance can be estimated.  In addition to filtering, the code removes the
% mean from each channel, and rescales the data (globally) to have a max abs value of
% one.  Filtering parameters:
%   - freq : range of frequencies (in Hz) for designing filter
%            Set to [] to turn off pre-filtering.
%   - type : type of filter for preprocessing. Currently supports
%            "fir1" and "butter"
%   - pad  : number of constant-value samples to pad
%   - order : order of the filter

filtdata = FilterData(data, params);

% Diagnostics for filtering: 
% Fig 1 shows filtered data.  In the next step, noise covariance will be estimated
% from below-threshold regions, which are indicated in red. There should be no spikes
% in these regions.  Fig 2 shows Fourier amplitude (effects of filtering should be
% visible.  Fig 3 shows histogram of the cross-channel magnitudes.  Below-threshold
% portion is colored red, and should look like a chi distribution with fitted
% variance (green).  If spikes appear to be included in the noise segments, reduce
% params.whitening.noise_threshold, or modify the filtering parameters in
% params.filtering, and re-run the filtering step.

%% -----------------------------------------------------------------
% Preprocessing Step 2: Estimate noise covariance and whiten

% Estimate and whiten the noise, assuming channel/time separability. This makes the
% L2-norm portion of the CBP objective into a sum of squares, simplifying the
% computation and improving computational efficiency.
%   - find "noise zones" (regions of length min_zone_len or more 
%     samples whose cross-channel L2-norm are below noise_threshold).
%   - compute noise auto-correlation function from these zones
%   - whiten each channel in time with the inverse matrix sqrt of the auto-correlation
%   - whiten across channels with inverse matrix sqrt of the covariance matrix.
% Whitening parameters:
%   - whitening.noise_threshold : rms signal magnitude (over channels) is 
%     compared to this value to determine noise regions.
%   - reg_const : a regularization constant, added to diagonal of measured covariance.

data_pp = WhitenNoise(filtdata, params);

% Diagnostics for whitening:
% Fig 4: original vs. whitened autocorrelation(s), which should be close to a delta
%   function (1 at 0, 0 elsewhere).  If not, try increasing
%   params.whitening.num_acf_lags.  If auto-corrlation is noisy, there may not be
%   enough data samples for estimation.  This can be improved by a. increasing
%   params.whitening.noise_threshold (allow more samples) b. decreasing
%   params.whitening.num_acf_lags c. decreasing params.whitening.min_zone_len (allow
%   shorter noise zones).
% Fig 5 (for multi-electrodes only): cross-channel correlation, which should be
%   diagonal. If not, a. increase params.whitening.num_acf_lags or b. increase
%   params.whitening.min_zone_len .  Note that this trades off with the quality the
%   estimates (see prev).
% Fig 1: Highlighted segments of whitened data (green) will be used to estimate
%   waveforms in the next step.  These should contain spikes (and non-highlighted
%   regions should contain background noise).  
% Fig 3, Top: Histograms of whitened channels - central portion should look
% Gaussian. Bottom: Histogram of across-channel magnitude, with magnitudes of
% highlighted segments in green.  If spikes are in noise regions, reduce params.whitening.threshold

%% -----------------------------------------------------------------
% Preprocessing Step 3: Estimate initial spike waveforms

% Initialize spike waveforms using clustering:
%  - collect data windows with L2-norm larger than params.clustering.spike_threshold
%  - align peaks of these windows
%  - Perform PCA on these segments
%  - Perform K-means clustering in the subspace of the principal components
%    accounting for params.clustering.percent_variance portion of the total variance.
% Parameters:
%  - num_waveforms : number of cells to be recovered from data.
%  - waveform_len : length (in samples) of spike waveforms.
%  - threshold : threshold used to pick spike-containing data segments (in stdevs)
%  - align_mode : During the clustering stage used to initialize spike waveforms,
%       segments of the trace are identified using a threshold. These segments are
%       then aligned based on this parameter. Supported values are:
%            maxrms : align wrt max L2 norm (RMS) value in segment
%            maxabs : align wrt max L1 norm value in segment
%            max : max signed sum across electrodes
%            min : min signed sum across electrodes
%       For example, if the spike waveforms always have large positive (negative)
%       peaks, use max (min).

[centroids, assignments, X, XProj, PCs, snippet_centers_cl] = ...
    EstimateInitialWaveforms(data_pp, params);

% For later comparisons, also compute spike times corresponding to the segments
% assigned to each cluster:
spike_times_cl = GetSpikeTimesFromAssignments(snippet_centers_cl, assignments);

VisualizeClustering(XProj, assignments, X, data_pp.nchan, ...
		    params.plotting.first_fig_num+3, ...
		    params.clustering.spike_threshold);

% Diagnostics for clustering waveform initialization: 
% At this point, waveforms of all potential cells should be identified (note that
% missed or false positive spikes are irrelevant - only the WAVEFORMS matter).  If
% not, may need to adjust params.clustering.num_waveforms and re-run the clustering
% to identify more/fewer cells.  May also wish to adjust the
% params.general.waveform_len, increasing it if the waveforms (Fig 5) are being
% chopped off, or shortening it if there is a substantial boundary region of silence.
% If you do this, you should go back and re-run starting from the whitening step,
% since the waveform_len affects the identification of noise regions.

%% -----------------------------------------------------------------
% CBP setup

% Partition data into snippets, separated by "noise zones" in which the RMS of the
% waveforms does not surpass "threshold" for at least "min_separation_len"
% consecutive samples. Snippets are constrained to have duration between
% min/max_snippet_len.  NOTE: This is for computational speedup only, so choose a
% conservative (low) threshold to avoid dropping spikes!

[snippets, breaks, snippet_lens, snippet_centers, snippet_idx] = ...
    PartitionSignal(data_pp.data, params.partition);
fprintf('Partitioned signal into %d chunks\n', length(snippets));

%*** Close all preprocessing-related figures
%figNums = params.plotting.first_fig_num+[0:2];
%close(figNums(ishghandle(figNums)));

% Fill in missing parameters
if (isempty (params.cbp.num_features))
  %num_waveforms = params.clustering.num_waveforms;
  num_waveforms = size(centroids, 2);
  params.cbp.num_features = num_waveforms;
else
  num_waveforms = params.cbp.num_features;
end

% Set initial estimates of spike waveforms to the clustering centroids
params.cbp_outer.init_features = cell(num_waveforms, 1);
params.cbp_outer.num_chan = data_pp.nchan;
for i = 1 : num_waveforms
    params.cbp_outer.init_features{i} = ...
        reshape(centroids(:, i), [], params.cbp_outer.num_chan);
end 

% Prior expected firing rates
params.cbp.firing_rate = 1e-3 .* ones(num_waveforms, 1);

% Try single-spike soln first, accepting if better than tolerance.  If most spikes
% are isolated (no overlap with other spikes), this makes the code run 
% much faster 
params.cbp.compare_greedy = false; 
params.cbp.greedy_p_value = 0; % tolerance to accept greedy soln
% params.cbp.greedy_p_value = 1 - 1e-5;

% Corr. threshold below which atoms will not be used during CBP.
% For speedup only; set to 0 to disable
params.cbp.prefilter_threshold = 0; %0.01, ... 

% -----------------------------------------------------------------
% Pick solver and reweighting parameters

% Setup CVX if needed
% addpath(fullfile(sstpath, '../cvx/'));
% cvx_setup;
% cvx_solver sedumi; 
% params.cbp.solve_fn = @cbp_cvx;
% params.cbp.solve_fn = @cbp_qcml_sumsquare;
% reweight_exp = 25 * [1 1 1];

params.cbp.solve_fn = @cbp_ecos_2norm;
reweight_exp = 1.5 * ones(1, num_waveforms);
params.cbp.lambda = reweight_exp(:); % multiplier for sparsity weight

% FIXME: Move function definition inside CBP, but leave explanation of
% reweight_exp for users.
%
% Set the reweighting function for the Iteration Reweighted L1 optimization
% for inferring the spikes. Theoretically, the new weight for a coefficient
% x should be set to -d/dz(log(P(z)))|z=x where P(z) is the prior density
% on the spike amplitudes. Here we employ a power-law distribution for 
% 0 <= x <= M with exponent=reweight_exp and offset=eps1
params.cbp.reweight_fn = cell(num_waveforms, 1);
for i = 1 : num_waveforms
    params.cbp.reweight_fn{i} = @(x) reweight_exp(i) ./ (eps + abs(x));    
end

%% -----------------------------------------------------------------
% CBP step 1: use CBP to estimate spike times of all cells

params.cbp.progress = true; % Set false if progress bar causes Java errors

starttime = tic;
[spike_times, spike_amps, recon_snippets] = ...
    SpikesortCBP(snippets, ...
                 snippet_centers, ...
                 params.cbp_outer, ...
                 params.cbp);
toc(starttime);

% Histogram of windowed norm for data, whitened data, residual
%** a mess for mulitple electrodes
if (params.general.plot_diagnostics)
  data_recon = cell(size(snippets));
  for i = 1:numel(snippets)
    data_recon{i} = snippets{i} - recon_snippets{i};
  end
  wnsnip = cellfun(@(s) windowed_norm(s', params.general.waveform_len), snippets, 'UniformOutput', false);
  wnbreak = windowed_norm(cell2mat(breaks')', params.general.waveform_len);
  wnresid = cellfun(@(s) windowed_norm(s', params.general.waveform_len), data_recon, 'UniformOutput', false);
  
  figure(params.plotting.first_fig_num+2);
  [N,X] = hist(cell2mat(wnsnip')', 100);
  % Nbreak = hist(wnbreak', X);
  Nresid = hist(cell2mat(wnresid')', X);
  chi = 2*X.*chi2pdf(X.^2, nchan*params.general.waveform_len);
  bar(X,N); set(gca,'Yscale','log'); yrg= get(gca, 'Ylim');
  hold on;
  % plot(X,(max(N)/max(Nbreak))*Nbreak,'r','LineWidth', 2);
  plot(X,(max(N)/max(Nresid))*Nresid, 'g', 'LineWidth', 2);
  plot(X,(max(N)/max(chi))*chi, 'c','LineWidth', 2);
  hold off; set(gca,'Ylim', yrg);
  legend('all snippets', 'snippet post-CBP residuals', 'expected noise (Chi)');
  title('Histogram of windowed 2-norms');
end

% Initial CBP diagnostics: (Unnecessary?)
% Residual should not have high-amplitude regions (if it does, increase
% params.cbp.firing_rates, and re-run CBP).
    
%% -----------------------------------------------------------------
% Pick amplitude thresholds and interactively visualize effect on ACorr/XCorr

% Allow this much slack in time bins of spike locations, for live updating 
% ground truth feedback
% **TODO: Why is slack parameter needed here?
% **TODO: Just make this a default value, either 30 or 40 inside relevant
% functions.  But first have to decide whether to try to reshift spike
% times to remove filtering-induced shifts, and shifts in ground truth
% data (e.g. intracellular electrode time lead).
spike_location_slack = 30; 

% Interactive selection of final thresholds for spike detection.  Top
% row of Figure 9 shows the distribution of amplitudes for each
% waveform (normalized, average spike has amplitude 1).  Bottom row
% shows the spike-conditioned raster for each cell, used to check for
% refractory violations.  Middle rows show spike-conditioned rasters
% across pairs of cells, and is used to check for dropped synchronous
% spikes (very common with clustering methods).  The figure comes up with 
% automatic values chosen.  You can adjust the
% thresholds by dragging red lines for each cell independently. Quit
% by closing the figure window.
[atgf amp_threshold] = AmplitudeThresholdGUI(spike_amps, spike_times, 'dt', data_pp.dt, 'location_slack', spike_location_slack);

% Diagnostics for CBP spike-detection: 
% Amplitudes should be around 1, Auto-correlations should show zero spikes in the
% refractory period (roughly 1-3msec), and Cross-correlations should not have a
% narrow notch around zero (common in clustering methods).

%% -----------------------------------------------------------------
% CBP Step 2: re-estimate waveforms

% Compute waveforms using regression, with interpolation (defaults to cubic spline)
nlrpoints = (params.general.waveform_len-1)/2;
waveforms = cell(size(spike_times));
for i = 1:numel(spike_times)
    % Have to shift spiketimes by 1 because existing code treats data as 1 indexed.
    sts = spike_times{i}(spike_times{i} > amp_threshold(i)) - 1;
    waveforms{i} = CalcSTA(data_pp.data', sts, [-nlrpoints nlrpoints]);
end

% Compare updated waveforms to initial estimates
if (params.general.plot_diagnostics)
    figure(params.plotting.first_fig_num+1);
    nc = ceil(sqrt(num_waveforms));
    nr = ceil(num_waveforms / nc);

    for i = 1:numel(waveforms)
        subplot(nr, nc, i);
        inith = plot(params.cbp_outer.init_features{i}, 'b');
        hold on
        finalh = plot(waveforms{i}, 'r');
        hold off
        err = norm(params.cbp_outer.init_features{i} - waveforms{i})/...
              norm(waveforms{i});
        title(sprintf('Waveform %d, Rel error=%.2f', i, err))
        legend([inith(1) finalh(1)], {'Initial', 'New'});
    end
end

% Diagnostics for CBP waveforms: 
% If recovered waveforms differ significantly from initial waveforms, then algorithm
% has not yet converged.  Execute this and go back to re-run CBP:
%     params.cbp_outer.init_features = waveforms

%**TODO: Visualize waveforms/spikes in PC space, compared to clustering result.
% Allow user to modify or increase/decrease number of waveforms.


%% -----------------------------------------------------------------
% Post-analysis: Comparison to clustering results, and to ground truth (if available)

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

if isfield(ground_truth, 'true_spike_times') && isfield(ground_truth, 'true_spike_class')
  % Since we already permuted ground truth to match clustering, this is true by definition
  best_ordering_cl = 1:length(spike_times_cl);
  best_ordering = 1:length(spike_times);

  % Evaluate clustering sorting
  [total_misses_cl, total_false_positives_cl, misses_cl, false_positives_cl] = ...
      evaluate_sorting(spike_times_cl, ground_truth.true_sp, ground_truth.spike_location_slack);
  fprintf('Clustering: %s', SortingEvaluationStr(ground_truth.true_sp, spike_times_cl, total_misses_cl, total_false_positives_cl));

  % Evaluate CBP sorting
  [total_misses, total_false_positives, prune_est_times, misses, false_positives] = ...
     EvaluateSorting(spike_times, spike_amps, ground_truth.true_sp, 'threshold', amp_threshold, 'location_slack', spike_location_slack);
  fprintf('       CBP: %s', SortingEvaluationStr(ground_truth.true_sp, prune_est_times, total_misses, total_false_positives));

end

% -----------------------------------------------------------------
% Plot various snippet subpopulations
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


%% -----------------------------------------------------------------
% Visualize true spike assignments in PC-space

if isfield(ground_truth, 'true_spike_class') && isfield(ground_truth, 'true_spike_times')
    cluster_pars = params.clustering;
    if isempty(cluster_pars.window_len), cluster_pars.window_len = params.general.waveform_len; end
    cluster_pars.align_mode = data_pp.polarity;

    Xstar = ConstructSnippetMatrix(data_pp.data, ground_truth.true_spike_times, cluster_pars);
    % Remove mean component and project onto leading PC's
    XProjstar = (Xstar - repmat(mean(Xstar, 2), 1, size(Xstar, 2)))' * PCs;
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

%% -----------------------------------------------------------------
% Get greedy spike matches and plot RoC-style
% NB: Much faster if mex greedymatchtimes.c is compiled
PlotCBPROC(spike_times, spike_amps, ground_truth.true_sp, ground_truth.spike_location_slack);
