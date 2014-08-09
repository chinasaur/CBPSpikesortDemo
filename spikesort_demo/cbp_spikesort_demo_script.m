%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Demonstration script for spike sorting using Continuous Basis Pursuit (CBP).  
%%% Please see the file README.md for a description of this code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the current working directory to the directory containing this file:
% cd spikesort_demo

% Run the setup function, which sets paths and prints warnings or errors if
% there are issues detected (for example, mex/C files that need to be compiled
% for your system).
spikesort_demo_setup(pwd());

%% ----------------------------------------------------------------------------------
% Step 0: Load raw electrode data

% Load an example data set, including raw data, the timestep, and (optionally) ground
% truth spike times.

% Simulated data example: Single electrode, from: Quiroga et. al., Neural
% Computation, 16:1661-1687, 2004
datasetName = 'Quiroga1';

% Real data example: Tetrode + one ground-truth intracellular electrode, rat
% hippocampus, from: Harris et. al., J. Neurophysiology, 84:401-414, 2000 
% datasetName = 'Harris1';

params = load_default_parameters();
[data, params] = load_raw_data(datasetName, params);

% Fig 1a shows the raw data.  
% Fig 2a plots the Fourier amplitude (averaged across channels).

%% ----------------------------------------------------------------------------------
% Preprocessing Step 1: Temporal filtering

% Remove low and high frequencies - purpose is to eliminate non-signal parts of the
% frequency spectrum, and enable crude removal of segments containing spikes via
% local amplitude thresholding, so that background noise covariance can be estimated.
% In addition to filtering, the code removes the mean from each channel, and rescales
% the data (globally) to have a max abs value of one.  params.filtering includes:
%   - freq : range of frequencies (in Hz) for designing filter
%            Set to [] to turn off pre-filtering.
%   - type : type of filter for preprocessing. Currently supports
%            "fir1" and "butter"
%   - pad  : number of constant-value samples to pad
%   - order : order of the filter

filtdata = FilterData(data, params);

% Diagnostics for filtering: 
% Fig 1b shows filtered data.  In the next step, noise covariance will be estimated
% from below-threshold regions, which are indicated in red. There should be no spikes
% in these regions.  Fig 2 shows Fourier amplitude (effects of filtering should be
% visible).  Fig 3 shows histogram of the cross-channel magnitudes.  Below-threshold
% portion is colored red, and should look like a chi distribution with fitted
% variance (green curve).  If spikes appear to be included in the noise segments,
% reduce params.whitening.noise_threshold before proceeding, or modify the filtering
% parameters in params.filtering, and re-run the filtering step.

%% ----------------------------------------------------------------------------------
% Preprocessing Step 2: Estimate noise covariance and whiten

% Estimate and whiten the noise, assuming channel/time separability. This makes the
% L2-norm portion of the CBP objective into a sum of squares, simplifying the
% computation and improving computational efficiency.
%   - find "noise zones"
%   - compute noise auto-correlation function from these zones
%   - whiten each channel in time with the inverse matrix sqrt of the auto-correlation
%   - whiten across channels with inverse matrix sqrt of the covariance matrix.
% params.whitening includes:
%   - noise_threshold: noise zones must have cross-channel L2-norm less than this.
%   - min_zone_len: noise zones must have duration of at least this many samples.
%   - num_acf_lags: number of samples over which auto-correlation is estimated.

data_pp = WhitenNoise(filtdata, params);

% Diagnostics for whitening:
% Fig 4: original vs. whitened autocorrelation(s), should be close to a delta
%   function (1 at 0, 0 elsewhere).  If not, try increasing
%   params.whitening.num_acf_lags.  If auto-correlation is noisy, there may not be
%   enough data samples for estimation.  This can be improved by a. increasing
%   params.whitening.noise_threshold (allow more samples) b. decreasing
%   params.whitening.num_acf_lags c. decreasing params.whitening.min_zone_len (allow
%   shorter noise zones).
% Fig 5 (multi-electrodes only): cross-channel correlation, should look like the
%   identity matrix. If not, a. increase params.whitening.num_acf_lags or b. increase
%   params.whitening.min_zone_len .  Note that this trades off with the quality of
%   the estimates (see prev).
% Fig 1: Highlighted segments of whitened data (green) will be used to estimate
%   waveforms in the next step.  These should contain spikes (and non-highlighted
%   regions should contain background noise).  Don't worry about getting all the
%   spikes: these are only used to initialize the waveforms!
% Fig 3, Top: Histograms of whitened channels - central portion should look
%   Gaussian. Bottom: Histogram of across-channel magnitude, with magnitudes of
%   highlighted segments in green.  If lots of spikes are in noise regions, reduce
%   params.whitening.noise_threshold

%% ----------------------------------------------------------------------------------
% Preprocessing Step 3: Estimate initial spike waveforms

% Initialize spike waveforms, using clustering:
%  - collect data segments with L2-norm larger than params.clustering.spike_threshold
%  - align peaks of waveforms within these segments
%  - Perform PCA on these segments, select a subspace containing desired percent of variance
%  - Perform K-means clustering in this subspace
% params.clustering includes:
%  - num_waveforms : number of cells to be recovered
%  - spike_threshold : threshold used to pick spike-containing data segments (in stdevs)
%  - percent_variance : used to determine number of principal components to use for clustering

[centroids, assignments, X, XProj, PCs, segment_centers_cl] = ...
    EstimateInitialWaveforms(data_pp, params);

init_waveforms = waveformMat2Cell(centroids, params.general.waveform_len, ...
                                  data_pp.nchan, params.clustering.num_waveforms);

if (params.general.plot_diagnostics)
  VisualizeClustering(XProj, assignments, X, data_pp.nchan, ...
                      params.plotting.first_fig_num+3, ...
                      params.clustering.spike_threshold);
end

% For later comparisons, also compute spike times corresponding to the segments
% assigned to each cluster:
spike_times_cl = GetSpikeTimesFromAssignments(segment_centers_cl, assignments);

% Diagnostics for waveform initialization: 
% Fig 5 shows the data segments projected onto the first two principal components,
% and the identified clusters.  Fig 6 shows the waveforms associated with each
% cluster.  The visualization function also prints out a table of distances between
% waveforms, and each of their distances to the origin (i.e., their norm).  All
% distances are relative to the noise amplitude.  These numbers provide some
% indication of how likely it is that waveforms could be confused with each other, or
% with background noise.

% At this point, waveforms of all potential cells should be identified (NOTE:
% spike identification errors are irrelevant - only the WAVEFORMS matter).  If
% not, may need to adjust params.clustering.num_waveforms and re-run the clustering
% to identify more/fewer cells.  May also wish to adjust the
% params.general.waveform_len, increasing it if the waveforms (Fig 5) are being
% chopped off, or shortening it if there is a substantial boundary region of silence.
% If you do this, you should go back and re-run starting from the whitening step,
% since the waveform_len affects the identification of noise regions.

%% -----------------------------------------------------------------------------------
% CBP preprocessing

closeIfOpen(params.plotting.first_fig_num+[1,3]); % close Fourier and clustering figures

% To speed up computation, partition data into "snippets", which will be processed
% independently. Snippets have duration between min/max_snippet_len and are separated
% by "noise zones" in which the RMS of the waveforms does not surpass "threshold" for
% at least "min_separation_len" consecutive samples. Choose a conservative (low)
% threshold to avoid dropping spikes!

[snippets, breaks, snippet_lens, snippet_centers, snippet_idx] = ...
    PartitionSignal(data_pp.data, params.partition);


%% ----------------------------------------------------------------------------------
% CBP step 1: use CBP to estimate spike times

% Turn off the Java progress bar if it causes errors
% params.cbp.progress = false; 

[spike_times, spike_amps, recon_snippets] = ...
    SpikesortCBP(snippets, snippet_centers, init_waveforms, params.cbp_outer, params.cbp);

if (params.general.plot_diagnostics)
    DisplaySortedSpikes(data_pp, spike_times, spike_amps, init_waveforms, ...
                        snippets, recon_snippets, params, 'CBP results');
end

% Diagnostics for CBP results:
% Fig1: visually compare whitened data, recovered spikes
% Fig2: residual histograms (raw, and cross-channel magnitudes) - compare to Fig3

%% ----------------------------------------------------------------------------------
% CBP step 2: Identify spikes by thresholding amplitudes of each cell

% Calculate default thresholds.  This is done by fitting the amplitude
% distribution (using a Gaussian kernel density estimator) and then choosing the
% largest local minimum.
amp_thresholds = cellfun(@(sa) ComputeKDEThreshold(sa, params.amplitude), spike_amps);

% Visualize thresholds (only works if image processing toolbox is installed).
if (params.general.plot_diagnostics)
    AmplitudeThresholdGUI(spike_amps, spike_times, amp_thresholds,    ...
                          'dt', data_pp.dt, ...
                          'f', figure(params.plotting.first_fig_num+6),    ...
                          'wfnorms', cellfun(@(wf) norm(wf), init_waveforms), ...
                          'location_slack', params.postproc.spike_location_slack);
end

% Fig7 allows interactive adjustment of waveform amplitudes, while visualizing effect
% on ACorr/XCorr.  Top row shows amplitude distribution (typical spikes should have
% amplitude 1), with thresholds indicated by red lines.  These can be dragged right
% or left.  Bottom row shows spike train autocorrelation that would result from
% chosen threshold, and can be examined for refractory violations.  Middle rows show
% spike train cross-correlations across pairs of cells, and can be examined for
% dropped synchronous spikes (very common with clustering methods).  Click the
% "Use waveforms" button to proceed with adjusted values.  Click the "Revert"
% button to revert to the automatically-chosen default values.

%% ----------------------------------------------------------------------------------
% CBP Step 3: Re-estimate waveforms

% Compute waveforms using regression, with interpolation (defaults to cubic spline)
nlrpoints = (params.general.waveform_len-1)/2;
waveforms = cell(size(spike_times));
for i = 1:numel(spike_times)
    sts = spike_times{i}(spike_amps{i} > amp_thresholds(i)) - 1;
    waveforms{i} = CalcSTA(data_pp.data', sts, [-nlrpoints nlrpoints]);
end

% Compare updated waveforms to initial estimates
%*** Hide this stuff somewhere else!
if (params.general.plot_diagnostics)
  num_waveforms = length(waveforms);
  cols= hsv(num_waveforms);
  nchan=size(data_pp.data,1);
    figure(params.plotting.first_fig_num+4);
    nc = ceil(sqrt(num_waveforms));
    nr = ceil(num_waveforms / nc);
    chSpace = 13; %**magic number, also in VisualizeClustering
    spacer = ones(size(waveforms{1},1), 1) * ([1:nchan]-1)*chSpace;
    for i = 1:num_waveforms
        subplot(nr, nc, i);
        inith = plot(reshape(init_waveforms{i},[],nchan)+spacer, 'k');
        hold on
        finalh = plot(reshape(waveforms{i},[],nchan)+spacer, 'Color', cols(i,:));
        hold off
        set(gca,'Xlim',[1, size(spacer,1)]);
        err = norm(init_waveforms{i} - waveforms{i})/...
              norm(waveforms{i});
        title(sprintf('cell %d, change=%.0f%%', i, 100*err))
        legend([inith(1) finalh(1)], {'Initial', 'New'});
    end
end

% Diagnostics for CBP waveforms:
% If recovered waveforms differ significantly from initial waveforms, then algorithm
% has not yet converged.  Execute this and go back to re-run CBP:
%     init_waveforms = waveforms;


%% ----------------------------------------------------------------------------------
% Post-analysis: Comparison of CBP to clustering results, and to ground truth (if
% available)

%** indicate which cells match ground truth.

ground_truth = load(data_pp.filename, 'true_spike_times', 'true_spike_class', 'dt');
ground_truth.filename = data_pp.filename;

if isfield(ground_truth, 'true_spike_times') && isfield(ground_truth, 'true_spike_class')
    % Reformat as 1-cellarr per cell of spike times.
    ground_truth.true_sp = GetSpikeTimesFromAssignments(ground_truth.true_spike_times, ground_truth.true_spike_class);
    
    % Reorder to match cell numbering from clustering.
    [ground_truth.true_sp, ...
     ground_truth.true_spike_times, ... 
     ground_truth.true_spike_class] = ReorderCells( ...
        ground_truth.true_sp, spike_times_cl, params.postproc.spike_location_slack);
end

if isfield(ground_truth, 'true_spike_times') && isfield(ground_truth, 'true_spike_class')
  % Since we already permuted ground truth to match clustering, this is true by definition
  best_ordering_cl = 1:length(spike_times_cl);
  best_ordering = 1:length(spike_times);

  % Evaluate clustering sorting
  [total_misses_cl, total_false_positives_cl, misses_cl, false_positives_cl] = ...
      evaluate_sorting(spike_times_cl, ground_truth.true_sp, params.postproc.spike_location_slack);
  fprintf('Clustering: %s', SortingEvaluationStr(ground_truth.true_sp, spike_times_cl, total_misses_cl, total_false_positives_cl));

  % Evaluate CBP sorting
  [total_misses, total_false_positives, prune_est_times, misses, false_positives] = ...
     EvaluateSorting(spike_times, spike_amps, ground_truth.true_sp, 'threshold', amp_thresholds, 'location_slack', params.postproc.spike_location_slack);
  fprintf('       CBP: %s', SortingEvaluationStr(ground_truth.true_sp, prune_est_times, total_misses, total_false_positives));

end

ground_truth = load(data_pp.filename, 'true_spike_times', 'true_spike_class', 'dt');
ground_truth.filename = data_pp.filename;

if isfield(ground_truth, 'true_spike_times') && isfield(ground_truth, 'true_spike_class')
    % Reformat as 1-cellarr per cell of spike times.
    ground_truth.true_sp = GetSpikeTimesFromAssignments(ground_truth.true_spike_times, ground_truth.true_spike_class);
    
    % Reorder to match cell numbering from clustering.
    [ground_truth.true_sp, ...
     ground_truth.true_spike_times, ... 
     ground_truth.true_spike_class] = ReorderCells( ...
        ground_truth.true_sp, spike_times_cl, params.postproc.spike_location_slack);
end

if isfield(ground_truth, 'true_spike_times') && isfield(ground_truth, 'true_spike_class')
  % Since we already permuted ground truth to match clustering, this is true by definition
  best_ordering_cl = 1:length(spike_times_cl);
  best_ordering = 1:length(spike_times);

  % Evaluate clustering sorting
  [total_misses_cl, total_false_positives_cl, misses_cl, false_positives_cl] = ...
      evaluate_sorting(spike_times_cl, ground_truth.true_sp, params.postproc.spike_location_slack);
  fprintf('Clustering: %s', SortingEvaluationStr(ground_truth.true_sp, spike_times_cl, total_misses_cl, total_false_positives_cl));

  % Evaluate CBP sorting
  [total_misses, total_false_positives, prune_est_times, misses, false_positives] = ...
     EvaluateSorting(spike_times, spike_amps, ground_truth.true_sp, 'threshold', amp_thresholds, 'location_slack', params.postproc.spike_location_slack);
  fprintf('       CBP: %s', SortingEvaluationStr(ground_truth.true_sp, prune_est_times, total_misses, total_false_positives));

end

%% ----------------------------------------------------------------------------------
% Plot various snippet subpopulations
[est_matches true_matches] = GreedyMatchTimes(spike_times, ground_truth.true_sp, params.postproc.spike_location_slack);

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
% desiredspiketimes = cell2mat(completemisses); % Pick one or the other...
desiredspiketimes = cell2mat(spike_times);

% Now find and plot the relevant snippets
snipindices = FindSnippets(desiredspiketimes, snippet_centers, snippets); % snippets just given as shorthand for calculating widths
ScrollSnippets(snippets, snippet_centers, ...
    'snipindices',  unique(snipindices(snipindices > 0)),  ...
    'cbp',          spike_times,        ...
    'cbpamp',       spike_amps,         ...
... %     'cbpampthresh', amp_thresholds,      ... % Could use amp_thresholds if we used that to pick snippets...
    'clust',        spike_times_cl,     ...
... %     'recons',       recon_snippets,     ...
    'true',         ground_truth.true_sp);


%% ----------------------------------------------------------------------------------
% Visualize true spike assignments in PC-space

if isfield(ground_truth, 'true_spike_class') && isfield(ground_truth, 'true_spike_times')
    cluster_pars = params.clustering;
    if isempty(cluster_pars.window_len), cluster_pars.window_len = params.general.waveform_len; end
    cluster_pars.align_mode = data_pp.polarity;

    Xstar = ConstructSnippetMatrix(data_pp.data, ground_truth.true_spike_times, cluster_pars);
    % Remove mean component and project onto leading PC's
%    XProjstar = (Xstar - repmat(mean(Xstar, 2), 1, size(Xstar, 2)))' * PCs;
XProjstar=Xstar'*PCs;
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

%% ----------------------------------------------------------------------------------
% Get greedy spike matches and plot RoC-style
% NB: Much faster if mex greedymatchtimes.c is compiled
%*** show chosen threshold in top plot
%*** also show log # spikes found?
PlotCBPROC(spike_times, spike_amps, ground_truth.true_sp, params.postproc.spike_location_slack);
