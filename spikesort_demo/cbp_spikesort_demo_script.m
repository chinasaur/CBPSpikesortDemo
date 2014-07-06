%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test script for CBP spike sorting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please see the file README.md for a description of this code.

%% -----------------------------------------------------------------
% Setup

% Set the current working directory to the directory containing this file:
% cd spikesort_demo

% Run the setup function, which sets paths and prints warnings or errors if there
% are issues detected (for example, mex/C files that need to be
% compiled for your system).
spikesort_demo_setup(pwd());

%% -----------------------------------------------------------------
% Step 0: Load raw electrode data

% Load default parameter set
params = load_default_parameters();

% Load an example data set, which includes the raw data, the timestep, and 
% (optionally) ground truth spike times.
switch 1 
    case 1
       % Simulated data: single electrode.
       % From: Quiroga et. al., Neural Computation, 16:1661-1687, 2004.
       [data, params] = load_raw_data('Quiroga1', params);
    case 2
       % Real data: tetrode + one ground-truth intracellular
       % electrode, rat hippocampus
       % From: Harris et. al., J. Neurophysiology, 84:401-414, 2000.
       [data, params] = load_raw_data('Harris1', params);
end

% Check that sampling rate is reasonable:
if (data.dt > 1/5000)
    warning('Sampling rate is less than recommended minimum of 5kHz');
end

% Plot a bit of raw data, to make sure it looks as expected:
if (params.general.plot_diagnostics)
    plotDur = min(2400,size(data.data,2));   %** magic number: should be a parameter
    plotT0 = round((size(data.data,2)-plotDur)/2);
    inds = plotT0+[1:plotDur];

    figure(params.plotting.first_fig_num); clf; subplot(3,1,1); 
    plot([inds(1), inds(end)]*data.dt, [0 0], 'k'); 
    hold on; plot((inds-1)*data.dt, data.data(:,inds)'); hold off
    axis tight; xlabel('time (sec)');  ylabel('voltage'); 
    title(sprintf('Partial raw data, nChannels=%d, dt=%.2fmsec', size(data.data,1), 1000*data.dt));

    figure(params.plotting.first_fig_num+1); clf; subplot(3,1,1); 
    noiseCol=[1 0.3 0.3];
    dftMag = abs(fft(data.data,[],2));
    if (size(dftMag,1) > 1.5), dftMag = sqrt(mean(dftMag.^2)); end;
    maxDFTind = floor(size(data.data,2)/2);
    maxDFTval = 1.2*max(dftMag(2:maxDFTind));
    hold on;
    if (~isempty(params.filtering.freq))
      yr = [min(dftMag(2:maxDFTind)), maxDFTval];
      xr = [0 params.filtering.freq(1)]; 
      patch([xr xr(2) xr(1)], [yr(1) yr yr(2)], noiseCol);
      if (length(params.filtering.freq) >  1)
          f2 = params.filtering.freq(2);
          if (f2 < 1/(data.dt*2))
            xr = [f2 1/(data.dt*2)];
            patch([xr xr(2) xr(1)], [yr(1) yr yr(2)], noiseCol);
          end
      end
      legend('Frequencies to be filtered');
    end
    plot(([1:maxDFTind]-1)/(maxDFTind*data.dt*2), dftMag(1:maxDFTind));
    hold off;
    axis tight; set(gca,'Ylim', [0 maxDFTval]);  set(gca, 'Yscale','log');
    xlabel('frequency (Hz)'); ylabel('amplitude'); 
    title('Fourier amplitude, averaged over channels');
end

%% -----------------------------------------------------------------
% Preprocessing Step 1: Temporal filtering

% Remove low and high frequencies - purpose is to eliminate non-signal
% parts of the frequency spectrum, and enable crude identification and
% removal of segments containg spikes via local amplitude
% thresholding, after which the background noise covariance can be
% estimated.
%   - pad with constant values (to avoid border effects when filtering)
%   - filter with band/highpass filter
%   - trim padding
%   - remove mean from each channel
%   - scale globally (across all channels) to be in range[-1, 1]
%
% The pre-filtering depends on values in params.filtering, which need
% to be adjusted for each data set:
%   - freq : range of frequencies (in Hz) for designing filter
%            Set to [] to turn off pre-filtering.
%   - type : type of filter for preprocessing. Currently supports
%            "fir1" and "butter"
%   - pad  : number of constant-value samples to pad
%   - order : order of the filter

filtdata = FilterData(data, params.filtering);

%Plot filtered data, Fourier amplitude, and histogram of magnitudes
if (params.general.plot_diagnostics)
    % copied from PreprocessTrace:
    dataMag = sqrt(sum(filtdata.data .^ 2, 1));
    nchan = size(filtdata.data,1);
    thresh = params.whitening.noise_threshold;
    minZoneLen = params.whitening.min_zone_len;
    if isempty(minZoneLen), minZoneLen = params.general.waveform_len/2; end
    noiseZones = GetNoiseZones(dataMag, thresh, minZoneLen);
    noiseZoneInds = cell2mat(cellfun(@(c) c', noiseZones, 'UniformOutput', false));
    zonesL = cellfun(@(c) c(1), noiseZones);  zonesR = cellfun(@(c) c(end), noiseZones);
                                                      
    % Plot filtered data, indicating noise regions that will be used in next
    % step for estimation of noise covariance:
    figure(params.plotting.first_fig_num); subplot(3,1,2);
    noiseCol = [1 0.4 0.4];
    visibleInds = find(((inds(1) < zonesL) & (zonesL < inds(end))) |...
                       ((inds(1) < zonesR) & (zonesR < inds(end))));
    nh=patch(filtdata.dt*[[1;1]*zonesL(visibleInds); [1;1]*zonesR(visibleInds)],...
             thresh*[-1;1;1;-1]*ones(1,length(visibleInds)), noiseCol,...
             'EdgeColor', noiseCol);
    hold on; 
    plot([inds(1), inds(end)]*filtdata.dt, [0 0], 'k'); 
    dh = plot((inds-1)*filtdata.dt, filtdata.data(:,inds)');
    hold off; 
    set(gca, 'Xlim', ([inds(1),inds(end)]-1)*filtdata.dt);
    set(gca,'Ylim',[-1 1]);  legend('noise regions, to be whitened');
    title('Filtered data, w/ regions to be used for noise covariance estimation');

    figure(params.plotting.first_fig_num+1); subplot(3,1,2);
    dftMag = abs(fft(filtdata.data,[],2));
    if (nchan > 1.5), dftMag = sqrt(sum(dftMag.^2)); end;
    plot(([1:maxDFTind]-1)/(maxDFTind*filtdata.dt*2), dftMag(1:maxDFTind));
    set(gca,'Yscale','log'); axis tight; 
    xlabel('frequency (Hz)'); ylabel('amplitude');
    title('Fourier amplitude of filtered data');
    
    % Plot histogram of magnitudes, indicating threshold, and chi distribution
    figure(params.plotting.first_fig_num+2); clf
    [N,X] = hist(dataMag, 100); 
    [Nnoise] = hist(dataMag(noiseZoneInds), X);
    sd = sqrt(sum(cellfun(@(c) sum(dataMag(c).^2), noiseZones)) / ...
              (nchan*sum(cellfun(@(c) length(c), noiseZones))));
    chi = 2*(X/sd).*chi2pdf((X/sd).^2, nchan);
    bar(X,N); set(gca,'Yscale','log'); yrg= get(gca, 'Ylim');
    hold on; 
    dh= bar(X,Nnoise); set(dh, 'FaceColor', noiseCol);
    ch= plot(X, (max(N)/max(chi))*chi, 'g'); 
    hold off; set(gca, 'Ylim', yrg);
    xlabel('voltage rms magnitude (over all channels)'); 
    legend([dh, ch], 'noise', 'fitted chi2');
    title('Histogram of magnitudes');
end

% Diagnostics for filtering: 
% In the next step, noise covariance will be
% estimated from regions below threshold (fig1, red).  There should be no
% spikes in these regions.  As an additional check, the
% below-threshold portion of the histogram (fig 3, red) should look
% like a chi-square distribution (green).  If spikes appear to be
% included in the noise segments, reduce
% params.whitening.noise_threshold, or modify the filtering parameters
% in params.filtering, and re-run the filtering step.

%% -----------------------------------------------------------------
% Preprocessing Step 2: Estimate noise covariance and whiten

% Estimate and whiten the noise, assuming space/time separability. This makes the L2-norm 
% portion of the objective into a simple sum of squares, greatly improving 
% computational efficiency. 
%   - estimate "noise zones" (regions of length min_zone_len or more 
%     samples whose Lp-norm are under noise_threshold).
%   - compute noise ACF from these zones
%   - whiten each channel in time by taking inverse matrix sqrt of ACF
%   - whiten across channels by left-multiplying each time slice 
%     by inverse matrix sqrt of the covariance matrix.
%
%  The whitening process depends on elements of params.whitening::
%   - whitening.noise_threshold : this threshold is used to crudely separate noise 
%     from signal, for purposes of estimating the noise covariance. The 
%     root-mean-squared signal value (across channels) is compared to this 
%     threshold.  
%   - reg_const : regularization constant.  This multiple of the identity 
%     matrix is added to the measured covariance matrix.

data_pp = PreprocessTrace(filtdata, params);

%Plot whitened data, Fourier amplitude, and histogram of magnitudes
if (params.general.plot_diagnostics)
    %taken from EstimateInitialWaveforms:
    dataMag = sqrt(sum(data_pp.data .^ 2, 1));
    nchan = size(filtdata.data,1);
    chiMean = sqrt(2)*gamma((nchan+1)/2)/gamma(nchan/2);
    chiVR = nchan - chiMean^2;
    thresh = chiMean + 6*sqrt(chiVR) %**Magic number - yuk!

    peakInds = dataMag(inds)> thresh;  
    peakLen = params.clustering.peak_len;
    if (isempty(peakLen)), peakLen=floor(params.general.waveform_len/2); end;
    for i = -peakLen : peakLen
        peakInds = peakInds & dataMag(inds) >= dataMag(inds+i);
    end
    peakInds = inds(peakInds);

    figure(params.plotting.first_fig_num); subplot(3,1,3); cla
    sigCol = [0.4 1 0.5];
    hold on;
    sh = patch(data.dt*[[1;1]*(peakInds-peakLen); [1;1]*(peakInds+peakLen)], ...
          thresh*[-1;1;1;-1]*ones(1,length(peakInds)), sigCol,'EdgeColor',sigCol);
    plot((inds-1)*data.dt, data_pp.data(:,inds)');
    hold off
    axis tight
    title('Filtered & noise-whitened data')
    legend('Segments exceeding threshold (putative spikes)');
    
    figure(params.plotting.first_fig_num+1); subplot(3,1,3);
    dftMag = abs(fft(data_pp.data,[],2));
    if (size(dftMag,1) > 1.5), dftMag = sqrt(mean(dftMag.^2)); end;
    plot(([1:maxDFTind]-1)/(maxDFTind*data.dt*2), dftMag(1:maxDFTind));
    set(gca, 'Yscale', 'log'); axis tight; 
    xlabel('frequency (Hz)'); ylabel('amplitude');
    title('Fourier amplitude of filtered & noise-whitened data');
    
    figure(params.plotting.first_fig_num+2); clf
    [N,X] = hist(dataMag, 100);
    Nspikes = hist(dataMag(dataMag>thresh), X);
    chi = 2*X.*chi2pdf(X.^2, nchan);
    bar(X,N); set(gca,'Yscale','log'); 
    yrg= get(gca, 'Ylim'); xrg= get(gca,'Xlim');
    hold on; 
    dh= bar(X,Nspikes); set(dh, 'FaceColor', sigCol);
    ch= plot(X, (max(N)/max(chi))*chi, 'r', 'LineWidth', 2);
    hold off; set(gca, 'Ylim', yrg); 
    title('Histogram of magnitudes of filtered & noise-whitened data');
    legend([dh, ch], 'putative spikes', 'fitted chi2');
end    

% Diagnostics for noise-whitening:
% Fig 4 shows original vs. whitened autocorrelation, which should be close
%   to a delta function (1 at 0, 0 elsewhere).  If not, try increasing
%   params.whitening.num_acf_lags.  If ACF is noisy, there may not be enough 
%   data samples for estimation.  This can be improved by 
%   a. increasing params.whitening.noise_threshold (allow more samples)
%   b. decreasing params.whitening.num_acf_lags 
%   c. decreasing params.whitening.min_zone_len (allow shorter noise zones)
% Fig 5 (for multi-electrodes only): shows correlation across electrodes,
%   which should be diagonal. If not, try the following:
%   a. increasing params.whitening.num_acf_lags
%   b. increasing params.whitening.min_zone_len
%   Note that this trades off with the quality the estimates (see prev).
% Fig 6 shows distributiobn of data in noise regions, which by default is 
%   assumed to be Gaussian.  If it doesn't match, may need to change the 
%   value of params.whitening.p_norm (set to 1 for Laplacian).
% Fig 1: Highlighted segments of whitened data (green) will be used to
%   estimate waveforms in the next step.  These should contain spikes, and
%   non-highlighted regions should contain background noise.  Threshold is
%   currently hard-coded to 3.9 (mean + 4*stdev of chi distribution).
% Fig 3: Histogram of magnitudes, and threshold used to identify
%   segments. 

%% -----------------------------------------------------------------
% Preprocessing step 3: Estimate initial spike waveforms

% Initialize spike waveforms:
%  - gather all data windows with L2-norm larger than threshold
%  - align spikes in these windows
%  - Compute PCA on these segments
%  - Perform K-means clustering in the subspace of the principal components
%    accounting for params.clustering.percent_variance portion of the total 
%    variance.
%
%  Parameters:
%  - num_waveforms : number of cells to be recovered from data.
%  - waveform_len : length (in samples) of spike waveforms.
%  - align_mode : During the clustering stage used to initialize spike
%                waveforms, segments of the trace are identified using a
%                threshold. These segments are then aligned based on this
%                parameter. Supported values are:
%                   maxrms : align wrt max L2 norm (RMS) value in segment
%                   maxabs : align wrt max L1 norm value in segment
%                   max : max signed sum across electrodes
%                   min : min signed sum across electrodes
%
%                For example, if the spike waveforms always have large
%                positive (negative) peaks, use max (min).

clear filtdata; % Recover some memory: don't need this any more

[centroids, assignments, X, XProj, PCs, snippet_centers_cl] = ...
    EstimateInitialWaveforms(data_pp, params);

% For later comparison to the CBP results, also grab the spike times
% corresponding to the segments assigned to each cluster:
spike_times_cl = GetSpikeTimesFromAssignments(snippet_centers_cl, assignments);

% Plot the data segments projected into the space of the 2 leading 
% principal components.
fignum = params.plotting.first_fig_num+3;
close(fignum+2);
VisualizeClustering(XProj, assignments, X, data_pp.nchan, fignum,fignum+1);

% Diagnostics for waveform initialization: 
% At this point, waveforms of all potential cells should
% be identified.  If not, may need to adjust
% params.general.num_waveforms and re-run the clustering to
% identify more/fewer cells.  May also wish to adjust the
% params.general.waveform_len, increasing it if the waveforms (Fig 5) are
% being chopped off, or shortening it if there is a substantial
% boundary region of silence.  If you do this, you should go back and
% re-run starting from the whitening step, since the waveform_len
% affects the identification of noise regions.

%% -----------------------------------------------------------------
% Preprocessing Step 4: partition data into snippets for improved efficiency of CBP

% Partition data into snippets, separated by "noise zones" in which
% the RMS of the waveforms does not surpass "threshold" for at least
% "min_separation_len" consecutive samples. Snippets are constrained 
% to have duration between min/max_snippet_len.
%
% NOTE: This is for computational speedup only, so choose a  
% conservative (low) threshold to avoid missing spikes!

%**TODO: move these into load_default_parameters
data_rms = sqrt(sum(data_pp.data .^ 2, 1));% root-mean-squared across channels
%cluster_threshold = 4 * median(data_rms) ./ 0.6745; % robust
%threshold = 4 * std(data_rms);
%*** Fix: make consistent with other choices  of threshold...
dof = size(data_pp.data,1);
chiMean = sqrt(2)*gamma((dof+1)/2)/gamma(dof/2);
chiVR = dof - chiMean^2;
threshold = chiMean + 3*sqrt(chiVR); 

partition_pars = struct('threshold', threshold, ...
                        'smooth_len', 1, ...
                        'min_separation_len', ...
                        floor(params.general.waveform_len/2), ...
                        'min_snippet_len', ...
                        params.general.waveform_len, ...
                        'max_snippet_len', ...
                        1001, ... % not enforced, only warnings
                        'min_pad_size', 5);
                    
[snippets, breaks, snippet_lens, snippet_centers, snippet_idx] = ...
    PartitionSignal(data_pp.data, partition_pars);
fprintf('Partitioned signal into %d chunks\n', length(snippets));

% Histogram of windowed norm of snippets versus "silences"
%**TODO: fix plotting for multi-electrodes
%**TODO: why windowed norm here?  The silences are detected using the norm
%smoothed by smooth_len....!
if (params.general.plot_diagnostics)
    figure(params.plotting.first_fig_num+2); clf; 
    len=params.general.waveform_len;
    wnAll = windowed_norm(data_rms, len);
    wnSnip = cellfun(@(s) windowed_norm(s', len), snippets, 'UniformOutput', false);
    wnBreak = windowed_norm(cell2mat(breaks')', len);
    [N,X] = hist(wnAll, 100);
    snipHist = hist(cell2mat(wnSnip')', X);
    breakHist = hist(wnBreak, X);
    bar(X,N);
    hold on;
    plot(X, 10*snipHist, 'g', 'LineWidth', 2);  %**TODO: why the 10* ??
    plot(X, breakHist, 'r', 'LineWidth', 2);
    hold off
    title('Histogram of windowed 2-norms');
    legend('all data', 'snippets', 'noise gaps');
    
    figure(params.plotting.first_fig_num+5); clf;
    snipLengths = cellfun(@(s) size(s,1), snippets);
    [N, X] = hist(snipLengths,100);
    bar(X,N); set(gca,'Yscale', 'log');
    title('Histogram of partitioned segment durations');
end

% Diagnostics for data partitioning: **TODO**
% Shouldn't see any large-amplitude content in the noise gaps.

%% -----------------------------------------------------------------
% CBP setup

% Close all preprocessing-related figures
figNums = params.plotting.first_fig_num+[0:5];
close(figNums(ishghandle(figNums)));

%**TODO: move default params into load_default_parameters
% Should be able to leave most of these defaults
num_waveforms = params.general.num_waveforms;

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
cbp_outer_pars.num_chan = data_pp.nchan;
for i = 1 : num_waveforms
    cbp_outer_pars.init_features{i} = ...
        reshape(centroids(:, i), [], cbp_outer_pars.num_chan);
end 

% cbp_pars are parameters for doing sparse inference.
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

% -----------------------------------------------------------------
% CBP parameters that should be adjusted by user

% For picking template delta
cbp_pars.accuracy = 0.1;

% Try single-spike soln first?
cbp_pars.compare_greedy = false; 
cbp_pars.greedy_p_value = 0; % tolerance to accept greedy soln
% cbp_pars.greedy_p_value = 1 - 1e-5;

% Corr. threshold below which atoms will not be used during CBP.
% For speedup only; set to 0 to disable
cbp_pars.prefilter_threshold = 0; %0.01, ... 

% -----------------------------------------------------------------
% Pick solver and reweighting parameters

% Setup CVX if needed
% addpath(fullfile(sstpath, '../cvx/'));
% cvx_setup;
% cvx_solver sedumi; 

% cbp_pars.solve_fn = @cbp_cvx;
% cbp_pars.solve_fn = @cbp_qcml_sumsquare;
% reweight_exp = 25 * [1 1 1];

cbp_pars.solve_fn = @cbp_ecos_2norm;
reweight_exp = 1.5 * ones(1, num_waveforms);

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

%% -----------------------------------------------------------------
% CBP step 1: use CBP to estimate spike times of all cells

starttime = tic;
cbp_pars.progress = true; % Set false if getting Java errors from progress bar
[spike_times, spike_amps, recon_snippets] = ...
    SpikesortCBP(snippets, ...
                 snippet_centers, ...
                 cbp_outer_pars, ...
                 cbp_pars);
toc(starttime);

%% -----------------------------------------------------------------
% Pick amplitude thresholds and interactively visualize effect on ACorr/XCorr

% NB: Much faster if mex trialevents.c is compiled

% Allow this much slack in time bins of spike locations, for live updating 
% ground truth feedback
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
%***Why is slack parameter needed here?
[atgf amp_threshold] = AmplitudeThresholdGUI(spike_amps, spike_times, 'dt', data_pp.dt, 'location_slack', spike_location_slack);

% -----------------------------------------------------------------
% Histogram of windowed norm for data, whitened data, residual

data_recon = cell(size(snippets));
for i = 1:numel(snippets)
    data_recon{i} = snippets{i} - recon_snippets{i};
end
wnsnip = cellfun(@(s) windowed_norm(s', params.general.waveform_len), snippets, 'UniformOutput', false);
wnbreak = windowed_norm(cell2mat(breaks')', params.general.waveform_len);
wnresid = cellfun(@(s) windowed_norm(s', params.general.waveform_len), data_recon, 'UniformOutput', false);

if (params.general.plot_diagnostics)
    figure(params.plotting.first_fig_num);
    [N,X] = hist(cell2mat(wnsnip')', 100);
    Nbreak = hist(wnbreak', X);
    Nresid = hist(cell2mat(wnresid')', X);
    chi = 2*X.*chi2pdf(X.^2, nchan*params.general.waveform_len);
    bar(X,N); set(gca,'Yscale','log'); yrg= get(gca, 'Ylim');
    hold on;
    plot(X,(max(N)/max(chi))*chi, 'c','LineWidth', 2);
    plot(X,(max(N)/max(Nbreak))*Nbreak,'r','LineWidth', 2);
    plot(X,(max(N)/max(Nresid))*Nresid, 'g', 'LineWidth', 2);
    hold off; set(gca,'Ylim', yrg);
    legend('all snippets', 'expected (Chi)', 'silences', 'snippet post-CBP residuals');
    title('Histogram of windowed 2-norms');
end

%**TODO: Plot [ do this for clustering also ]
% for each cell: waveform, total spikes and spike rate, amplitude distribution with thredhold indicated, autocorrelationPSTH 
% for each pair: cross-correlationPSTH, in time bins roughly the width of the main portion of the waveforms
%    ( for easy discernment of synchrony artifacts).
    
% Diagnostics for CBP spike-finding: ***TODO***
% amplitudes should be around 1.  autocorr should show zeros spikes in refractory period (roughly *** msec)
% crosscorr should not have a notch around zero (common in clustering methods)
% Residual should not have high-amplitude regions (if it does, increase cbp_pars.firing_rates, and re-run CBP).

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
        inith = plot(cbp_outer_pars.init_features{i}, 'b');
        hold on
        finalh = plot(waveforms{i}, 'r');

        err = norm(cbp_outer_pars.init_features{i} - waveforms{i})/...
              norm(waveforms{i});
        title(sprintf('Waveform %d, Rel error=%.2f', i, err))
        legend([inith(1) finalh(1)], {'Initial', 'New'});
    end
end

% Diagnostics for waveforms: 
%**TODO: Visualize spikes in PC space, compared to clustering result.  Allow user to increase/decrease 
% number of waveforms.

% if recovered waveforms differ significantly from initial waveforms, execute this and re-run CBP:
%      cbp_outer_pars.init_features = waveforms


%% -----------------------------------------------------------------
% Post-analysis: Compare to ground truth, if available

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
