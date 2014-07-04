function params = load_default_parameters()

general.num_waveforms    = 3;    % Number of cells
general.waveform_len    = 81;    % number of samples in spike waveform
general.plot_diagnostics = true;

filtering.freq  = [100];       % Low/high cutoff in Hz, empty default means no filtering
filtering.type  = 'fir1';      % "fir1", "butter"
filtering.pad   = 5e3;         % padding (in samples) to avoid border effects
filtering.order = 50;          % See Harris 2000

whitening.noise_threshold = 0.15;  % threshold, rel. to 1, used to detect and remove spikes, estimating
                                   % covariance of the remaining noise
whitening.p_norm       = 2;   % Sliding-window Lp-norm is computed and thresholded
whitening.num_acf_lags = 120; % # lags for which to est. ACF
whitening.reg_const    = 0;   % Regularization const for making the ACF PSD
whitening.min_zone_len = [];  % Minimum duration of a noise zone
                              % used for covariance estimation. Empty => general.waveform_len/2

clustering.percent_variance = 90;   % criteria for choosing # of PCs to retain
clustering.upsample_fac     = 5;    % upsampling factor
clustering.smooth_len       = 5;    % smoothing factor (in upsampled space)
clustering.window_len       = [];   % Empty means just use general.waveform_len
clustering.peak_len         = [];   % Empty means use half the window_len
clustering.downsample_after_align = true; % downsample after aligning

plotting.first_fig_num = 1;
plotting.font_size = 12;

params.general      = general;
params.filtering    = filtering;
params.whitening    = whitening;
params.clustering   = clustering;
params.plotting     = plotting;