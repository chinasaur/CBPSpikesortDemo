function params = load_default_parameters()

general.waveform_len    = 81;
general.noise_threshold = 0.15;
general.plot_diagnostics = true;

filtering.freq  = [100];       % Low/high cutoff in Hz, empty default means no filtering
filtering.type  = 'fir1';      % "fir1", "butter"
filtering.pad   = 5e3;         % padding (in samples) to avoid border effects
filtering.order = 50;          % See Harris 2000

whitening.p_norm       = 2;   % Sliding-window Lp-norm is computed and thresholded
whitening.num_acf_lags = 120; % # lags for which to est. ACF
whitening.reg_const    = 0;   % Regularization const for making the ACF PSD
whitening.min_zone_len = [];  % Minimum duration of a noise zone
                              % used for covariance estimation. Empty => general.waveform_len/2

clustering.num_waveforms    = 3;    % Number of clusters
clustering.percent_variance = 90;   % criteria for choosing # of PCs to retain
clustering.upsample_fac     = 5;    % upsampling factor
clustering.smooth_len       = 5;    % smoothing factor (in upsampled space)
clustering.window_len       = [];   % Empty means just use general.waveform_len
clustering.peak_len         = [];   % Empty means use half the window_len
clustering.downsample_after_align = true; % downsample after aligning

plotting.font_size = 16;

params.general      = general;
params.filtering    = filtering;
params.whitening    = whitening;
params.clustering   = clustering;
params.plotting     = plotting;