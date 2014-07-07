function params = load_default_parameters()

% Default parameters, assembled into a hierarchical structure to clarify which are used where. 

general.waveform_len    = 81;    % number of samples in spike waveform
general.plot_diagnostics = true;

plotting.first_fig_num = 1;  % resetting this will make figure display start from a new base, allowing 
                             % comparison between different parameter conditions
plotting.font_size = 12;

filtering.freq  = [100];       % Low/high cutoff in Hz, empty default means no filtering
filtering.type  = 'fir1';      % "fir1", "butter"
filtering.pad   = 5e3;         % padding (in samples) to avoid border effects
filtering.order = 50;          % See Harris 2000

whitening.noise_threshold = 0.15;  % threshold, rel. to 1, used to detect and remove spikes, estimating
                                   % covariance of the remaining noise
whitening.p_norm          = 2;     % sliding-window Lp-norm is computed and thresholded
whitening.num_acf_lags    = 120;   % duration over which to estimate ACF
whitening.reg_const       = 0;     % regularization const for making the ACF PSD
whitening.min_zone_len    = [];    % minimum duration of a noise zone
				   % used for covariance estimation. Empty => use general.waveform_len/2

clustering.num_waveforms    = 3;    % Number of cells
clustering.threshold        = 8;    % threshold for picking spike-containing intervals (in stdevs))
clustering.window_len       = [];   % empty =>  use general.waveform_len
clustering.peak_len         = [];   % empty =>  use 0.5*window_len
clustering.percent_variance = 90;   % criteria for choosing # of PCs to retain
clustering.upsample_fac     = 5;    % upsampling factor
clustering.smooth_len       = 5;    % smoothing factor (in upsampled space)
clustering.downsample_after_align = true; % downsample after aligning

partition.threshold          = 3;  % threshold for silent 'break' regions (in stdevs)
partition.smooth_len         = 1;
partition.min_separation_len = floor(general.waveform_len/2);
partition.min_snippet_len    = general.waveform_len;
partition.max_snippet_len    = 1001;   % not enforced, only warnings
partition.min_pad_size       =  5;

params.general      = general;
params.plotting     = plotting;
params.filtering    = filtering;
params.whitening    = whitening;
params.clustering   = clustering;
params.partition    = partition;