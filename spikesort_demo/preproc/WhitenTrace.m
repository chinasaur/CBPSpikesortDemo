function [ECwhitened, ...
          old_acfs, whitened_acfs, ...
          NoiseCovMtx, WhitenedNoiseCovMtx] = ...
    WhitenTrace(EC, noise_comp_idx, num_acf_lags, reg_const, plot_flag)

%function ECwhitened = ...
%    WhitenTrace(EC, noise_comp_idx, num_acf_lags, reg_const, plot_flag)
%
% Utility function for whiteining voltage traces, given estimated "noise
% regions" within the trace.
% 
% EC : time x chan voltage trace matrix
%
% noise_comp_idx : cell array of time indices corresponding to "noise
%                  regions"
%
% num_acf_lags : number of lags at which to estimate temporal
%                autocorrelation
%
% reg_const : regularization constant to force ACF matrix to be PSD
%
% plot_flag : for "before-and-after" visualization of whitening process
%
% Returns the whitened trace (matrix same size as EC) ad well as the
% temporal and spatial covariance functions pre- and post-whitening

fprintf('Whitening trace...\n');
ECwhitened = zeros(size(EC));
num_channels = size(EC, 2);    
        
if (~exist('reg_const', 'var'))
	reg_const = 0;
end
nr = ceil(sqrt(num_channels));
old_acfs = cell(num_channels, 1);
whitened_acfs = cell(size(old_acfs));

for channel_num = 1 : num_channels
    
    % Construct noise samples for this channel.    
	noise = cell(size(noise_comp_idx));
    for comp_num = 1 : length(noise)
        noise{comp_num} = EC(noise_comp_idx{comp_num}, channel_num);
    end
           
    % Estimate the noise ACF for this channel
    noise_comps = GetSnippetsFromCellRanges(noise_comp_idx, ...
                                            EC(:, channel_num));
    noise_acf = EstimateACFFromCells(noise_comps, num_acf_lags);
    old_acfs{channel_num} = noise_acf;
    
    % Whiten each channel temporally    
	ECwhitened(:, channel_num) = ...
        WhitenTraceInTime(EC(:, channel_num), noise_acf, reg_const);

    whitened_noise_comps = ...
        GetSnippetsFromCellRanges(noise_comp_idx, ...
                                  ECwhitened(:, channel_num));
    whitened_noise_acf = ...
        EstimateACFFromCells(whitened_noise_comps, num_acf_lags);
    whitened_acfs{channel_num} = whitened_noise_acf;
       
end 
   
% Now whiten in space
NoiseCovMtx = EstimateNoiseCovMtx(noise_comp_idx, ECwhitened);
ECwhitened = WhitenTraceInSpace(ECwhitened, NoiseCovMtx, reg_const);
WhitenedNoiseCovMtx = EstimateNoiseCovMtx(noise_comp_idx, ECwhitened);
ECwhitened = ECwhitened';
end

function [acf acfs] = EstimateACFFromCells(noise_comps, num_acf_lags)
fprintf('Averaging ACF over %d noise regions.\n', length(noise_comps));    
acfs = zeros(num_acf_lags + 1, length(noise_comps));
for comp_num = 1 : length(noise_comps)
    acfs(:, comp_num) = ...
        EstimateACFFromSamples(noise_comps{comp_num}, num_acf_lags);
end
acf = sum(acfs, 2);
acf = acf - acf(end);
acf = acf ./ max(abs(acf));
end

function snippets = GetSnippetsFromCellRanges(idx_list, signal)

snippets = cell(length(idx_list), 1);
for comp_num = 1 : length(idx_list)
    snippets{comp_num} = signal(idx_list{comp_num}, :);
end
end

function NoiseCovMtx = EstimateNoiseCovMtx(noise_comp_idx, EC)
    noise_comps = GetSnippetsFromCellRanges(noise_comp_idx, EC);
    noise = cell2mat(noise_comps(:));
    noise = noise - repmat(mean(noise, 1), size(noise, 1), 1);
    NoiseCovMtx = noise' * noise ./ size(noise, 1);
end

function [ECwhitened, whitening_filter] = ...
    WhitenTraceInTime(EC, noise_acf, reg_const)

num_acf_lags = length(noise_acf) - 1;
% Estimate the noise autocorrelation function.
if (mod(num_acf_lags, 2) ~= 0)
    error('num_acf_lags must be even!');
end

% Compute a whitening filter by inverting the autocovariance matrix
T = toeplitz(noise_acf);
% Make sure this is positive definite!!!
T = sqrtm(T' * T);
if (norm(T - toeplitz(noise_acf)) > 1e-5)
    fprintf('Warning: toeplitz(noise_acf) was not positive definite. Forcing it.');
end
NoiseCovMtx = sqrtm(inv(T + reg_const * eye(length(noise_acf))));

whitening_filter = NoiseCovMtx(:, ceil(size(NoiseCovMtx, 2) / 2));
% Whiten the trace
ECwhitened = conv(EC, whitening_filter);
ECwhitened = ECwhitened(ceil(num_acf_lags / 2) : ...
             ceil(num_acf_lags / 2) + size(EC, 1) - 1, :);       
end

function ECwhitened = WhitenTraceInSpace(EC, NoiseCovMtx, reg_const)              
                   
% Estimate noise covariance matrix from noise samples.
WhiteningMtx = sqrtm(inv(NoiseCovMtx + ...
                         reg_const * eye(size(NoiseCovMtx, 1))));

ECwhitened = EC * WhiteningMtx';
end


                   