function whitedatastruct = WhitenNoise(datastruct, params)
% Preprocess extracellularly recorded voltage trace by estimating noise
% and whitening if desired.
% FIXME: Appears to clip rawdata
% FIXME: In process cleaning up

%% NO longer used (previously used opts.plot instead of gen_pars.plot_diagnostics)
%opts = inputParser();
%opts.addParamValue('plot', true);
%opts.parse(varargin{:});
%opts = opts.Results;

gen_pars = params.general;
white_pars = params.whitening;
data = datastruct.data;

% For output
log.operation = 'whitening';
log.params.general = gen_pars;
log.params.whitening = white_pars;

min_zone_len = white_pars.min_zone_len;
if isempty(min_zone_len)
    min_zone_len = floor(gen_pars.waveform_len / 2);
end


% Standardization [EPS: moved to FilterData]

% Remove CHANNEL-WISE means
%data = data - repmat(mean(data, 2), 1, size(data, 2));

% Scale GLOBALLY across all channels
%data = data ./ max(abs(data(:)));



% Noise zone estimation and whitening

% Root-mean-squared of data
data_rms = sqrt(sum(data .^ 2, 1)); 

% Estimate noise zones
noise_zone_idx = GetNoiseZones(data_rms, ...
                               white_pars.noise_threshold, ...
                               min_zone_len);

if (0)
%**TODO: Doesn't make sense to show this, since it's the marginals
% BEFORE whitening.  All that matters is the marginals AFTER
% whitening (which the objective function for CBP relies on).
    PlotNoiseDbn(noise_zone_idx, data, params.plotting.first_fig_num+5, params.plotting.font_size);
% Fig 6: distribution of data in noise regions, which by default is assumed to
%   be Gaussian.  If it doesn't match, may need to change the value of
%   params.whitening.p_norm (set to 1 for Laplacian).
end

% Whiten trace if desired
[data_whitened, old_acfs, whitened_acfs, old_cov, whitened_cov] = ...
        WhitenTrace(data', ...
                    noise_zone_idx, ...
                    white_pars.num_acf_lags, ...
                    white_pars.reg_const, ...
                    gen_pars.plot_diagnostics);
noise_sigma = 1;
%fprintf('noise_sigma=%f\n', noise_sigma);


% Standardization AGAIN

if 0
% Remove CHANNEL-WISE means
data_whitened = data_whitened - ...
    repmat(mean(data_whitened, 2), 1, size(data_whitened, 2));
% Scale GLOBALLY across all channels
data_whitened = data_whitened ./ max(abs(data_whitened(:)));
end

% Output
whitedatastruct = datastruct;
whitedatastruct.data = data_whitened;
whitedatastruct.noise_sigma = noise_sigma;
whitedatastruct.processing{end+1} = log;


% Visualization of noise zone estimation
dt = datastruct.dt;
font_size = params.plotting.font_size;
nchan = size(data_whitened, 1);

if (0) %Old display code: unnecessary since this is plotted in FilterData.m
    figure; clf;
    plot_len = 5e3;
    middle_idx = ceil(size(data, 2) / 2) + (-plot_len : plot_len);
    xax = (0 : dt : (length(middle_idx) - 1) * dt) .* 1e3;
    h(1) = plot(xax, data_rms(middle_idx), 'k-'); hold on;
    for comp = 1 : length(noise_zone_idx)
        if (noise_zone_idx{comp}(1) <= middle_idx(end) && ...
            noise_zone_idx{comp}(end) >= middle_idx(1))
            sub_xax = (noise_zone_idx{comp} - middle_idx(1)) * dt * 1e3;
            h(2) = plot(sub_xax, data_rms(noise_zone_idx{comp}), 'r-');
        end
    end
    plot([xax(1) xax(end)], white_pars.noise_threshold * [1 1], 'Color', 0.8 * [1 1 1]);
    xlim([xax(1) xax(end)]);
    ylim([0 max(abs(data_rms))]);
    set(gca, 'FontSize', font_size);
    xlabel('time (ms)');
    ylabel('RMS of trace');
    title(sprintf('Estimation of noise zones (thres=%0.3f)', ...
          white_pars.noise_threshold));
    lg = legend(h, {'est. signal', 'est. noise'});
    set(lg, 'FontSize', font_size);
end


% Visualization of whitening effects
if (gen_pars.plot_diagnostics)
    figure(params.plotting.first_fig_num+3); clf
    nr = ceil(sqrt(nchan));    
    tax = (0 : dt : (length(old_acfs{1}) - 1) * dt)' .* 1e3;
    for chan = 1 : nchan
        subplot(nr, nr, chan);        
        plot(tax, [old_acfs{chan}, whitened_acfs{chan}], ...
              '.-', 'LineWidth', 1, 'MarkerSize', 14);        
        set(gca, 'FontSize', font_size);
        xlabel('time lag (ms)');
        ylabel('autocorrelation');
        legend('original', 'whitened');
        hold on; plot([tax(1), tax(end)], [0 0], 'k-');
        title(sprintf('Channel %d', chan));
    end
    if (nchan > 1.5)
      figure(params.plotting.first_fig_num+4); clf
      subplot(1, 2, 1), imagesc(old_cov); 
      colormap(gray); axis equal; axis tight;
      set(gca, 'FontSize', font_size); 
      title('Orig. cross-channel covariance');
      xlabel('channel');    ylabel('channel');
      set(gca, 'XTick', 1 : nchan, 'YTick', 1 : nchan);
      subplot(1, 2, 2), imagesc(whitened_cov); 
      colormap(gray); axis equal; axis tight;
      set(gca, 'FontSize', font_size);
      xlabel('channel');
      title('Whitened cross-channel covariance');    
      set(gca, 'XTick', 1 : nchan, 'YTick', 1 : nchan);
    end
end

% visualization of whitened data, FFT, and histograms
if (gen_pars.plot_diagnostics)
    inds = params.plotting.dataPlotInds;
    % copied from preproc/EstimateInitialWaveforms.m:
    dataMag = sqrt(sum(whitedatastruct.data .^ 2, 1));
    thresh = params.clustering.spike_threshold;
    peakInds = dataMag(inds)> thresh;  
    peakLen = params.clustering.peak_len;
    if (isempty(peakLen)), peakLen=floor(gen_pars.waveform_len/2); end;
    for i = -peakLen : peakLen
        peakInds = peakInds & dataMag(inds) >= dataMag(inds+i);
    end
    peakInds = inds(peakInds);

    figure(params.plotting.first_fig_num); subplot(3,1,3); cla
    sigCol = [0.4 1 0.5];
    hold on;
    sh = patch(whitedatastruct.dt*[[1;1]*(peakInds-peakLen); [1;1]*(peakInds+peakLen)], ...
          thresh*[-1;1;1;-1]*ones(1,length(peakInds)), sigCol,'EdgeColor',sigCol);
    plot((inds-1)*whitedatastruct.dt, whitedatastruct.data(:,inds)');
    hold off
    axis tight
    title('Filtered & noise-whitened data')
    legend('putative spike segments (for initial waveform estimation)');
    
    figure(params.plotting.first_fig_num+1); subplot(3,1,3);
    maxDFTind = floor(whitedatastruct.nsamples/2);
    dftMag = abs(fft(whitedatastruct.data,[],2));
    if (nchan > 1.5), dftMag = sqrt(mean(dftMag.^2)); end;
    plot(([1:maxDFTind]-1)/(maxDFTind*whitedatastruct.dt*2), dftMag(1:maxDFTind));
    set(gca, 'Yscale', 'log'); axis tight; 
    xlabel('frequency (Hz)'); ylabel('amplitude');
    title('Fourier amplitude, filtered & noise-whitened data');
    
    figure(params.plotting.first_fig_num+2); clf
    subplot(2,1,1);
    mx = max(abs(whitedatastruct.data(:)));
    X=linspace(-mx,mx,100);
    N=hist(whitedatastruct.data',X);
    plot(X,N); set(gca,'Yscale','log'); rg=get(gca,'Ylim');
    hold on;
    gh=plot(X, max(N(:))*exp(-(X.^2)/2), 'r', 'LineWidth', 2); 
    plot(X,N); set(gca,'Ylim',rg);
    hold off; 
    if (nchan < 1.5)
      title('Histogram, filtered/whitened data');
    else
      title(sprintf('Histograms, filtered/whitened data, %d channels', nchan));
    end
    legend(gh, 'univariate Gaussian');
    subplot(2,1,2);
    [N,X] = hist(dataMag, 100);
    Nspikes = hist(dataMag(dataMag>thresh), X);
    chi = 2*X.*chi2pdf(X.^2, nchan);
    bar(X,N); set(gca,'Yscale','log'); 
    yrg= get(gca, 'Ylim'); xrg= get(gca,'Xlim');
    hold on; 
    dh= bar(X,Nspikes); set(dh, 'FaceColor', sigCol);
    ch= plot(X, (max(N)/max(chi))*chi, 'r', 'LineWidth', 2);
    hold off; set(gca, 'Ylim', yrg); 
    title('Histogram, magnitude over filtered/whitened channel(s)');
    legend([dh, ch], 'putative spike segments', 'chi-distribution, univariate Gaussian');
end    


%% -----------------------------------------------------------------
%% Whitening routines

function PlotNoiseDbn(noise_zone_idx, data, fignum, font_size)

noise_samples = zeros(numel(data), 1);
offset = 0;
for i = 1 : length(noise_zone_idx)
    tmp = reshape(data(:, noise_zone_idx{i}), [], 1);
    
    noise_samples(offset + 1 : offset + length(tmp)) = tmp;
    offset = offset + length(tmp);
end
noise_samples = noise_samples(1 : offset);

%noise_sigma = median(abs(noise_samples)) / 0.6745; % [Harris]
noise_sigma = std(noise_samples);

% Compute gaussian and laplacian average likelihoods of noise samples
gauss_ll = -0.5 * sum(log(2 * pi * noise_sigma .^ 2)) - ...
           sum(mean(noise_samples .^ 2, 1)' ./ ...
           (2 * pi * noise_sigma .^ 2));
beta_hat = noise_sigma ./ sqrt(2);
laplace_ll = -sum(log(2 * beta_hat)) - ...
             sum(mean(abs(noise_samples), 1)' ./ beta_hat);

figure(fignum), clf;
min_val = quantile(noise_samples, 0);
max_val = quantile(noise_samples, 1);
x = linspace(min_val, max_val, 1e2)';
h = histc(noise_samples(:), x);
h = h./(sum(h))./ mean(diff(x));
plot(x, 1 / (2 * beta_hat) * exp(-abs(x) ./ beta_hat), 'g', ...
     'LineWidth', 2);
xlim([x(1) x(end)]);
hold on;
plot(x, normpdf(x, 0, noise_sigma), 'r', 'LineWidth', 2);
plot(x, h, '.'); set(gca, 'YScale', 'log');
set(gca, 'FontSize', font_size);
title(sprintf('Noise dbn (GaussLL=%0.3f LaplaceLL=%0.3f', ...
      gauss_ll, laplace_ll));
l = legend('laplacian', 'gaussian', 'data');
set(l, 'FontSize', font_size, 'Location', 'SouthWest');