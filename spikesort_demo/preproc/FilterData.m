function datastruct = FilterData(datastruct, params)
datastruct = datastruct;
filt_pars = params.filtering;

% For output
log.operation = 'filtering';
log.params = filt_pars;

% filter frequenices specified?
if ~isempty(filt_pars.freq)
  % Pad with 0's
  paddata = PadTrace(datastruct.data, filt_pars.pad);
  
  % Highpass/Bandpass filtering
  data = FilterTrace(paddata, filt_pars, 1/(2*datastruct.dt));

  % Remove padding
  data = data(:, filt_pars.pad + 1 : end - filt_pars.pad);
else
   fprintf(1, 'No filtering performed (params.filtering.freq was empty).\n');
   data = datastruct.data;
end

% Remove CHANNEL-WISE means
data = data - repmat(mean(data, 2), 1, size(data, 2));

% Scale GLOBALLY across all channels
data = data ./ max(abs(data(:)));

% Output
datastruct.data = data;
datastruct.processing{end+1} = log;


% Plot filtered data, Fourier amplitude, and histogram of magnitudes
if (params.general.plot_diagnostics)
    % copied from WhitenNoise:
    dataMag = sqrt(sum(datastruct.data .^ 2, 1));
    nchan = datastruct.nchan;
    thresh = params.whitening.noise_threshold;
    minZoneLen = params.whitening.min_zone_len;
    if isempty(minZoneLen), minZoneLen = params.general.waveform_len/2; end
    noiseZones = GetNoiseZones(dataMag, thresh, minZoneLen);
    noiseZoneInds = cell2mat(cellfun(@(c) c', noiseZones, 'UniformOutput', false));
    zonesL = cellfun(@(c) c(1), noiseZones);  zonesR = cellfun(@(c) c(end), noiseZones);
                                                      
    figure(params.plotting.first_fig_num); subplot(3,1,2);
    noiseCol = [1 0.4 0.4];
    inds = params.plotting.dataPlotInds;
    visibleInds = find(((inds(1) < zonesL) & (zonesL < inds(end))) |...
                       ((inds(1) < zonesR) & (zonesR < inds(end))));
    nh=patch(datastruct.dt*[[1;1]*zonesL(visibleInds); [1;1]*zonesR(visibleInds)],...
             thresh*[-1;1;1;-1]*ones(1,length(visibleInds)), noiseCol,...
             'EdgeColor', noiseCol);
    hold on; 
    plot([inds(1), inds(end)]*datastruct.dt, [0 0], 'k'); 
    dh = plot((inds-1)*datastruct.dt, datastruct.data(:,inds)');
    hold off; 
    set(gca, 'Xlim', ([inds(1),inds(end)]-1)*datastruct.dt);
    set(gca,'Ylim',[-1 1]);  legend('noise regions, to be whitened');
    title('Filtered data, w/ regions to be used for noise covariance estimation');

    figure(params.plotting.first_fig_num+1); subplot(3,1,2);
    maxDFTind = floor(datastruct.nsamples/2);
    dftMag = abs(fft(datastruct.data,[],2));
    if (nchan > 1.5), dftMag = sqrt(sum(dftMag.^2)); end;
    plot(([1:maxDFTind]-1)/(maxDFTind*datastruct.dt*2), dftMag(1:maxDFTind));
    set(gca,'Yscale','log'); axis tight; 
    xlabel('frequency (Hz)'); ylabel('amplitude');
    title('Fourier amplitude of filtered data');
  
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
    xlabel('rms magnitude (over all channels)'); 
    legend([dh, ch], 'noise regions', 'chi2, fitted to noise regions');
    title('Histogram, cross-channel magnitude of filtered data');
end

%%-----------------------------------------------------------------
%% Auxilliary functions:

% Pad with constant values on left/right side to avoid border effects from
% filtering
function raw_data = PadTrace(raw_data, npad)
raw_data = padarray(raw_data, [0 npad], 0, 'both');
raw_data(:, 1 : npad) = repmat(raw_data(:, npad + 1), 1, npad);
raw_data(:, end - npad + 1 : end) = ...
    repmat(raw_data(:, end - npad), 1, npad);

% Filter trace
function data = FilterTrace(raw_data, filt_pars, rate)
if (isempty(filt_pars.freq))
  data = raw_data;
  return;
end
% Preprocessing parameters (see Harris et al 2000)
Wn = filt_pars.freq  / rate;

if ((length(Wn)==1) || (Wn(2) >= 1))
    switch(filt_pars.type)
        case 'butter'
            [fb, fa] = butter(filt_pars.order, Wn(1), 'high');
            fprintf('Highpass butter filtering with cutoff %f Hz\n', ...
                filt_pars.freq(1));
        case 'fir1'
            fb = fir1(filt_pars.order, Wn(1), 'high');
            fa = 1;
            fprintf('Highpass fir1 filtering with cutoff %f Hz\n', ...
                filt_pars.freq(1));
    end
else
    switch(filt_pars.type)
        case 'butter'
            [fb, fa] = butter(filt_pars.order, Wn, 'bandpass');
            fprintf('Bandpass butter filtering with cutoff %s Hz\n', ...
                mat2str(filt_pars.freq));
        case 'fir1'
            fb = fir1(filt_pars.order, Wn, 'bandpass');
            fa = 1;
            fprintf('Bandpass fir1 filtering with cutoff %s Hz\n', ...
                mat2str(filt_pars.freq));
    end
end

data = zeros(size(raw_data));
parfor chan = 1 : size(raw_data, 1)
    data(chan, :) = filter(fb, fa, raw_data(chan, :));

end
