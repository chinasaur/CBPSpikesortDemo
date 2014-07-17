function [d, params] = load_raw_data(identifier, params)

% Load raw electrode data from a file, and adjust any parameters
% that need to be specialized for that data set.

% Data are assumed to be in a matlab (.mat) file, containing:
%   data: channel x time matrix of voltage traces
%   dt : the temporal sampling interval (in seconds)
% The file may optionally also contain:
%   polarity : 'max' or 'min', indicating whether primary voltage
%              deflections are positive or negative.
%   true_spike_times : vector of ground truth spike times, if known
%   true_spike_class : vector of numerical cell labels for each time

switch identifier
    case 'Quiroga1'
        % Simulated data: single electrode.
        % From: Quiroga et. al., Neural Computation, 16:1661-1687, 2004.
        fprintf(1,'Loading Quiroga dataset 1...\n');
        filename = '../example_data/C_Easy1_noise015.mat';
        d = load(filename, 'data', 'dt');
        d.polarity = 'max';
  
        params.filtering.freq = [];  % Low freq cutoff in Hz    

    case 'Harris1'
        % Real data: tetrode + one ground-truth intracellular electrode, rat hippocampus
        % From: Harris et. al., J. Neurophysiology, 84:401-414, 2000.
        fprintf(1,'Loading Harris dataset 1...\n');
        filename = '../example_data/harris_d533101.mat';
        d = load(filename, 'data', 'dt');
        d.polarity = 'min';
   
        params.filtering.freq = [800 1e4];  % Low/high cutoff in Hz
  
    otherwise
        error('Unrecognized dataset identifier', identifier);
end

d.filename = filename;
d.nchan = size(d.data, 1);
d.nsamples = size(d.data, 2);
d.processing = {};

fprintf (1, '  ... %d sec (%.1f min) of voltage data at %.1f kHz on %d channel(s).\n', ...
	 round(d.nsamples*d.dt), ...
	 (d.nsamples*d.dt/60), ...
	 1/(d.dt*1000), ...
	 d.nchan);

if (d.dt > 1/5000)
  warning('Sampling rate is %.1f kHz, but recommended minimum is 5kHz', 1/(1000*d.dt)); 
end

if (params.general.plot_diagnostics)
    plotDur = min(2400,d.nsamples);   %** magic number: should be a parameter
    plotT0 = round((d.nsamples-plotDur)/2);
    inds = plotT0+[1:plotDur];
    params.plotting.dataPlotInds = inds;

    figure(params.plotting.first_fig_num); clf; subplot(3,1,1); 
    plot([inds(1), inds(end)]*d.dt, [0 0], 'k'); 
    hold on; plot((inds-1)*d.dt, d.data(:,inds)'); hold off
    axis tight; xlabel('time (sec)');  ylabel('voltage'); 
    title(sprintf('Partial raw data, nChannels=%d, dt=%.2fmsec', d.nchan, 1000*d.dt));

    figure(params.plotting.first_fig_num+1); clf; subplot(3,1,1); 
    noiseCol=[1 0.3 0.3];
    dftMag = abs(fft(d.data,[],2));
    if (d.nchan > 1.5), dftMag = sqrt(mean(dftMag.^2)); end;
    maxDFTind = floor(d.nsamples/2);
    maxDFTval = 1.2*max(dftMag(2:maxDFTind));
    hold on;
    if (~isempty(params.filtering.freq))
      yr = [min(dftMag(2:maxDFTind)), maxDFTval];
      xr = [0 params.filtering.freq(1)]; 
      patch([xr xr(2) xr(1)], [yr(1) yr yr(2)], noiseCol);
      if (length(params.filtering.freq) >  1)
          f2 = params.filtering.freq(2);
          if (f2 < 1/(d.dt*2))
            xr = [f2 1/(d.dt*2)];
            patch([xr xr(2) xr(1)], [yr(1) yr yr(2)], noiseCol);
          end
      end
      legend('Frequencies to be filtered');
    end
    plot(([1:maxDFTind]-1)/(maxDFTind*d.dt*2), dftMag(1:maxDFTind));
    hold off;
    axis tight; set(gca,'Ylim', [0 maxDFTval]);  set(gca, 'Yscale','log');
    xlabel('frequency (Hz)'); ylabel('amplitude'); 
    title('Fourier amplitude, averaged over channels');
end
