function datastruct = FilterData(datastruct, filt_pars)
datastruct = datastruct;

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
   warning('FilterData:nofiltering', 'FILT_PARS.freq was empty; no filtering performed.');
   data = datastruct.data;
end

% Remove CHANNEL-WISE means
data = data - repmat(mean(data, 2), 1, size(data, 2));

% Scale GLOBALLY across all channels
data = data ./ max(abs(data(:)));

% Output
datastruct.data = data;
datastruct.processing{end+1} = log;



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
