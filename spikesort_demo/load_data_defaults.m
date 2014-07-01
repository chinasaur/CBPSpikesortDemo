function [datastruct, params] = load_data_defaults(datastruct, params)

datastruct.nchan = size(datastruct.data, 1);
datastruct.nsamples = size(datastruct.data, 2);
datastruct.processing = {};

if strncmpi(datastruct.filename, 'C_', 2)
    % Quiroga
    datastruct.polarity = 'max';
    params.filtering.freq = [200];  % Low freq cutoff in Hz    

elseif strncmpi(datastruct.filename, 'harris_', 7)
    % Harris
    'harris'
    datastruct.polarity = 'min';
    params.filtering.freq = [800 1e4];  % Low/high cutoff in Hz
end