function [datastruct, params] = load_raw_data(identifier, params)

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
        fprintf(1,'Loading Quiroga dataset 1...\n');
        filename = 'C_Easy1_noise015.mat';
        datastruct = load(filename, 'data', 'dt');
        datastruct.polarity = 'max';
  
        params.filtering.freq = [100];  % Low freq cutoff in Hz    

    case 'Harris1'
        fprintf(1,'Loading Harris dataset 1...\n');
        filename = 'harris_d533101_v2.mat';
        datastruct = load(filename, 'data', 'dt');
   
        datastruct.polarity = 'min';
   
        params.filtering.freq = [800 1e4];  % Low/high cutoff in Hz
  
    otherwise
        error('Unrecognized dataset identifier', identifier);
end

if (datastruct.dt > 1/5000)
    warning('Sampling rate is less than recommended minimum of 5kHz');
end

datastruct.filename = filename;
datastruct.nchan = size(datastruct.data, 1);
datastruct.nsamples = size(datastruct.data, 2);
datastruct.processing = {};
