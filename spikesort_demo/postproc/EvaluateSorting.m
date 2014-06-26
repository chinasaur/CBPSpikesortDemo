function [total_misses, total_false_positives, ...
          est_times, misses, false_positives] = ...
    EvaluateSorting(est_times, est_amps, true_times, varargin)
% Given a set of estimated spike times and amplitudes and an amplitude
% threshold, together with a set of true spike times, determine the # of 
% misses and false positives
%
% est_times/amps : cell array of spike times/amplitudes for each waveform
%
% true_times : cell array of true spike times for each waveform
%              NOTE: this assumes that the est_times(1:m) are matched with
%              true_times(1:m) where 
%              m = min(length(est_times), length(true_times))
%
% opts : parameters with fields:
%        threshold : vector of thresholds for each waveform
%        location_slack : maximum number of samples between est/true spike
%                         to be considered a match
%
% Returns:
%   total_misses/false_positives : total numbers of misses/false positives
%                                  PER WAVEFORM
%   pruned_est_times : cell array of pruned est. spike times (by threshold)
%   misses : cell array of same size as true_times indicating misses
%   false_postives : cell array of same size as est_times (thresholded)

opts = inputParser();
opts.addParamValue('threshold', 0);
opts.addParamValue('location_slack', 0);
opts.parse(varargin{:});
opts = opts.Results;

% If no amplitudes are supplied assume they are all Inf (i.e. exceed any
% threhsold)
if (isempty(est_amps)) || all(cellfun(@isempty, est_amps))
    est_amps = cell(size(est_times));
    for i = 1:length(est_amps)
        est_amps{i} = Inf .* ones(size(est_times{i}));
    end
end

if (isscalar(opts.threshold))
    num_wfs = length(est_times);
    %fprintf('Assuming threshold=%0.3f for all %d waveforms.\n', ...
    %        pars.threshold, num_wfs);
    opts.threshold = opts.threshold .* ones(num_wfs, 1);
end

% Prune spikes according to thresholds
est_times = PruneSpikes(est_times, est_amps, opts.threshold);

[total_misses, total_false_positives, misses, ...
    false_positives] = ...
    evaluate_sorting(est_times, true_times, opts.location_slack);

% Include only those estimated spike times that exceed the threshold
function times = PruneSpikes(times, amps, threshold)
for i = 1 : length(times)
    times{i} = times{i}(amps{i} > threshold(i));
end