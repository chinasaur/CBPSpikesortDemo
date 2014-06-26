function [spike_times, spike_mags] = ...
    ConvertSpikeTimesFromCell(spike_times_cell, spike_mags_cell, ...
                              snippet_centers, threshold)

% Convert a cell array of spikes for each snippet into a cell array of
% spike times.

% Arguments
% spike_times_cell: num_snippet x 1 cell array. Each element is a
% num_features x 1 cell array with spike times relative to a snippet center
%
% spike_mags_cell : cell array with same structure as spike_times_cell
% specifying spike magnitudes
%
% snippet_centers : num_snippets x 1 vector of snippet_centers
%
% threshold : (optional) threshold on spike magnitudes (default = 0)



if (length(spike_times_cell) ~= length(snippet_centers))
    error('cell array must have as many elements as snippet centers!');
end

if (~exist('threshold', 'var'))
    threshold = 0;
end

num_snippets = length(spike_times_cell);
num_features = length(spike_times_cell{1});
spike_counts = zeros(num_features, 1);
% First figure out how many spike of each feature
for snippet_num = 1 : num_snippets
    for feature_num = 1 : num_features
            spikes = spike_times_cell{snippet_num}{feature_num};
            spike_mags = spike_mags_cell{snippet_num}{feature_num};
            spikes = spikes(spike_mags > threshold);
            spike_counts(feature_num) = spike_counts(feature_num) + ...
                                        length(spikes);
    end
end

fprintf('Spike counts: %s\n', mat2str(spike_counts));

% Allocate the memory
spike_times = cell(num_features, 1);
spike_mags = cell(num_features, 1);
for feature_num = 1 : num_features
    spike_times{feature_num} = zeros(spike_counts(feature_num), 1);
    spike_mags{feature_num} = zeros(spike_counts(feature_num), 1);
end

% Now fill up the cell array
offsets = zeros(num_features, 1);
for snippet_num = 1 : num_snippets
    for feature_num = 1 : num_features
            times = snippet_centers(snippet_num) + ...
                    spike_times_cell{snippet_num}{feature_num};
            amplitudes = spike_mags_cell{snippet_num}{feature_num};
            times = times(amplitudes > threshold);
            if isempty(times), continue; end
                
            spike_times{feature_num}(offsets(feature_num) + 1 : ...
                offsets(feature_num) + length(times)) = ....
                    times;
            spike_mags{feature_num}(offsets(feature_num) + 1 : ...
                offsets(feature_num) + length(times)) = ....
                    amplitudes(amplitudes > threshold);                
            offsets(feature_num) = offsets(feature_num) + length(times);
    end
end

                            