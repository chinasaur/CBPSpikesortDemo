function str = SortingEvaluationStr(true_times, est_times, total_misses, total_false_positives)

% Get total number of true spikes                
total_true = 0;
for i = 1:length(true_times)
    total_true = total_true + length(true_times{i});
end

% Get total number of est. spikes
total_est = 0;
for i = 1 : length(est_times)
    total_est = total_est + length(est_times{i});
end

str = sprintf('Misses=%d/%d (%0.2f%%) FPs=%d/%d (%0.2f%%)\n', ...
        sum(total_misses), total_true, ...
        sum(total_misses) / total_true * 100, ...
        sum(total_false_positives), total_est, ...
        sum(total_false_positives) / total_est * 100);