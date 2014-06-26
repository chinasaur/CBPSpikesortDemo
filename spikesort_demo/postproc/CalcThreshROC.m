function [ampthreshes cumfps cummisses ...
    cumfprate cummissrate] = ...
    CalcThreshROC(true_times, est_times, est_amps, location_slack)

allsortamps = sort(est_amps(:));
ampthreshes = [0; allsortamps(1:end-1)];
cumfps    = zeros(size(ampthreshes));
cummisses = zeros(size(ampthreshes));
for i = 1:length(ampthreshes)
    pruned_est_times = est_times(est_amps > ampthreshes(i));
    [est_matches, true_matches] = ...
        greedymatchtimes(pruned_est_times, true_times, -location_slack, location_slack);
    cummisses(i) = sum(true_matches == 0);
    cumfps(i)    = sum(est_matches == 0);
end

cumfprate = cumfps ./ (length(est_times):-1:1)';
cummissrate = cummisses / length(true_times);