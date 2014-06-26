function [misses, false_positives] = ComputeErrors(est_times, true_times, slack)
% Given a set of true/est spike times, compute # errors

[est_matches true_matches] = GreedyMatchTimes(est_times, true_times, slack);

% Do a greedy matching by marching through true spikes.
misses          = sum(true_matches == 0);
false_positives = sum(est_matches  == 0);