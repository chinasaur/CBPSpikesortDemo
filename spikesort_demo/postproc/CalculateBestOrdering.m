function best_ordering = CalculateBestOrdering(est_times, true_times, slack)
% Wrapper around EvaluateSorting that tries all permuations of waveforms to
% find best match between computed waveforms and true units.

num_est_wfs = length(est_times);
num_true_wfs = length(true_times);

% Enumerate all possible matchings of est. waveforms to true waveforms
orderings = perms(1 : num_est_wfs);
if (num_true_wfs < num_est_wfs)
    sk = factorial(num_est_wfs - num_true_wfs);
    orderings = orderings(1 : sk : end, 1 : num_true_wfs);
end

% Evaluate misses/false positives for each ordering
tm = zeros(size(orderings, 1), 1);
tfp = zeros(size(tm));
parfor i = 1 : size(orderings, 1) % for each permutation
    [m, fp] = evaluate_sorting(est_times(orderings(i, :)), true_times, slack);
	% total misses/false positives across waveforms
    tm(i) = sum(m);
    tfp(i) = sum(fp);
end

% Pick the best ordering.
[blah, min_idx] = min(tm + tfp);
best_ordering = orderings(min_idx, :);