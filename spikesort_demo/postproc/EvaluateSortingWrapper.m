function [total_misses, total_false_positives, ...
          prune_est_times, misses, false_positives, ...
          best_ordering] = ...
          EvaluateSortingWrapper(est_times, est_amps, true_times, pars)
% Wrapper around EvaluateSorting that tries all permuations of waveforms.
% Returns same things as EvaluateSorting for the best permuation.

best_ordering = CalculateBestOrdering(est_times, true_times, pars.location_slack);

if (isempty(est_amps))
    est_amps = cell(size(est_times));
end

% Re-evaluate with this ordering
pars2 = pars;
pars2.threshold = pars.threshold(best_ordering);
[total_misses, total_false_positives, ...
 prune_est_times, misses, false_positives] = ...
	EvaluateSorting(est_times(best_ordering), ...
                    est_amps(best_ordering), ...
                    true_times, pars2);
% NOTE: total_misses (total_false_positives) are VECTORS specifying missed
% (false positive) spikes for each true (estimated) cell.

fprintf('Best ordering: %s ', mat2str(best_ordering));
fprintf(SortingEvaluationStr(true_times, est_times, prune_est_times, total_misses, total_false_positives));