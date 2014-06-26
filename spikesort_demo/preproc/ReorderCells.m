function [cells times classes] = ReorderCells(cells, ordercells, spiketimeslack)
% Reorder cells defined by cell arrays of spike times to match based on
% spike timings.
%
% FIXME: Fairly kludgy...
%

bestorder = CalculateBestOrdering(ordercells, cells, spiketimeslack);

% Invert permutation
true_perm = zeros(1, length(bestorder));
true_perm(bestorder) = 1:length(bestorder);

% Recover spike time and class vectors
times = cell2mat(cells);
[times timeidx] = sort(times);
classes = cell(size(cells));
for i = 1:length(classes)
    classes{i} = i * ones(size(cells{i}));
end
classes = cell2mat(classes);
classes = classes(timeidx);

% Permute
classes = PermuteAssignments(classes, true_perm);
cells = GetSpikeTimesFromAssignments(times, classes);