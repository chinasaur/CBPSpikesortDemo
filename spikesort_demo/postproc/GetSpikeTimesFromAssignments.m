function spike_times = GetSpikeTimesFromAssignments(times, assignments)

% Given a set of times and assignments, create a cell array with times
% corresponding to each class.

% Create a cell array of true spike times
classes = 1:max(assignments);
numclasses = length(classes);
spike_times = cell(numclasses, 1);
for i = 1:numclasses
    spike_times{i} = reshape(times(assignments == classes(i)), [], 1);
end
