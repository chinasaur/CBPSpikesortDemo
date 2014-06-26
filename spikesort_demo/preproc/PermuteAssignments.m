function permassign = PermuteAssignments(assign, perm)
% PERMUTEASSIGNMENTS
% usage: assign = PermuteAssignments(assign, perm)
%

permassign = zeros(size(assign));
for i = 1:length(perm)
    permassign(assign == perm(i)) = i;
end