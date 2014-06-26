function [m1 m2] = GreedyMatchTimesAmps(t1, a1, t2, d)
% Pretty tricky to get this right actually.  Best way seems to be to pick
% the spike within the matching window whose amplitude is closest to 1.

if ~iscell(t1)
    [m1 m2] = GreedyMatchTimes1(t1, a1, t2, d);
    return
end

m1 = cell(size(t1));
m2 = cell(size(t1));
for i = 1:length(t1)
    if i > length(t1) || i  > length(t2) || isempty(t1{i}) || isempty(t2{i})
        continue;
    end
    
    thisd = d(1);
    if length(d) > 1
        thisd = d(i);
    end
    
    [m1{i} m2{i}] = GreedyMatchTimes1(t1{i}, a1{i}, t2{i}, thisd);
end


function [m1 m2] = GreedyMatchTimes1(t1, a1, t2, d)

m1 = zeros(size(t1));
m2 = zeros(size(t2));
unused = true(size(t1));

for i = 1:length(t2)
    idx = abs(t1 - t2(i)) < d;

    % Don't use t2 times that have already been matched
    idx = idx & unused;
    
    if any(idx)
        % At least one match; now find the match with the largest
        % amplitude.  This is a weird thing to do, but is the correct thing
        % to do for calculating cumulative errors as amplitude threshold
        % increases, e.g. for ROC-style plots
        [~, iamp] = sort(a1(idx));
        idx = find(idx);
        
        match = idx(iamp(end));
        unused(match) = false;
        m1(match) = i;
        m2(i) = match;
    end
end