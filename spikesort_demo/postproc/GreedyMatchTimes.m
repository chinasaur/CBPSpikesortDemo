function [m1 m2] = GreedyMatchTimes(t1, t2, d)

if ~iscell(t1)
    [m1 m2] = greedymatchtimes(t1, t2, -d, d);
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
    
    [m1{i} m2{i}] = greedymatchtimes(t1{i}, t2{i}, -thisd, thisd);
end
