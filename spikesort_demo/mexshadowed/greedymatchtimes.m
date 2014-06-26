function [m1 m2] = greedymatchtimes(t1, t2, w1, w2)

warning('TRIALEVENTS:no_mex', 'Running slow uncompiled version; please compile mex/greedymatchtimes.c');

m1 = zeros(size(t1));
m2 = zeros(size(t2));
unused = true(size(t1));

for i = 1:length(t2)
    diff = t1 - t2(i);
    idx = (w1 < diff) & (diff < w2);

    % Don't use t2 times that have already been matched
    idx = idx & unused;
    
    if any(idx) % matched
        match = find(idx, 1);
        unused(match) = false;
        m1(match) = i;
        m2(i) = match;
    end
end