function s = cell_numel(S)

s = zeros(size(S));
1;
for i=1:numel(S)
    if (numel(S{i}) == 0)
        s(i) = 0;
    else
        s(i) = numel(S{i});
    end
end
    