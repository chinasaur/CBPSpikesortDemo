function s = cell_length(S)

s = zeros(size(S));
1;
for i=1:numel(S)
    if (numel(S{i}) == 0)
        s(i) = 0;
    else
        s(i) = size(S{i},1);
    end
end
    