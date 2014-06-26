function d = cell_mse(x1,x2)

if (length(x1) ~= length(x2))
    error('cell arrays must be same length');
end
d = zeros(size(x1));
for i=1:length(x1)
    if (sum(size(x1{i}) ~= size(x2{i})) > 0)
        error('cell array elements must be the same size');
    end
    d(i) = mean((x2{i}(:)-x1{i}(:)).^2) ./ var(x1{i}(:)) * 100;
end