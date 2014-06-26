function xvec = vectorize_cell(X)

xvec = zeros(sum(cell_numel(X)),1);
offset = 0;
for i=1:length(X)
    xvec(offset+1:offset+numel(X{i})) = X{i}(:);
    offset = offset + numel(X{i});
end