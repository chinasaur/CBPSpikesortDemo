function W = vec2cell(w,lens,m)

if (nargin < 3)
    m =1 ;
end
if (sum(lens) * m ~= numel(w))
    error('Incompatible sizes for vector to cell conversion.');
end
offset = 0;
W = cell(length(lens),1);
for i=1:length(lens)
    W{i} = reshape(w(offset+1:offset+lens(i)*m),lens(i),m);
    offset = offset + lens(i)*m;
end