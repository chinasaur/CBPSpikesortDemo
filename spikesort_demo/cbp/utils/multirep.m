function x2 = multirep(x,lens)

x2 = zeros(sum(lens),1);
offset = 0;
for i=1:length(lens)
    x2(offset+1:offset+lens(i)) = x(i);
    offset = offset + lens(i);
end
