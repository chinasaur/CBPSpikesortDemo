function v2 = zshift(v,k)

if (k > 0)
    v2 = [zeros(k,size(v,2));v(1:size(v,1)-k,:)];
else
    v2 = [v(-k+1:end,:);zeros(-k,size(v,2))];
end
    