function [C D data_idx] = hist_coeffs(M,R1out)

C = cell(length(M{1}),1);
if (nargout > 1 && nargin > 1)
    D = cell(size(C));
end
if (nargout > 2)
    data_idx = cell(size(C));
end

for i=1:length(M)    
    offset = 0;
    for j=1:length(M{i})
        C{j} = [C{j}; M{i}{j}];
        if (nargout > 1 && nargin > 1)
            D{j} = [D{j}; R1out{i}(offset+1:offset+size(M{i}{j},1))];
            data_idx{j} = [data_idx{j};i.*ones(size(M{i}{j},1),1)];
            offset = offset + size(M{i}{j},1);
        end        
    end                    
end


