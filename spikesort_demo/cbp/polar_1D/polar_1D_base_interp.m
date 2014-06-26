function [CUV r theta F] = polar_1D_base_interp(w,spacing,dim,F)
% Base polar 1D interpolator
% Chaitanya Ekanadham 06/07/2011

% Arguments : 

% w : matrix feature (temporal)
% spacing : integer specifying amount of shift (assumed to be doubled)
% dim : desired dimensions the output
% F : (optional) already gives the 3 points to operate

% Returns : 
% CUV : size(w,1) x size(w,2) x 3 matrix containing the C,U,V interpolators in each slice
% r,theta : interpolator parameters


% FIXME: (PHLI)
% This section is done in polar_1D_get_radii_angles too; ideally call there
% instead of having duplicate code
if (~exist('F','var'))
    F = fbshift(w, spacing / 2);
end
theta = 2*acos( dot(F(:,3)-F(:,2),F(:,3)-F(:,1)) / norm(F(:,3)-F(:,1)) / norm(F(:,2)-F(:,3)) );
r = norm(F(:,2)-F(:,3)) / (2 * cos((pi-theta)/2) ); % use cosine rule


CUV0 = F/[1 1 1;r*[cos(theta) 1 cos(theta)];r*[-sin(theta) 0 sin(theta)]];
1;
if (exist('dim','var') && prod(dim) ~= size(CUV0,1))
    if (mod(size(CUV0,1)/size(w,2),2) == 0 || mod(dim(1),2) == 0)
        error('dimension of output and signal should be odd');
    end
    CUV0 = reshape(CUV0,size(w,1),size(w,2),3);
    CUV = zeros(dim(1),size(w,2),3);   
    for i=1:size(CUV,3)
        if (size(CUV0,1) > dim(1))
            idx = ceil(size(CUV0, 1) / 2) + ...
                  (-floor(dim(1) / 2) : floor(dim(1) / 2));
            CUV(:, :, i) = CUV0(idx, :, i);
            % error('dimension of output should always be greater than feature');
        else
            CUV(:,:,i) = padarray(CUV0(:,:,i),(dim(1)-size(CUV0,1))/2,'both');
        end
    end    
else
    1;
    CUV = reshape(CUV0,[size(w,1) size(w,2) 3]);
end