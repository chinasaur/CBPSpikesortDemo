function [Locations Magnitudes reconstructed_data] = ...
    polar_1D_extract_fn_cvx(coeffs, num_blocks, ...
                            Dictionary, radii_per_block, ...
                            grid_points, ...
                            spacings, thetas, magnitude_threshold)
% Given a vector of raw coefficients from CVX, 
% extract the locations and magnitudes.

% Arguments:
% coeffs : raw coefficient vector (output of CVX / ECOS)
% num_blocks : vector specifying number of blocks per feature
% Dictionary : dictionary matrix.
% radii_per_block : vector of size sum(num_blocks) x 1 specifying the polar
%                   radius associated with each block.
% grid_points : cell array of grid points on which dictionary lives.
% spacings : spacings for each feature's grid points
% thetas : angles for each feature's polar construction

% Returns:
% Locations/Magnitudes : cell array of each feature's locations/magnitudes
% coeffs : portion of raw_coeffs vector corresponding to the feature coeffs
% reconstructed_data : dictionary's approximation of the data.

num_features = length(num_blocks);
Locations = cell(num_features, 1);
Magnitudes = cell(num_features, 1);

offset = 0;
block_size = 3;  % Size of the interpolator group.
for feature_num = 1 : num_features    
    idx = (offset + 1) : (offset + block_size * num_blocks(feature_num));
    cx = coeffs(idx(1 : block_size : end));
    ux = coeffs(idx(2 : block_size : end));
    vx = coeffs(idx(3 : block_size :  end));    
    if (num_blocks(feature_num) > 0)
        Locations{feature_num} = grid_points{feature_num} + ...
                                 atan( vx ./ ux ) ./ ...
                                 thetas(feature_num) * ...
                                 spacings(feature_num) / 2;                
        delete_idx = cx < magnitude_threshold;
        cx(delete_idx) = 0;
        ux(delete_idx) = 0;
        vx(delete_idx) = 0;        
        Magnitudes{feature_num} = cx;                        
    end       
    coeffs(idx) = reshape([cx ux vx]', [], 1)';
    offset = offset + block_size * num_blocks(feature_num);
end

if (nargout > 2)
    reconstructed_data = Dictionary * coeffs;
end


