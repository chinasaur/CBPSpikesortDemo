function grid_points = create_grid_points_uniform(data_dims, spacings)
% Create uniformly spaced grid-points and place them in a cell array,
% given a set of spacings
% Chaitanya Ekanadham 06/07/2011

% Arguments :
% data_dims : dimensions of a data example
% spacings : N x M matrix of spacings (N=#features,M=#transforms)

% Returns :
% grid_points : Nx1 cell array with each element a  k x M matrix of 
%               equally-spaced gridpoints for the Nth feature,
%               where k is determined by the data dims and Nth spacing.

[num_features num_transforms] = size(spacings);
data_dims = data_dims(data_dims ~= 1);
grid_points = cell(num_features, 1);
for feature_num = 1 : num_features
    if (num_transforms == 1) % 1-D just create equally spaced gridpoints
        grid_points{feature_num} = ...
            (-floor(data_dims(1) / 2) : ...
            spacings(feature_num) : ...
            floor(data_dims(1) / 2))';
    elseif (num_transforms == 2)
        [x y] = meshgrid((-floor(data_dims(1) / 2) : ...
            spacings(feature_num, 1) : ...
            floor(data_dims(1) / 2))', ...
            (-floor(data_dims(2) / 2) : ...
            spacings(feature_num, 2) : ...
            floor(data_dims(2) / 2))');
        grid_points{feature_num} =  [x(:) y(:)];
    else
        error('more than 2 data dimensions not supported');
    end
end
