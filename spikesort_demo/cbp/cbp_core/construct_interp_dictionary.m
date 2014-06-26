function Dictionary = construct_interp_dictionary(BaseFeatures, ...
                                                  desired_dims, ...
                                                  interp_group_size, ...
                                                  base_interp_fn, ...
                                                  spacings, ...
                                                  transform_operator, ...
                                                  grid_points)
% Construct a dictionary matrix Dictionary.

% Arguments:
% BaseFeatures :  base features (N x 1 cell array of vectors/matrices)
% desired_dims : desired_dimension of the dictionary elements
% interp_group_size  : number of interpolators in each group
% base_interp_fn : interp-group constructor of form:
%               g = @(w,spacing)
%               that takes a feature and spacings and returns an 
%               interpolation group w.r.t. transformations transform_operators.
% spacings : spacings used to construct interpolators (N X M) matrix
% transform_operator : function handle of the form:
%                      x2 = @(x,  transform_params)
%                      That exactly transforms x.
% grid_points :  Nx1 cell array of matrices. 
%                Each element is a T_n x length(transform_operators) matrix
%                of grid points for that feature

% TODO: (PHLI) Would be nice to generate the dictionary in a simpler
% format, e.g. blockwise instead of interleaved

% TODO: (PHLI) Could use conv instead of zshift to do all the shifting in 
% one C call?

% Returns:
% Dictionary matrix.

num_features = length(BaseFeatures);

% Init a cell array that will hold interpolator groups for each feature.
BaseGroups = cell(num_features, 1); 
if (size(spacings, 1) == 1)
    spacings = repmat(spacings(:)', length(BaseFeatures), 1);
end

for feature_num =  1 : num_features
    BaseGroups{feature_num} = base_interp_fn(BaseFeatures{feature_num}, ...
                                             spacings(feature_num, :)', ...
                                             desired_dims);

    % The last dimension should index the group component.
    % The preceding dimensions should be desired_dims.
    base_group_size = size(BaseGroups{feature_num});        
    base_group_element_size = base_group_size(1 : end - 1);
    %base_group_element_dims = ...
    %    base_group_element_size(base_group_element_size ~= 1);
    if (norm(base_group_element_size - desired_dims) > 0 || ...
        base_group_size(end) ~= interp_group_size)
        error('base interpolator must have output size %s x %d', ...
              mat2str(desired_dims, 3), interp_group_size);
    end
end

% I believe for our setup features must always be the same length, but not
% willing to assume this.
maxfeaturelen = max(cellfun(@numel, BaseFeatures));

% How much memory do we need?
num_groups_per_feature = cellfun(@(C) size(C,1), grid_points);
n_dict_entries = sum(num_groups_per_feature) * interp_group_size;
spsize = maxfeaturelen * n_dict_entries;

Dictionary = spalloc(prod(desired_dims), n_dict_entries, spsize);

% Dictionary columns are ordered as follows:
%  ---- Feature 1 ---- | ---- Feature 2 ---- | ....
%[ group1, group2, ... | group1, group2, ... | ....]

% Each group is a submatrx with dimensions:
% prod(desired_dims) x interp_group_size
for col_num = 1 : size(Dictionary,2)
    % Get the group index
    group_num = ceil(col_num / interp_group_size);
    
    % Get the feature index.
    feature_num =  find(cumsum(num_groups_per_feature) >= ...
                        group_num, 1);            
    
    % Get index of interpolator component (within the group)
    interp_group_num = mod(col_num - 1, interp_group_size) + 1; 
    
    % Get the index of the grid point associated with this column in 
    % grid_points{feature_num}
    grid_idx = group_num - ...
               sum(num_groups_per_feature(1 : feature_num - 1));    
    grid_point = grid_points{feature_num}(grid_idx, :);
    
    % Transform basis function according to this gridpoint
    basis_fn = BaseGroups{feature_num}(:, :, interp_group_num);
    transformed_basis_fn = transform_operator(basis_fn, grid_point);

    Dictionary(:, col_num) = transformed_basis_fn(:);
end