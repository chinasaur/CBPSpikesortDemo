function [prefilter_idx, grid_point_corrs] = ...
    PrefilterDictionary(data, Features, grid_points, threshold)

% Correlate a set of convolutional features (placed at some specified
% grid points) with a data example. Return the indices at which the
% correlation exceeds the desired threshold.

num_features = length(Features);
prefilter_cell_idx = cell(num_features, 1);
%if (size(Features{1}, 2) > 1)
%    error('Prefiltering multiple measurements unsupported!');
%end
grid_point_corrs = zeros(sum(cellfun(@(C) size(C,1), grid_points)), 1);
counter = 1;
for feature_num = 1 : num_features
    num_grid_points = size(grid_points{feature_num}, 1);
    prefilter_cell_idx{feature_num} = true(num_grid_points, 1);
    if (num_grid_points < 2)
        continue;
    end
        
    % Convolve the feature template with the data
    corrs = zeros(size(data));
    for channel_num = 1 : size(data, 2)
        corrs_i = conv(data(:, channel_num), ...
                     flipud(Features{feature_num}(:, channel_num)));
        corrs_i = corrs_i(ceil(size(Features{feature_num}, 1) / 2) : ...
                  (ceil(size(Features{feature_num}, 1) / 2) + ...
                   size(data, 1) - 1));
        corrs(:, channel_num) = corrs_i;
    end
    corrs = sum(corrs, 2) ./ norm(Features{feature_num}, 'fro') ./ ...
            norm(data, 'fro');
    corr_times = -floor(size(data, 1) / 2) : floor(size(data, 1) / 2);
    spacing = grid_points{feature_num}(2) - grid_points{feature_num}(1);
        
    for grid_point_num = 1 : num_grid_points
        idx = abs(corr_times - ...
              grid_points{feature_num}(grid_point_num, 1)) < ...
              spacing / 2;
        if (sum(corrs(idx) > threshold) == 0)
            prefilter_cell_idx{feature_num}(grid_point_num) = false;
        end
        grid_point_corrs(counter) = max(corrs(idx));
        counter = counter + 1;
    end
        
    % Keep only blocks whose range has a significant correlation with data        
    %grid_idx = floor((corr_times  + floor(size(data, 1) / 2)) / spacing) + 1;
    %flags = true(size(grid_points{feature_num}, 1), 1);
    %flags(grid_idx(corr_times 
    %corrs > threshold;
    
end
prefilter_idx = cell2mat(prefilter_cell_idx);